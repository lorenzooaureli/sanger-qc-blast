"""Ambiguous base calling for Sanger chromatograms with IUPAC support.

This module implements advanced base calling that can detect and report
mixed signals (heterozygous positions, contamination, etc.) using peak
intensity analysis from AB1 chromatogram files.
"""

import logging
import numpy as np
from typing import Tuple, List, Dict, Any, Optional
from pathlib import Path
from Bio import SeqIO
from dataclasses import dataclass

logger = logging.getLogger(__name__)


# IUPAC ambiguity codes for two-base mixtures
IUPAC_CODES = {
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
}

# Reverse mapping for IUPAC codes
IUPAC_REVERSE = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'N': ['A', 'C', 'G', 'T'],
}


@dataclass
class AmbiguousCallingConfig:
    """Configuration parameters for ambiguous base calling."""

    # Quality and signal thresholds
    q_min_noise: int = 12          # Minimum quality for non-N calls
    snr_min: float = 4.0           # Minimum signal-to-noise ratio
    spr_noise_max: float = 0.20    # Max SPR to consider secondary as noise
    spr_het_low: float = 0.33      # Lower bound for heterozygous mixture
    spr_het_high: float = 0.67     # Upper bound for heterozygous mixture
    spr_unbalanced: float = 0.85   # Threshold for unbalanced mixtures
    q_confident: int = 30          # Quality for confident single-base calls
    q_ambig: int = 20              # Minimum quality for ambiguous calls

    # Context-specific behavior
    clonal_context: bool = True    # True for clonal/plasmid, False for mixed
    peak_spacing_tolerance: float = 0.20  # Max deviation for peak spacing

    # Intensity measurement parameters
    peak_window: int = 3           # Window size for peak integration (±bases)


@dataclass
class BaseCall:
    """Result of ambiguous base calling at a single position."""

    position: int                  # 0-based position in sequence
    called_base: str               # Final base call (A/C/G/T/IUPAC/N)
    primary_base: str              # Primary base by intensity
    secondary_base: str            # Secondary base by intensity
    primary_intensity: float       # H1
    secondary_intensity: float     # H2
    spr: float                     # Secondary-to-Primary Ratio
    snr: float                     # Signal-to-Noise Ratio
    quality: int                   # Phred quality score
    call_mode: str                 # 'single', 'ambiguous', 'N'
    allele_fraction: float         # H1 / (H1 + H2) for heterozygous
    flags: List[str]               # Additional flags/warnings


class PeakIntensityExtractor:
    """Extract peak intensities from AB1 chromatogram files."""

    @staticmethod
    def extract_traces(ab1_path: Path) -> Optional[Dict[str, np.ndarray]]:
        """
        Extract raw trace data from AB1 file.

        Args:
            ab1_path: Path to .ab1 file

        Returns:
            Dictionary with keys 'A', 'C', 'G', 'T' mapping to intensity arrays,
            or None if extraction fails
        """
        try:
            record = SeqIO.read(str(ab1_path), "abi")

            # Access raw trace data from AB1 annotations
            # Channel mapping: DATA9='A', DATA10='C', DATA11='G', DATA12='T'
            traces = {}

            if hasattr(record, 'annotations') and 'abif_raw' in record.annotations:
                abif_raw = record.annotations['abif_raw']

                # Extract trace data for each channel
                traces['A'] = np.array(abif_raw.get('DATA9', []), dtype=float)
                traces['C'] = np.array(abif_raw.get('DATA10', []), dtype=float)
                traces['G'] = np.array(abif_raw.get('DATA11', []), dtype=float)
                traces['T'] = np.array(abif_raw.get('DATA12', []), dtype=float)

                # Get base positions (peak locations in trace)
                peak_locations = abif_raw.get('PLOC2', [])
                if not peak_locations:
                    peak_locations = abif_raw.get('PLOC1', [])

                traces['peak_locations'] = np.array(peak_locations, dtype=int)

                # Validate we got data
                if all(len(traces[base]) > 0 for base in ['A', 'C', 'G', 'T']):
                    return traces

            logger.warning(f"Could not extract trace data from {ab1_path}")
            return None

        except Exception as e:
            logger.warning(f"Failed to extract traces from {ab1_path}: {e}")
            return None

    @staticmethod
    def get_intensities_at_position(
        traces: Dict[str, np.ndarray],
        position: int,
        window: int = 3
    ) -> Dict[str, float]:
        """
        Get peak intensities for all four bases at a specific position.

        Args:
            traces: Dictionary of trace arrays for A/C/G/T
            position: Base position (0-based)
            window: Integration window size (±bases)

        Returns:
            Dictionary mapping base -> integrated peak intensity
        """
        peak_locations = traces.get('peak_locations', [])

        if len(peak_locations) == 0 or position >= len(peak_locations):
            # No peak location data, return zeros
            return {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}

        # Get trace position for this base call
        trace_pos = peak_locations[position]

        # Integrate intensities in a window around the peak
        intensities = {}
        for base in ['A', 'C', 'G', 'T']:
            trace = traces[base]

            # Define window bounds
            start = max(0, trace_pos - window)
            end = min(len(trace), trace_pos + window + 1)

            # Sum intensities in window (area under peak)
            intensities[base] = float(np.sum(trace[start:end]))

        return intensities


class AmbiguousBaseCaller:
    """Perform ambiguous base calling using peak intensity analysis."""

    def __init__(self, config: Optional[AmbiguousCallingConfig] = None):
        """
        Initialize ambiguous base caller.

        Args:
            config: Configuration parameters (uses defaults if None)
        """
        self.config = config if config else AmbiguousCallingConfig()
        self.extractor = PeakIntensityExtractor()

    def call_bases(
        self,
        ab1_path: Path,
        sequence: str,
        qualities: List[int]
    ) -> List[BaseCall]:
        """
        Perform ambiguous base calling on an AB1 file.

        Args:
            ab1_path: Path to .ab1 file
            sequence: Called sequence string
            qualities: List of Phred quality scores

        Returns:
            List of BaseCall objects, one per position
        """
        # Extract trace data
        traces = self.extractor.extract_traces(ab1_path)

        if traces is None:
            # Fall back to simple quality-based calling
            logger.warning(f"Using fallback calling for {ab1_path} (no trace data)")
            return self._fallback_calling(sequence, qualities)

        # Process each position
        base_calls = []
        for pos in range(len(sequence)):
            base_call = self._call_position(
                position=pos,
                original_base=sequence[pos],
                quality=qualities[pos] if pos < len(qualities) else 0,
                traces=traces
            )
            base_calls.append(base_call)

        return base_calls

    def _call_position(
        self,
        position: int,
        original_base: str,
        quality: int,
        traces: Dict[str, np.ndarray]
    ) -> BaseCall:
        """
        Call a single base position using the 6-rule criteria.

        Args:
            position: 0-based position
            original_base: Original base call from sequencer
            quality: Phred quality score
            traces: Trace data dictionary

        Returns:
            BaseCall object with full annotation
        """
        # Get intensities at this position
        intensities = self.extractor.get_intensities_at_position(
            traces, position, self.config.peak_window
        )

        # Sort bases by intensity
        sorted_bases = sorted(
            intensities.items(),
            key=lambda x: x[1],
            reverse=True
        )

        b1, H1 = sorted_bases[0]  # Primary base
        b2, H2 = sorted_bases[1]  # Secondary base

        # Calculate metrics
        spr = H2 / H1 if H1 > 0 else 0.0
        snr = self._calculate_snr(traces, position, H1)

        # Apply the 6 default calling criteria
        called_base, call_mode, allele_frac, flags = self._apply_calling_criteria(
            b1, b2, H1, H2, spr, snr, quality
        )

        return BaseCall(
            position=position,
            called_base=called_base,
            primary_base=b1,
            secondary_base=b2,
            primary_intensity=H1,
            secondary_intensity=H2,
            spr=spr,
            snr=snr,
            quality=quality,
            call_mode=call_mode,
            allele_fraction=allele_frac,
            flags=flags
        )

    def _calculate_snr(
        self,
        traces: Dict[str, np.ndarray],
        position: int,
        signal: float
    ) -> float:
        """
        Calculate signal-to-noise ratio at a position.

        Args:
            traces: Trace data dictionary
            position: Base position
            signal: Primary signal intensity (H1)

        Returns:
            SNR value
        """
        peak_locations = traces.get('peak_locations', [])

        if len(peak_locations) == 0 or position >= len(peak_locations):
            return 0.0

        trace_pos = peak_locations[position]

        # Estimate noise from local baseline
        # Use median of minimum values in surrounding regions
        noise_samples = []
        for base in ['A', 'C', 'G', 'T']:
            trace = traces[base]

            # Sample regions before and after peak
            before_start = max(0, trace_pos - 20)
            before_end = max(0, trace_pos - 5)
            after_start = min(len(trace), trace_pos + 5)
            after_end = min(len(trace), trace_pos + 20)

            if before_end > before_start:
                noise_samples.extend(trace[before_start:before_end])
            if after_end > after_start:
                noise_samples.extend(trace[after_start:after_end])

        if noise_samples:
            noise = float(np.median(noise_samples))
            return signal / noise if noise > 0 else 0.0

        return 0.0

    def _apply_calling_criteria(
        self,
        b1: str,
        b2: str,
        H1: float,
        H2: float,
        spr: float,
        snr: float,
        quality: int
    ) -> Tuple[str, str, float, List[str]]:
        """
        Apply the 6 default calling criteria to determine base call.

        Args:
            b1: Primary base
            b2: Secondary base
            H1: Primary intensity
            H2: Secondary intensity
            spr: Secondary-to-Primary Ratio
            snr: Signal-to-Noise Ratio
            quality: Phred quality score

        Returns:
            Tuple of (called_base, call_mode, allele_fraction, flags)
        """
        flags = []

        # Calculate allele fraction
        total_signal = H1 + H2
        allele_frac = H1 / total_signal if total_signal > 0 else 0.0

        # Rule 1: Too noisy -> call N
        if quality < self.config.q_min_noise or snr < self.config.snr_min:
            flags.append('low_quality')
            return 'N', 'N', allele_frac, flags

        # Rule 2: Clear single base -> call b1
        if quality >= self.config.q_confident and spr < self.config.spr_noise_max:
            return b1, 'single', allele_frac, flags

        # Rule 3: Minor mixture (borderline) -> context dependent
        if (self.config.spr_noise_max <= spr < self.config.spr_het_low and
            quality >= self.config.q_ambig):

            if self.config.clonal_context:
                # Clonal sample: call primary with flag
                flags.append('minor_secondary')
                return b1, 'single', allele_frac, flags
            else:
                # Mixed context: call IUPAC
                iupac = self._get_iupac_code(b1, b2)
                return iupac, 'ambiguous', allele_frac, flags

        # Rule 4: Heterozygous or mixed signal -> call IUPAC
        if (self.config.spr_het_low <= spr <= self.config.spr_het_high and
            quality >= self.config.q_ambig and
            snr >= self.config.snr_min):

            iupac = self._get_iupac_code(b1, b2)
            flags.append('heterozygous')
            return iupac, 'ambiguous', allele_frac, flags

        # Rule 5: Dominant + real minor (unbalanced het) -> IUPAC or b1 with flag
        if (self.config.spr_het_high < spr < self.config.spr_unbalanced and
            quality >= self.config.q_ambig):

            if not self.config.clonal_context:
                # Mixed samples: prefer IUPAC
                iupac = self._get_iupac_code(b1, b2)
                flags.append('unbalanced_mixture')
                return iupac, 'ambiguous', allele_frac, flags
            else:
                # Clonal samples: call b1 with flag
                flags.append('uncertain_mixture')
                return b1, 'single', allele_frac, flags

        # Rule 6: Default fallback
        # If SPR is very high (>0.85) or other edge cases
        if spr >= self.config.spr_unbalanced:
            # Secondary almost as strong as primary - unclear which to call
            flags.append('unclear_primary')
            return 'N', 'N', allele_frac, flags

        # Final fallback: return primary base
        return b1, 'single', allele_frac, flags

    @staticmethod
    def _get_iupac_code(base1: str, base2: str) -> str:
        """
        Get IUPAC ambiguity code for two bases.

        Args:
            base1: First base
            base2: Second base

        Returns:
            IUPAC ambiguity code or 'N' if not found
        """
        base_set = frozenset([base1.upper(), base2.upper()])
        return IUPAC_CODES.get(base_set, 'N')

    def _fallback_calling(
        self,
        sequence: str,
        qualities: List[int]
    ) -> List[BaseCall]:
        """
        Fallback to quality-based calling when trace data unavailable.

        Args:
            sequence: Called sequence string
            qualities: List of Phred quality scores

        Returns:
            List of BaseCall objects with limited information
        """
        base_calls = []

        for pos, base in enumerate(sequence):
            qual = qualities[pos] if pos < len(qualities) else 0

            # Simple quality-based decision
            if qual < self.config.q_min_noise:
                called = 'N'
                mode = 'N'
                flags = ['low_quality', 'no_trace_data']
            elif qual >= self.config.q_confident:
                called = base
                mode = 'single'
                flags = ['no_trace_data']
            else:
                called = base
                mode = 'single'
                flags = ['moderate_quality', 'no_trace_data']

            base_calls.append(BaseCall(
                position=pos,
                called_base=called,
                primary_base=base,
                secondary_base='N',
                primary_intensity=0.0,
                secondary_intensity=0.0,
                spr=0.0,
                snr=0.0,
                quality=qual,
                call_mode=mode,
                allele_fraction=1.0,
                flags=flags
            ))

        return base_calls

    @staticmethod
    def base_calls_to_sequence(base_calls: List[BaseCall]) -> str:
        """
        Convert list of BaseCall objects to sequence string.

        Args:
            base_calls: List of BaseCall objects

        Returns:
            Sequence string with IUPAC codes
        """
        return ''.join(bc.called_base for bc in base_calls)

    @staticmethod
    def base_calls_to_annotations(base_calls: List[BaseCall]) -> List[Dict[str, Any]]:
        """
        Convert BaseCall objects to annotation dictionaries for output.

        Args:
            base_calls: List of BaseCall objects

        Returns:
            List of annotation dictionaries
        """
        annotations = []

        for bc in base_calls:
            annotations.append({
                'position': bc.position,
                'called_base': bc.called_base,
                'primary_base': bc.primary_base,
                'secondary_base': bc.secondary_base,
                'primary_intensity': round(bc.primary_intensity, 2),
                'secondary_intensity': round(bc.secondary_intensity, 2),
                'spr': round(bc.spr, 4),
                'snr': round(bc.snr, 2),
                'quality': bc.quality,
                'call_mode': bc.call_mode,
                'allele_fraction': round(bc.allele_fraction, 4),
                'flags': ','.join(bc.flags) if bc.flags else ''
            })

        return annotations


def create_ambiguous_caller(
    clonal_context: bool = True,
    q_min_noise: int = 12,
    snr_min: float = 4.0,
    spr_noise_max: float = 0.20,
    spr_het_low: float = 0.33,
    spr_het_high: float = 0.67,
    **kwargs
) -> AmbiguousBaseCaller:
    """
    Factory function to create an AmbiguousBaseCaller with custom config.

    Args:
        clonal_context: True for clonal/plasmid samples, False for mixed
        q_min_noise: Minimum quality for non-N calls
        snr_min: Minimum signal-to-noise ratio
        spr_noise_max: Max SPR to consider secondary as noise
        spr_het_low: Lower bound for heterozygous mixture
        spr_het_high: Upper bound for heterozygous mixture
        **kwargs: Additional configuration parameters

    Returns:
        Configured AmbiguousBaseCaller instance
    """
    config = AmbiguousCallingConfig(
        clonal_context=clonal_context,
        q_min_noise=q_min_noise,
        snr_min=snr_min,
        spr_noise_max=spr_noise_max,
        spr_het_low=spr_het_low,
        spr_het_high=spr_het_high,
        **kwargs
    )

    return AmbiguousBaseCaller(config)
