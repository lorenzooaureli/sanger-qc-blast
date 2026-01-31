"""Tests for ambiguous base calling module."""

import pytest
import numpy as np
from pathlib import Path
from sanger_qc_trim.ambiguous_calling import (
    AmbiguousCallingConfig,
    AmbiguousBaseCaller,
    BaseCall,
    PeakIntensityExtractor,
    IUPAC_CODES,
    create_ambiguous_caller,
)


class TestIUPACCodes:
    """Test IUPAC ambiguity code mapping."""

    def test_iupac_r_code(self):
        """Test R code for A/G mixture."""
        assert IUPAC_CODES[frozenset(['A', 'G'])] == 'R'

    def test_iupac_y_code(self):
        """Test Y code for C/T mixture."""
        assert IUPAC_CODES[frozenset(['C', 'T'])] == 'Y'

    def test_iupac_s_code(self):
        """Test S code for G/C mixture."""
        assert IUPAC_CODES[frozenset(['G', 'C'])] == 'S'

    def test_iupac_w_code(self):
        """Test W code for A/T mixture."""
        assert IUPAC_CODES[frozenset(['A', 'T'])] == 'W'

    def test_iupac_k_code(self):
        """Test K code for G/T mixture."""
        assert IUPAC_CODES[frozenset(['G', 'T'])] == 'K'

    def test_iupac_m_code(self):
        """Test M code for A/C mixture."""
        assert IUPAC_CODES[frozenset(['A', 'C'])] == 'M'


class TestAmbiguousCallingConfig:
    """Test configuration dataclass."""

    def test_default_config(self):
        """Test default configuration values."""
        config = AmbiguousCallingConfig()
        assert config.q_min_noise == 12
        assert config.snr_min == 4.0
        assert config.spr_noise_max == 0.20
        assert config.spr_het_low == 0.33
        assert config.spr_het_high == 0.67
        assert config.q_confident == 30
        assert config.q_ambig == 20
        assert config.clonal_context is True

    def test_custom_config(self):
        """Test custom configuration values."""
        config = AmbiguousCallingConfig(
            q_min_noise=15,
            snr_min=5.0,
            spr_noise_max=0.15,
            clonal_context=False
        )
        assert config.q_min_noise == 15
        assert config.snr_min == 5.0
        assert config.spr_noise_max == 0.15
        assert config.clonal_context is False


class TestAmbiguousBaseCaller:
    """Test ambiguous base caller logic."""

    def test_get_iupac_code_ag(self):
        """Test IUPAC code retrieval for A/G."""
        caller = AmbiguousBaseCaller()
        assert caller._get_iupac_code('A', 'G') == 'R'
        assert caller._get_iupac_code('G', 'A') == 'R'

    def test_get_iupac_code_ct(self):
        """Test IUPAC code retrieval for C/T."""
        caller = AmbiguousBaseCaller()
        assert caller._get_iupac_code('C', 'T') == 'Y'

    def test_get_iupac_code_invalid(self):
        """Test IUPAC code for invalid combination returns N."""
        caller = AmbiguousBaseCaller()
        assert caller._get_iupac_code('X', 'Y') == 'N'

    def test_rule1_low_quality_returns_n(self):
        """Rule 1: Low quality should return N."""
        caller = AmbiguousBaseCaller()
        # Q=10 < q_min_noise=12
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=20, spr=0.20, snr=5.0, quality=10
        )
        assert called == 'N'
        assert mode == 'N'
        assert 'low_quality' in flags

    def test_rule1_low_snr_returns_n(self):
        """Rule 1: Low SNR should return N."""
        caller = AmbiguousBaseCaller()
        # SNR=3 < snr_min=4.0
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=20, spr=0.20, snr=3.0, quality=25
        )
        assert called == 'N'
        assert mode == 'N'
        assert 'low_quality' in flags

    def test_rule2_clear_single_base(self):
        """Rule 2: Clear single base with high quality and low SPR."""
        caller = AmbiguousBaseCaller()
        # Q=35 >= q_confident=30, SPR=0.10 < spr_noise_max=0.20
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=10, spr=0.10, snr=10.0, quality=35
        )
        assert called == 'A'
        assert mode == 'single'
        assert len(flags) == 0

    def test_rule3_minor_mixture_clonal(self):
        """Rule 3: Minor mixture in clonal context returns primary base."""
        config = AmbiguousCallingConfig(clonal_context=True)
        caller = AmbiguousBaseCaller(config)
        # SPR=0.25 in borderline range (0.20 - 0.33)
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=25, spr=0.25, snr=10.0, quality=25
        )
        assert called == 'A'
        assert mode == 'single'
        assert 'minor_secondary' in flags

    def test_rule3_minor_mixture_mixed(self):
        """Rule 3: Minor mixture in mixed context returns IUPAC."""
        config = AmbiguousCallingConfig(clonal_context=False)
        caller = AmbiguousBaseCaller(config)
        # SPR=0.25 in borderline range (0.20 - 0.33)
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=25, spr=0.25, snr=10.0, quality=25
        )
        assert called == 'R'
        assert mode == 'ambiguous'

    def test_rule4_heterozygous_balanced(self):
        """Rule 4: Balanced heterozygous mixture returns IUPAC."""
        caller = AmbiguousBaseCaller()
        # SPR=0.50 in heterozygous range (0.33 - 0.67)
        called, mode, allele_frac, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=50, spr=0.50, snr=10.0, quality=25
        )
        assert called == 'R'
        assert mode == 'ambiguous'
        assert 'heterozygous' in flags
        # allele_frac should be H1/(H1+H2) = 100/150 = 0.6667
        assert 0.66 < allele_frac < 0.68

    def test_rule5_unbalanced_mixture_mixed(self):
        """Rule 5: Unbalanced mixture in mixed context returns IUPAC."""
        config = AmbiguousCallingConfig(clonal_context=False)
        caller = AmbiguousBaseCaller(config)
        # SPR=0.75 in unbalanced range (0.67 - 0.85)
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=75, spr=0.75, snr=10.0, quality=25
        )
        assert called == 'R'
        assert mode == 'ambiguous'
        assert 'unbalanced_mixture' in flags

    def test_rule5_unbalanced_mixture_clonal(self):
        """Rule 5: Unbalanced mixture in clonal context returns primary base."""
        config = AmbiguousCallingConfig(clonal_context=True)
        caller = AmbiguousBaseCaller(config)
        # SPR=0.75 in unbalanced range (0.67 - 0.85)
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=75, spr=0.75, snr=10.0, quality=25
        )
        assert called == 'A'
        assert mode == 'single'
        assert 'uncertain_mixture' in flags

    def test_rule6_very_high_spr_returns_n(self):
        """Rule 6: Very high SPR (>0.85) returns N."""
        caller = AmbiguousBaseCaller()
        # SPR=0.90 > spr_unbalanced=0.85
        called, mode, _, flags = caller._apply_calling_criteria(
            b1='A', b2='G', H1=100, H2=90, spr=0.90, snr=10.0, quality=25
        )
        assert called == 'N'
        assert mode == 'N'
        assert 'unclear_primary' in flags

    def test_base_calls_to_sequence(self):
        """Test converting base calls to sequence string."""
        caller = AmbiguousBaseCaller()
        base_calls = [
            BaseCall(0, 'A', 'A', 'G', 100, 10, 0.10, 10, 35, 'single', 0.91, []),
            BaseCall(1, 'R', 'A', 'G', 100, 50, 0.50, 10, 25, 'ambiguous', 0.67, ['heterozygous']),
            BaseCall(2, 'C', 'C', 'T', 100, 15, 0.15, 10, 35, 'single', 0.87, []),
            BaseCall(3, 'N', 'A', 'G', 50, 45, 0.90, 3, 10, 'N', 0.53, ['unclear_primary']),
        ]
        seq = caller.base_calls_to_sequence(base_calls)
        assert seq == 'ARCN'

    def test_base_calls_to_annotations(self):
        """Test converting base calls to annotation dictionaries."""
        caller = AmbiguousBaseCaller()
        base_calls = [
            BaseCall(0, 'A', 'A', 'G', 100.5, 10.2, 0.10, 10.5, 35, 'single', 0.91, []),
            BaseCall(1, 'R', 'A', 'G', 100.0, 50.0, 0.50, 10.0, 25, 'ambiguous', 0.67, ['heterozygous']),
        ]
        annotations = caller.base_calls_to_annotations(base_calls)

        assert len(annotations) == 2
        assert annotations[0]['position'] == 0
        assert annotations[0]['called_base'] == 'A'
        assert annotations[0]['primary_base'] == 'A'
        assert annotations[0]['secondary_base'] == 'G'
        assert annotations[0]['spr'] == 0.1000
        assert annotations[0]['call_mode'] == 'single'
        assert annotations[0]['flags'] == ''

        assert annotations[1]['position'] == 1
        assert annotations[1]['called_base'] == 'R'
        assert annotations[1]['call_mode'] == 'ambiguous'
        assert annotations[1]['flags'] == 'heterozygous'

    def test_fallback_calling_high_quality(self):
        """Test fallback calling with high quality bases."""
        caller = AmbiguousBaseCaller()
        sequence = "ACGT"
        qualities = [35, 35, 35, 35]

        base_calls = caller._fallback_calling(sequence, qualities)

        assert len(base_calls) == 4
        for i, bc in enumerate(base_calls):
            assert bc.called_base == sequence[i]
            assert bc.call_mode == 'single'
            assert 'no_trace_data' in bc.flags
            assert bc.quality == 35

    def test_fallback_calling_low_quality(self):
        """Test fallback calling with low quality bases."""
        caller = AmbiguousBaseCaller()
        sequence = "ACGT"
        qualities = [10, 10, 10, 10]

        base_calls = caller._fallback_calling(sequence, qualities)

        assert len(base_calls) == 4
        for bc in base_calls:
            assert bc.called_base == 'N'
            assert bc.call_mode == 'N'
            assert 'low_quality' in bc.flags
            assert 'no_trace_data' in bc.flags

    def test_fallback_calling_moderate_quality(self):
        """Test fallback calling with moderate quality bases."""
        caller = AmbiguousBaseCaller()
        sequence = "ACGT"
        qualities = [20, 20, 20, 20]

        base_calls = caller._fallback_calling(sequence, qualities)

        assert len(base_calls) == 4
        for i, bc in enumerate(base_calls):
            assert bc.called_base == sequence[i]
            assert bc.call_mode == 'single'
            assert 'moderate_quality' in bc.flags
            assert 'no_trace_data' in bc.flags


class TestPeakIntensityExtractor:
    """Test peak intensity extraction (limited without real AB1 files)."""

    def test_get_intensities_no_peak_locations(self):
        """Test intensity extraction with missing peak location data."""
        extractor = PeakIntensityExtractor()

        # Mock trace data without peak locations
        traces = {
            'A': np.array([10, 20, 30, 20, 10]),
            'C': np.array([5, 10, 15, 10, 5]),
            'G': np.array([8, 16, 24, 16, 8]),
            'T': np.array([3, 6, 9, 6, 3]),
            'peak_locations': np.array([])
        }

        intensities = extractor.get_intensities_at_position(traces, 0, window=1)

        # Should return zeros when no peak locations
        assert intensities == {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}

    def test_get_intensities_with_peak_locations(self):
        """Test intensity extraction with peak location data."""
        extractor = PeakIntensityExtractor()

        # Mock trace data with peak at position 2
        traces = {
            'A': np.array([10, 20, 100, 20, 10]),
            'C': np.array([5, 10, 50, 10, 5]),
            'G': np.array([8, 16, 80, 16, 8]),
            'T': np.array([3, 6, 30, 6, 3]),
            'peak_locations': np.array([2, 10, 20])  # First base peak at trace position 2
        }

        intensities = extractor.get_intensities_at_position(traces, 0, window=1)

        # Should integrate around position 2 with window=1 (positions 1,2,3)
        assert intensities['A'] == 20 + 100 + 20  # 140
        assert intensities['C'] == 10 + 50 + 10   # 70
        assert intensities['G'] == 16 + 80 + 16   # 112
        assert intensities['T'] == 6 + 30 + 6     # 42


class TestCreateAmbiguousCaller:
    """Test factory function for creating callers."""

    def test_create_caller_default(self):
        """Test creating caller with default settings."""
        caller = create_ambiguous_caller()
        assert isinstance(caller, AmbiguousBaseCaller)
        assert caller.config.clonal_context is True
        assert caller.config.q_min_noise == 12
        assert caller.config.snr_min == 4.0

    def test_create_caller_custom(self):
        """Test creating caller with custom settings."""
        caller = create_ambiguous_caller(
            clonal_context=False,
            q_min_noise=15,
            snr_min=5.0,
            spr_noise_max=0.15,
            spr_het_low=0.30,
            spr_het_high=0.70
        )
        assert isinstance(caller, AmbiguousBaseCaller)
        assert caller.config.clonal_context is False
        assert caller.config.q_min_noise == 15
        assert caller.config.snr_min == 5.0
        assert caller.config.spr_noise_max == 0.15
        assert caller.config.spr_het_low == 0.30
        assert caller.config.spr_het_high == 0.70


class TestBaseCallDataclass:
    """Test BaseCall dataclass."""

    def test_base_call_creation(self):
        """Test creating a BaseCall object."""
        bc = BaseCall(
            position=10,
            called_base='R',
            primary_base='A',
            secondary_base='G',
            primary_intensity=100.0,
            secondary_intensity=50.0,
            spr=0.50,
            snr=10.0,
            quality=25,
            call_mode='ambiguous',
            allele_fraction=0.67,
            flags=['heterozygous']
        )

        assert bc.position == 10
        assert bc.called_base == 'R'
        assert bc.primary_base == 'A'
        assert bc.secondary_base == 'G'
        assert bc.spr == 0.50
        assert bc.snr == 10.0
        assert bc.call_mode == 'ambiguous'
        assert 'heterozygous' in bc.flags
