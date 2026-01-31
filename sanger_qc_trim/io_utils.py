"""I/O utilities for file discovery, format detection, and sequence parsing."""

import logging
import os
from pathlib import Path
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def discover_files(inputs: List[str], recursive: bool = False) -> List[Tuple[Path, str]]:
    """
    Discover all .ab1 and .phd.1 files from input paths.

    Deduplicates files with the same sample_id, preferring .ab1 over .phd.1.

    Args:
        inputs: List of file or directory paths
        recursive: Whether to recursively search directories

    Returns:
        List of tuples (file_path, format) where format is 'ab1' or 'phd.1'
    """
    files = []
    seen_samples = {}  # Track sample_id -> (path, format) for deduplication

    for input_path in inputs:
        path = Path(input_path)

        if path.is_file():
            fmt = detect_format(path)
            if fmt:
                files.append((path, fmt))
            else:
                logger.warning(f"Skipping file with unknown format: {path}")

        elif path.is_dir():
            if recursive:
                pattern = "**/*"
            else:
                pattern = "*"

            for file_path in path.glob(pattern):
                if file_path.is_file():
                    fmt = detect_format(file_path)
                    if fmt:
                        files.append((file_path, fmt))

        else:
            logger.warning(f"Input path does not exist: {path}")

    # Deduplicate by sample_id (prefer .ab1 over .phd.1)
    for file_path, fmt in files:
        sample_id = get_sample_id(file_path)

        if sample_id not in seen_samples:
            seen_samples[sample_id] = (file_path, fmt)
        else:
            # If we already have this sample_id, prefer .ab1 format
            existing_path, existing_fmt = seen_samples[sample_id]
            if fmt == "ab1" and existing_fmt == "phd.1":
                logger.info(f"Replacing {existing_path.name} with {file_path.name} (preferring .ab1)")
                seen_samples[sample_id] = (file_path, fmt)
            elif fmt == "phd.1" and existing_fmt == "ab1":
                logger.info(f"Skipping {file_path.name} (already have .ab1 version)")
            else:
                logger.warning(f"Duplicate sample_id '{sample_id}': {existing_path.name} and {file_path.name}")

    deduplicated_files = list(seen_samples.values())

    if len(deduplicated_files) < len(files):
        logger.info(f"Deduplicated {len(files)} files to {len(deduplicated_files)} unique samples")

    logger.info(f"Discovered {len(deduplicated_files)} files to process")
    return deduplicated_files


def detect_format(file_path: Path) -> Optional[str]:
    """
    Detect file format based on extension.

    Args:
        file_path: Path to file

    Returns:
        'ab1' for .ab1 files, 'phd.1' for .phd.1 files, None for unknown
    """
    suffix = file_path.suffix.lower()
    name = file_path.name.lower()

    if suffix == ".ab1":
        return "ab1"
    elif name.endswith(".phd.1"):
        return "phd.1"
    else:
        return None


def get_sample_id(file_path: Path) -> str:
    """
    Extract sample ID from filename (filename without extension).

    Args:
        file_path: Path to file

    Returns:
        Sample ID string
    """
    name = file_path.name

    # For .phd.1 files, remove both .phd and .1
    if name.lower().endswith(".phd.1"):
        return name[:-6]
    else:
        return file_path.stem


def parse_sequence_file(
    file_path: Path, file_format: str
) -> Optional[Tuple[str, List[int]]]:
    """
    Parse a sequence file and extract sequence and quality scores.

    Args:
        file_path: Path to sequence file
        file_format: File format ('ab1' or 'phd.1')

    Returns:
        Tuple of (sequence_string, quality_list) or None if parsing fails
    """
    try:
        if file_format == "ab1":
            return _parse_ab1(file_path)
        elif file_format == "phd.1":
            return _parse_phd(file_path)
        else:
            logger.warning(f"Unknown format '{file_format}' for file: {file_path}")
            return None
    except Exception as e:
        logger.warning(f"Failed to parse {file_path}: {e}")
        return None


def _parse_ab1(file_path: Path) -> Optional[Tuple[str, List[int]]]:
    """
    Parse an AB1 file.

    Args:
        file_path: Path to .ab1 file

    Returns:
        Tuple of (sequence, qualities) or None if no quality data
    """
    record = SeqIO.read(str(file_path), "abi")
    seq = str(record.seq)

    # Check for quality scores
    if "phred_quality" not in record.letter_annotations:
        logger.warning(f"No quality scores found in AB1 file: {file_path}")
        return None

    quals = record.letter_annotations["phred_quality"]

    return (seq, quals)


def _parse_phd(file_path: Path) -> Optional[Tuple[str, List[int]]]:
    """
    Parse a PHD.1 file.

    Args:
        file_path: Path to .phd.1 file

    Returns:
        Tuple of (sequence, qualities) or None if parsing fails
    """
    # PHD files usually contain one record, but we'll take the first
    records = list(SeqIO.parse(str(file_path), "phd"))

    if not records:
        logger.warning(f"No sequences found in PHD file: {file_path}")
        return None

    record = records[0]
    seq = str(record.seq)

    # Check for quality scores
    if "phred_quality" not in record.letter_annotations:
        logger.warning(f"No quality scores found in PHD file: {file_path}")
        return None

    quals = record.letter_annotations["phred_quality"]

    return (seq, quals)


def make_read_id(sample_id: str, trim_start: int, trim_end: int) -> str:
    """
    Create a read ID with trimming coordinates.

    Args:
        sample_id: Sample identifier
        trim_start: Trim start position
        trim_end: Trim end position

    Returns:
        Read ID string in format "sample_id/trim:start-end"
    """
    return f"{sample_id}/trim:{trim_start}-{trim_end}"
