"""Trimming algorithms for Sanger sequencing reads."""

from typing import Tuple


def trim_mott(quals: list[int], threshold: int) -> Tuple[int, int]:
    """
    Modified Mott/Kadane trimming algorithm.

    Finds the maximum-score contiguous subsequence where weights are q - threshold.
    Returns 0-based [start, end) indices. Returns (0, 0) if all qualities are below threshold.

    Args:
        quals: List of phred quality scores
        threshold: Quality threshold (T)

    Returns:
        Tuple of (trim_start, trim_end) as 0-based [start, end) indices

    Example:
        >>> trim_mott([10, 10, 30, 30, 30, 10], 20)
        (2, 5)
    """
    if not quals:
        return (0, 0)

    # Compute weights: positive if good, negative if bad
    weights = [q - threshold for q in quals]

    max_sum = 0
    cur_sum = 0
    best_start = 0
    best_end = 0
    cur_start = 0

    for i, val in enumerate(weights):
        cur_sum += val
        if cur_sum <= 0:
            # Reset the window
            cur_sum = 0
            cur_start = i + 1
        elif cur_sum > max_sum:
            # Found a better window
            max_sum = cur_sum
            best_start = cur_start
            best_end = i + 1

    return (best_start, best_end)


def trim_ends(quals: list[int], threshold: int) -> Tuple[int, int]:
    """
    Hard clip from 5' and 3' ends to first/last base with Q >= threshold.

    Returns 0-based [start, end) indices. Returns (0, 0) if no bases meet threshold.

    Args:
        quals: List of phred quality scores
        threshold: Quality threshold (T)

    Returns:
        Tuple of (trim_start, trim_end) as 0-based [start, end) indices

    Example:
        >>> trim_ends([15, 25, 25, 15], 20)
        (1, 3)
    """
    if not quals:
        return (0, 0)

    n = len(quals)

    # Find first base with Q >= threshold
    start = next((i for i, q in enumerate(quals) if q >= threshold), n)

    # Find last base with Q >= threshold
    end = next((i for i in range(n - 1, -1, -1) if quals[i] >= threshold), -1)

    # Check if valid range found
    if start >= n or end < 0 or end < start:
        return (0, 0)

    return (start, end + 1)  # Convert to [start, end) format


def apply_trim(seq: str, quals: list[int], method: str, threshold: int) -> Tuple[str, list[int], int, int]:
    """
    Apply trimming to a sequence and return trimmed results.

    Args:
        seq: DNA sequence string
        quals: List of phred quality scores
        method: Trimming method ('mott' or 'ends')
        threshold: Quality threshold

    Returns:
        Tuple of (trimmed_seq, trimmed_quals, trim_start, trim_end)
    """
    if method == "mott":
        trim_start, trim_end = trim_mott(quals, threshold)
    elif method == "ends":
        trim_start, trim_end = trim_ends(quals, threshold)
    else:
        raise ValueError(f"Unknown trimming method: {method}")

    trimmed_seq = seq[trim_start:trim_end]
    trimmed_quals = quals[trim_start:trim_end]

    return trimmed_seq, trimmed_quals, trim_start, trim_end
