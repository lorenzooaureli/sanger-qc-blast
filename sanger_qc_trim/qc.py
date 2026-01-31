"""QC statistics computation for Sanger sequencing reads."""

import numpy as np
from typing import Dict, Any


def compute_qc_metrics(
    sample_id: str,
    source_file: str,
    file_format: str,
    seq: str,
    quals: list[int],
    trim_start: int,
    trim_end: int,
    qthreshold: int,
    min_length: int,
) -> Dict[str, Any]:
    """
    Compute comprehensive QC statistics for a single read.

    Args:
        sample_id: Sample identifier (from filename without extension)
        source_file: Path to source file
        file_format: File format (ab1 or phd.1)
        seq: DNA sequence string
        quals: List of phred quality scores
        trim_start: Trimming start position (0-based)
        trim_end: Trimming end position (0-based, exclusive)
        qthreshold: Quality threshold used for trimming
        min_length: Minimum acceptable trimmed length

    Returns:
        Dictionary of QC metrics
    """
    q_array = np.array(quals, dtype=float)
    raw_length = len(seq)
    trimmed_length = trim_end - trim_start

    # Basic quality statistics
    mean_q = float(np.mean(q_array)) if len(q_array) > 0 else 0.0
    median_q = float(np.median(q_array)) if len(q_array) > 0 else 0.0

    # Quality thresholds
    pct_q20 = float(np.mean(q_array >= 20)) if len(q_array) > 0 else 0.0
    pct_q30 = float(np.mean(q_array >= 30)) if len(q_array) > 0 else 0.0

    # GC content (ignoring N)
    seq_upper = seq.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    n_count = seq_upper.count("N")
    non_n_count = raw_length - n_count
    gc_percent = (100.0 * gc_count / non_n_count) if non_n_count > 0 else 0.0

    # Expected errors
    expected_errors = float(np.sum(10 ** (-q_array / 10))) if len(q_array) > 0 else 0.0

    # Longest high-quality stretch
    hq_longest_stretch_len = _longest_hq_stretch(quals, qthreshold)

    # Pass/fail
    passed_minlen = "yes" if trimmed_length >= min_length else "no"

    return {
        "sample_id": sample_id,
        "source_file": source_file,
        "format": file_format,
        "raw_length": raw_length,
        "mean_q": round(mean_q, 2),
        "median_q": round(median_q, 2),
        "pct_q20": round(pct_q20, 4),
        "pct_q30": round(pct_q30, 4),
        "gc_percent": round(gc_percent, 2),
        "n_count": n_count,
        "expected_errors": round(expected_errors, 2),
        "hq_longest_stretch_len": hq_longest_stretch_len,
        "trim_start": trim_start,
        "trim_end": trim_end,
        "trimmed_length": trimmed_length,
        "passed_minlen": passed_minlen,
    }


def _longest_hq_stretch(quals: list[int], threshold: int) -> int:
    """
    Find the longest contiguous stretch of bases with Q >= threshold.

    Args:
        quals: List of phred quality scores
        threshold: Quality threshold

    Returns:
        Length of longest high-quality stretch
    """
    if not quals:
        return 0

    max_len = 0
    current_len = 0

    for q in quals:
        if q >= threshold:
            current_len += 1
            max_len = max(max_len, current_len)
        else:
            current_len = 0

    return max_len


def compute_summary_stats(metrics_list: list[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Compute aggregated summary statistics across all reads.

    Args:
        metrics_list: List of per-read metric dictionaries

    Returns:
        Dictionary of summary statistics
    """
    if not metrics_list:
        return {
            "total_reads": 0,
            "mean_raw_length": 0.0,
            "median_raw_length": 0.0,
            "mean_trimmed_length": 0.0,
            "median_trimmed_length": 0.0,
            "mean_mean_q": 0.0,
            "median_mean_q": 0.0,
            "mean_pct_q20": 0.0,
            "mean_pct_q30": 0.0,
            "mean_gc_percent": 0.0,
            "mean_expected_errors": 0.0,
            "reads_passed_minlen": 0,
            "reads_failed_minlen": 0,
            "pct_passed": 0.0,
        }

    # Extract arrays
    raw_lengths = np.array([m["raw_length"] for m in metrics_list])
    trimmed_lengths = np.array([m["trimmed_length"] for m in metrics_list])
    mean_qs = np.array([m["mean_q"] for m in metrics_list])
    pct_q20s = np.array([m["pct_q20"] for m in metrics_list])
    pct_q30s = np.array([m["pct_q30"] for m in metrics_list])
    gc_percents = np.array([m["gc_percent"] for m in metrics_list])
    expected_errors_list = np.array([m["expected_errors"] for m in metrics_list])
    passed = [m["passed_minlen"] == "yes" for m in metrics_list]

    reads_passed = sum(passed)
    reads_failed = len(metrics_list) - reads_passed
    pct_passed = 100.0 * reads_passed / len(metrics_list)

    return {
        "total_reads": len(metrics_list),
        "mean_raw_length": round(float(np.mean(raw_lengths)), 2),
        "median_raw_length": round(float(np.median(raw_lengths)), 2),
        "mean_trimmed_length": round(float(np.mean(trimmed_lengths)), 2),
        "median_trimmed_length": round(float(np.median(trimmed_lengths)), 2),
        "mean_mean_q": round(float(np.mean(mean_qs)), 2),
        "median_mean_q": round(float(np.median(mean_qs)), 2),
        "mean_pct_q20": round(float(np.mean(pct_q20s)), 4),
        "mean_pct_q30": round(float(np.mean(pct_q30s)), 4),
        "mean_gc_percent": round(float(np.mean(gc_percents)), 2),
        "mean_expected_errors": round(float(np.mean(expected_errors_list)), 2),
        "reads_passed_minlen": reads_passed,
        "reads_failed_minlen": reads_failed,
        "pct_passed": round(pct_passed, 2),
    }
