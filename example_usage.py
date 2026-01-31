#!/usr/bin/env python3
"""
Example usage of the Sanger QC & Trim tool programmatically.

This demonstrates how to use the tool's modules directly in Python code,
rather than through the CLI.
"""

from pathlib import Path
from sanger_qc_trim.io_utils import parse_sequence_file, get_sample_id
from sanger_qc_trim.trim import apply_trim
from sanger_qc_trim.qc import compute_qc_metrics, compute_summary_stats


def example_process_single_file():
    """Example: Process a single file programmatically."""

    print("=" * 60)
    print("Example 1: Processing a single file")
    print("=" * 60)

    # Path to a sample file
    file_path = Path("seq_data/Y19/Y19_ITS4.phd.1")

    if not file_path.exists():
        print(f"File not found: {file_path}")
        return

    # Parse the file
    result = parse_sequence_file(file_path, "phd.1")

    if result is None:
        print("Failed to parse file")
        return

    seq, quals = result
    sample_id = get_sample_id(file_path)

    print(f"Sample ID: {sample_id}")
    print(f"Sequence length: {len(seq)}")
    print(f"Quality scores: {len(quals)}")
    print(f"Mean quality: {sum(quals) / len(quals):.2f}")
    print()

    # Apply Mott trimming at Q20
    trimmed_seq, trimmed_quals, trim_start, trim_end = apply_trim(
        seq, quals, method="mott", threshold=20
    )

    print(f"Trim coordinates: [{trim_start}, {trim_end})")
    print(f"Trimmed length: {len(trimmed_seq)}")
    print(f"Trimmed sequence: {trimmed_seq[:50]}..." if len(trimmed_seq) > 50 else f"Trimmed sequence: {trimmed_seq}")
    print()

    # Compute QC metrics
    metrics = compute_qc_metrics(
        sample_id=sample_id,
        source_file=str(file_path),
        file_format="phd.1",
        seq=seq,
        quals=quals,
        trim_start=trim_start,
        trim_end=trim_end,
        qthreshold=20,
        min_length=50,
    )

    print("QC Metrics:")
    for key, value in metrics.items():
        print(f"  {key}: {value}")
    print()


def example_compare_trimming_methods():
    """Example: Compare Mott vs Ends trimming."""

    print("=" * 60)
    print("Example 2: Comparing trimming methods")
    print("=" * 60)

    # Synthetic example
    seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    quals = [10, 10, 15, 25, 30, 35, 30, 25, 20, 15, 10, 10, 15, 20, 25, 30, 35, 30, 25, 20, 15, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

    print(f"Original length: {len(seq)}")
    print(f"Quality scores: {quals}")
    print()

    # Mott trimming
    mott_seq, mott_quals, mott_start, mott_end = apply_trim(
        seq, quals, method="mott", threshold=20
    )

    print("Mott trimming (Q20):")
    print(f"  Coordinates: [{mott_start}, {mott_end})")
    print(f"  Length: {len(mott_seq)}")
    print(f"  Sequence: {mott_seq}")
    print()

    # Ends trimming
    ends_seq, ends_quals, ends_start, ends_end = apply_trim(
        seq, quals, method="ends", threshold=20
    )

    print("Ends trimming (Q20):")
    print(f"  Coordinates: [{ends_start}, {ends_end})")
    print(f"  Length: {len(ends_seq)}")
    print(f"  Sequence: {ends_seq}")
    print()


def example_batch_processing():
    """Example: Batch process multiple files."""

    print("=" * 60)
    print("Example 3: Batch processing")
    print("=" * 60)

    from sanger_qc_trim.io_utils import discover_files

    # Discover files
    files = discover_files(["seq_data/Y19/"], recursive=False)

    print(f"Found {len(files)} files")
    print()

    metrics_list = []

    for file_path, file_format in files:
        print(f"Processing: {file_path.name}")

        result = parse_sequence_file(file_path, file_format)
        if result is None:
            print(f"  Skipped (parsing failed)")
            continue

        seq, quals = result
        sample_id = get_sample_id(file_path)

        # Apply trimming
        _, _, trim_start, trim_end = apply_trim(seq, quals, "mott", 20)

        # Compute metrics
        metrics = compute_qc_metrics(
            sample_id=sample_id,
            source_file=str(file_path),
            file_format=file_format,
            seq=seq,
            quals=quals,
            trim_start=trim_start,
            trim_end=trim_end,
            qthreshold=20,
            min_length=50,
        )

        metrics_list.append(metrics)
        print(f"  Raw length: {metrics['raw_length']}")
        print(f"  Trimmed length: {metrics['trimmed_length']}")
        print(f"  Mean quality: {metrics['mean_q']}")
        print()

    # Compute summary
    if metrics_list:
        summary = compute_summary_stats(metrics_list)
        print("Summary Statistics:")
        for key, value in summary.items():
            print(f"  {key}: {value}")


if __name__ == "__main__":
    # Run examples
    example_process_single_file()
    print("\n")
    example_compare_trimming_methods()
    print("\n")
    example_batch_processing()

    print("\n" + "=" * 60)
    print("Examples complete!")
    print("=" * 60)
