"""Output writers for QC metrics and trimmed sequences."""

import gzip
import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd

logger = logging.getLogger(__name__)


def write_qc_metrics(
    metrics_list: List[Dict[str, Any]], output_dir: Path
) -> None:
    """
    Write per-read QC metrics to CSV and Parquet files.

    Args:
        metrics_list: List of per-read metric dictionaries
        output_dir: Output directory
    """
    if not metrics_list:
        logger.warning("No metrics to write")
        return

    # Create output directory
    qc_dir = output_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Convert to DataFrame
    df = pd.DataFrame(metrics_list)

    # Write CSV
    csv_path = qc_dir / "per_read_metrics.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Wrote per-read metrics CSV: {csv_path}")

    # Write Parquet
    try:
        parquet_path = qc_dir / "per_read_metrics.parquet"
        df.to_parquet(parquet_path, index=False, engine="pyarrow")
        logger.info(f"Wrote per-read metrics Parquet: {parquet_path}")
    except Exception as e:
        logger.warning(f"Failed to write Parquet file: {e}")


def write_summary_stats(summary: Dict[str, Any], output_dir: Path) -> None:
    """
    Write aggregated summary statistics to JSON file.

    Args:
        summary: Summary statistics dictionary
        output_dir: Output directory
    """
    # Create output directory
    qc_dir = output_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Write JSON
    json_path = qc_dir / "summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Wrote summary statistics: {json_path}")


def write_trimmed_fastq(
    sequences: List[Dict[str, Any]], output_path: Path
) -> None:
    """
    Write trimmed sequences to FASTQ file (gzipped).

    Args:
        sequences: List of sequence dictionaries with keys:
                   'read_id', 'seq', 'quals'
        output_path: Output FASTQ path (will be gzipped)
    """
    if not sequences:
        logger.warning("No sequences to write to FASTQ")
        return

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write gzipped FASTQ
    with gzip.open(output_path, "wt") as f:
        for seq_dict in sequences:
            read_id = seq_dict["read_id"]
            seq = seq_dict["seq"]
            quals = seq_dict["quals"]

            # Convert quality scores to ASCII
            qual_string = "".join(chr(q + 33) for q in quals)

            # Write FASTQ record
            f.write(f"@{read_id}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual_string}\n")

    logger.info(f"Wrote {len(sequences)} trimmed sequences to FASTQ: {output_path}")


def write_trimmed_fasta(
    sequences: List[Dict[str, Any]], output_path: Path
) -> None:
    """
    Write trimmed sequences to FASTA file (gzipped).

    Args:
        sequences: List of sequence dictionaries with keys:
                   'read_id', 'seq'
        output_path: Output FASTA path (will be gzipped)
    """
    if not sequences:
        logger.warning("No sequences to write to FASTA")
        return

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write gzipped FASTA
    with gzip.open(output_path, "wt") as f:
        for seq_dict in sequences:
            read_id = seq_dict["read_id"]
            seq = seq_dict["seq"]

            # Write FASTA record
            f.write(f">{read_id}\n")
            f.write(f"{seq}\n")

    logger.info(f"Wrote {len(sequences)} trimmed sequences to FASTA: {output_path}")


def write_base_call_annotations(
    sample_id: str,
    annotations: List[Dict[str, Any]],
    output_dir: Path
) -> None:
    """
    Write per-base call annotations to CSV file.

    Args:
        sample_id: Sample identifier
        annotations: List of per-base annotation dictionaries
        output_dir: Output directory
    """
    if not annotations:
        logger.warning(f"No base call annotations to write for {sample_id}")
        return

    # Create output directory
    base_calls_dir = output_dir / "base_calls"
    base_calls_dir.mkdir(parents=True, exist_ok=True)

    # Convert to DataFrame
    df = pd.DataFrame(annotations)

    # Add sample_id column
    df.insert(0, 'sample_id', sample_id)

    # Write CSV
    csv_path = base_calls_dir / f"{sample_id}_base_calls.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Wrote base call annotations: {csv_path}")


def write_all_base_call_annotations(
    annotations_dict: Dict[str, List[Dict[str, Any]]],
    output_dir: Path
) -> None:
    """
    Write all per-base call annotations to a single combined CSV file.

    Args:
        annotations_dict: Dictionary mapping sample_id -> annotation list
        output_dir: Output directory
    """
    if not annotations_dict:
        logger.warning("No base call annotations to write")
        return

    # Create output directory
    base_calls_dir = output_dir / "base_calls"
    base_calls_dir.mkdir(parents=True, exist_ok=True)

    # Combine all annotations
    all_annotations = []
    for sample_id, annotations in annotations_dict.items():
        for annotation in annotations:
            # Add sample_id to each annotation
            ann_copy = annotation.copy()
            ann_copy['sample_id'] = sample_id
            all_annotations.append(ann_copy)

    # Convert to DataFrame
    df = pd.DataFrame(all_annotations)

    # Reorder columns to have sample_id first
    cols = ['sample_id'] + [col for col in df.columns if col != 'sample_id']
    df = df[cols]

    # Write CSV
    csv_path = base_calls_dir / "all_base_calls.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Wrote combined base call annotations: {csv_path}")

    # Also write Parquet for efficient storage
    try:
        parquet_path = base_calls_dir / "all_base_calls.parquet"
        df.to_parquet(parquet_path, index=False, engine="pyarrow")
        logger.info(f"Wrote base call annotations Parquet: {parquet_path}")
    except Exception as e:
        logger.warning(f"Failed to write Parquet file: {e}")


def setup_logging(output_dir: Path, verbose: bool = False, quiet: bool = False) -> None:
    """
    Setup logging to both console and file.

    Args:
        output_dir: Output directory for log file
        verbose: Enable verbose (DEBUG) logging
        quiet: Suppress console output (only log to file)
    """
    # Create logs directory
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = logs_dir / "run.log"

    # Determine log level
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Remove existing handlers
    logger.handlers.clear()

    # File handler (always INFO or above)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Console handler (unless quiet)
    if not quiet:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_formatter = logging.Formatter("%(levelname)s: %(message)s")
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

    logger.info(f"Logging to: {log_file}")
