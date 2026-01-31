"""Command-line interface for Sanger QC and trimming tool."""

import sys
import logging
from pathlib import Path
from typing import List, Optional
import typer
from tqdm import tqdm

from .io_utils import discover_files, get_sample_id, parse_sequence_file, make_read_id
from .trim import apply_trim
from .qc import compute_qc_metrics, compute_summary_stats
from .writers import (
    write_qc_metrics,
    write_summary_stats,
    write_trimmed_fastq,
    write_trimmed_fasta,
    write_all_base_call_annotations,
    setup_logging,
)
from .ambiguous_calling import create_ambiguous_caller, AmbiguousBaseCaller
from .plots import (
    plot_sequence_trim,
    plot_sequence_trim_interactive,
    plot_multiple_sequences,
    plot_summary_histograms,
    plot_length_comparison,
    plot_ambiguous_calling_interactive,
)

app = typer.Typer(help="QC and trimming tool for Sanger sequencing reads")
logger = logging.getLogger(__name__)


@app.command()
def qc(
    inputs: List[str] = typer.Argument(..., help="Input files or directories"),
    output: str = typer.Option(..., "-o", "--output", help="Output directory"),
    qthreshold: int = typer.Option(20, "--qthreshold", help="Quality threshold for trimming"),
    method: str = typer.Option("mott", "--method", help="Trimming method (mott or ends)"),
    min_length: int = typer.Option(50, "--min-length", help="Minimum acceptable trimmed length"),
    recursive: bool = typer.Option(False, "--recursive", "-r", help="Recursively search directories"),
    plots: bool = typer.Option(False, "--plots", help="Generate quality and trimming plots"),
    ambiguous_calling: bool = typer.Option(False, "--ambiguous-calling", help="Enable ambiguous base calling with IUPAC codes"),
    clonal_context: bool = typer.Option(True, "--clonal-context/--mixed-context", help="Sample context (clonal vs mixed)"),
    spr_noise: float = typer.Option(0.20, "--spr-noise", help="Max SPR for noise threshold"),
    spr_het_low: float = typer.Option(0.33, "--spr-het-low", help="Lower SPR for heterozygous calls"),
    spr_het_high: float = typer.Option(0.67, "--spr-het-high", help="Upper SPR for heterozygous calls"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose logging"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress console output"),
):
    """
    Perform QC analysis on Sanger sequencing reads.
    """
    output_dir = Path(output)
    setup_logging(output_dir, verbose, quiet)

    logger.info("Starting QC analysis")
    logger.info(f"Parameters: qthreshold={qthreshold}, method={method}, min_length={min_length}")

    if ambiguous_calling:
        logger.info(f"Ambiguous calling enabled: clonal_context={clonal_context}, spr_noise={spr_noise}, spr_het={spr_het_low}-{spr_het_high}")

    # Discover files
    files = discover_files(inputs, recursive)

    if not files:
        logger.error("No valid input files found")
        raise typer.Exit(code=1)

    # Initialize ambiguous base caller if enabled
    caller = None
    if ambiguous_calling:
        caller = create_ambiguous_caller(
            clonal_context=clonal_context,
            spr_noise_max=spr_noise,
            spr_het_low=spr_het_low,
            spr_het_high=spr_het_high
        )

    # Process files
    metrics_list = []
    sequences_data = []  # For plotting
    base_call_annotations = {}  # For ambiguous calling
    skipped_count = 0

    for file_path, file_format in tqdm(files, desc="Processing files"):
        result = parse_sequence_file(file_path, file_format)

        if result is None:
            skipped_count += 1
            continue

        seq, quals = result
        sample_id = get_sample_id(file_path)

        # Perform ambiguous base calling if enabled (only for AB1 files)
        recalled_seq = seq
        if ambiguous_calling and caller and file_format == "ab1":
            try:
                base_calls = caller.call_bases(file_path, seq, quals)
                recalled_seq = caller.base_calls_to_sequence(base_calls)
                annotations = caller.base_calls_to_annotations(base_calls)
                base_call_annotations[sample_id] = annotations
                logger.debug(f"Ambiguous calling completed for {sample_id}: {len(base_calls)} bases")
            except Exception as e:
                logger.warning(f"Ambiguous calling failed for {sample_id}: {e}")
                recalled_seq = seq

        # Compute trimming coordinates (use recalled sequence)
        _, _, trim_start, trim_end = apply_trim(recalled_seq, quals, method, qthreshold)

        # Compute QC metrics (use recalled sequence)
        metrics = compute_qc_metrics(
            sample_id=sample_id,
            source_file=str(file_path),
            file_format=file_format,
            seq=recalled_seq,
            quals=quals,
            trim_start=trim_start,
            trim_end=trim_end,
            qthreshold=qthreshold,
            min_length=min_length,
        )
        metrics_list.append(metrics)

        # Store data for plotting
        if plots:
            sequences_data.append({
                'sample_id': sample_id,
                'quals': quals,
                'trim_start': trim_start,
                'trim_end': trim_end,
                'qthreshold': qthreshold,
            })

    logger.info(f"Processed {len(metrics_list)} reads successfully")
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} files due to errors")

    if not metrics_list:
        logger.error("No valid reads were processed")
        raise typer.Exit(code=1)

    # Write outputs
    write_qc_metrics(metrics_list, output_dir)
    summary = compute_summary_stats(metrics_list)
    write_summary_stats(summary, output_dir)

    # Write base call annotations if ambiguous calling was enabled
    if ambiguous_calling and base_call_annotations:
        write_all_base_call_annotations(base_call_annotations, output_dir)
        logger.info(f"Wrote base call annotations for {len(base_call_annotations)} samples")

    # Generate plots if requested
    if plots:
        logger.info("Generating plots...")
        plots_dir = output_dir / "plots"

        # Plot individual sequences (first 10)
        for seq_data in sequences_data[:10]:
            # Static matplotlib plot
            plot_path_png = plots_dir / f"{seq_data['sample_id']}_trim.png"
            plot_sequence_trim(
                sample_id=seq_data['sample_id'],
                quals=seq_data['quals'],
                trim_start=seq_data['trim_start'],
                trim_end=seq_data['trim_end'],
                qthreshold=seq_data['qthreshold'],
                output_path=plot_path_png,
            )

            # Interactive Plotly plot
            plot_path_html = plots_dir / f"{seq_data['sample_id']}_trim_interactive.html"
            plot_sequence_trim_interactive(
                sample_id=seq_data['sample_id'],
                quals=seq_data['quals'],
                trim_start=seq_data['trim_start'],
                trim_end=seq_data['trim_end'],
                qthreshold=seq_data['qthreshold'],
                output_path=plot_path_html,
            )

        # Plot multi-sequence overview
        if sequences_data:
            plot_multiple_sequences(sequences_data, plots_dir / "sequences_overview.png")

        # Plot summary histograms
        plot_summary_histograms(metrics_list, plots_dir / "summary_histograms.png")

        # Plot length comparison
        plot_length_comparison(metrics_list, plots_dir / "length_comparison.png")

        logger.info(f"Generated plots in: {plots_dir}")

    # Print summary
    logger.info("\n=== Summary Statistics ===")
    logger.info(f"Total reads: {summary['total_reads']}")
    logger.info(f"Mean raw length: {summary['mean_raw_length']:.1f}")
    logger.info(f"Mean trimmed length: {summary['mean_trimmed_length']:.1f}")
    logger.info(f"Mean quality: {summary['mean_mean_q']:.1f}")
    logger.info(f"Passed min length: {summary['reads_passed_minlen']} ({summary['pct_passed']:.1f}%)")

    logger.info("\n=== Output Files ===")
    logger.info(f"Per-read metrics: {output_dir}/qc/per_read_metrics.csv")
    logger.info(f"Per-read metrics: {output_dir}/qc/per_read_metrics.parquet")
    logger.info(f"Summary stats: {output_dir}/qc/summary.json")
    if ambiguous_calling and base_call_annotations:
        logger.info(f"Base call annotations: {output_dir}/base_calls/all_base_calls.csv")
        logger.info(f"Base call annotations: {output_dir}/base_calls/all_base_calls.parquet")
    if plots:
        logger.info(f"Plots: {output_dir}/plots/")
    logger.info(f"Log file: {output_dir}/logs/run.log")


@app.command()
def trim(
    inputs: List[str] = typer.Argument(..., help="Input files or directories"),
    output: str = typer.Option(..., "-o", "--output", help="Output directory"),
    qthreshold: int = typer.Option(20, "--qthreshold", help="Quality threshold for trimming"),
    method: str = typer.Option("mott", "--method", help="Trimming method (mott or ends)"),
    min_length: int = typer.Option(50, "--min-length", help="Minimum acceptable trimmed length"),
    out_fastq: Optional[str] = typer.Option(None, "--out-fastq", help="Output FASTQ file path"),
    out_fasta: Optional[str] = typer.Option(None, "--out-fasta", help="Output FASTA file path"),
    recursive: bool = typer.Option(False, "--recursive", "-r", help="Recursively search directories"),
    ambiguous_calling: bool = typer.Option(False, "--ambiguous-calling", help="Enable ambiguous base calling with IUPAC codes"),
    clonal_context: bool = typer.Option(True, "--clonal-context/--mixed-context", help="Sample context (clonal vs mixed)"),
    spr_noise: float = typer.Option(0.20, "--spr-noise", help="Max SPR for noise threshold"),
    spr_het_low: float = typer.Option(0.33, "--spr-het-low", help="Lower SPR for heterozygous calls"),
    spr_het_high: float = typer.Option(0.67, "--spr-het-high", help="Upper SPR for heterozygous calls"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose logging"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress console output"),
):
    """
    Trim Sanger sequencing reads and output trimmed sequences.
    """
    output_dir = Path(output)
    setup_logging(output_dir, verbose, quiet)

    logger.info("Starting trimming")
    logger.info(f"Parameters: qthreshold={qthreshold}, method={method}, min_length={min_length}")

    if ambiguous_calling:
        logger.info(f"Ambiguous calling enabled: clonal_context={clonal_context}, spr_noise={spr_noise}, spr_het={spr_het_low}-{spr_het_high}")

    # Discover files
    files = discover_files(inputs, recursive)

    if not files:
        logger.error("No valid input files found")
        raise typer.Exit(code=1)

    # Initialize ambiguous base caller if enabled
    caller = None
    if ambiguous_calling:
        caller = create_ambiguous_caller(
            clonal_context=clonal_context,
            spr_noise_max=spr_noise,
            spr_het_low=spr_het_low,
            spr_het_high=spr_het_high
        )

    # Process files
    trimmed_sequences = []
    base_call_annotations = {}  # For ambiguous calling
    skipped_count = 0

    for file_path, file_format in tqdm(files, desc="Trimming files"):
        result = parse_sequence_file(file_path, file_format)

        if result is None:
            skipped_count += 1
            continue

        seq, quals = result
        sample_id = get_sample_id(file_path)

        # Perform ambiguous base calling if enabled (only for AB1 files)
        recalled_seq = seq
        if ambiguous_calling and caller and file_format == "ab1":
            try:
                base_calls = caller.call_bases(file_path, seq, quals)
                recalled_seq = caller.base_calls_to_sequence(base_calls)
                annotations = caller.base_calls_to_annotations(base_calls)
                base_call_annotations[sample_id] = annotations
                logger.debug(f"Ambiguous calling completed for {sample_id}: {len(base_calls)} bases")
            except Exception as e:
                logger.warning(f"Ambiguous calling failed for {sample_id}: {e}")
                recalled_seq = seq

        # Apply trimming (use recalled sequence)
        trimmed_seq, trimmed_quals, trim_start, trim_end = apply_trim(
            recalled_seq, quals, method, qthreshold
        )

        # Create read ID
        read_id = make_read_id(sample_id, trim_start, trim_end)

        trimmed_sequences.append({
            "read_id": read_id,
            "seq": trimmed_seq,
            "quals": trimmed_quals,
        })

    logger.info(f"Trimmed {len(trimmed_sequences)} reads successfully")
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} files due to errors")

    if not trimmed_sequences:
        logger.error("No valid reads were processed")
        raise typer.Exit(code=1)

    # Write base call annotations if ambiguous calling was enabled
    if ambiguous_calling and base_call_annotations:
        write_all_base_call_annotations(base_call_annotations, output_dir)
        logger.info(f"Wrote base call annotations for {len(base_call_annotations)} samples")

    # Write outputs
    if out_fastq:
        fastq_path = Path(out_fastq)
        write_trimmed_fastq(trimmed_sequences, fastq_path)
    else:
        # Default FASTQ output
        fastq_path = output_dir / "trim" / "trimmed.fastq.gz"
        write_trimmed_fastq(trimmed_sequences, fastq_path)

    if out_fasta:
        fasta_path = Path(out_fasta)
        write_trimmed_fasta(trimmed_sequences, fasta_path)

    logger.info("\n=== Output Files ===")
    if ambiguous_calling and base_call_annotations:
        logger.info(f"Base call annotations: {output_dir}/base_calls/all_base_calls.csv")
        logger.info(f"Base call annotations: {output_dir}/base_calls/all_base_calls.parquet")
    if out_fastq:
        logger.info(f"Trimmed FASTQ: {out_fastq}")
    else:
        logger.info(f"Trimmed FASTQ: {fastq_path}")
    if out_fasta:
        logger.info(f"Trimmed FASTA: {out_fasta}")
    logger.info(f"Log file: {output_dir}/logs/run.log")


@app.command()
def all(
    inputs: List[str] = typer.Argument(..., help="Input files or directories"),
    output: str = typer.Option(..., "-o", "--output", help="Output directory"),
    qthreshold: int = typer.Option(20, "--qthreshold", help="Quality threshold for trimming"),
    method: str = typer.Option("mott", "--method", help="Trimming method (mott or ends)"),
    min_length: int = typer.Option(50, "--min-length", help="Minimum acceptable trimmed length"),
    out_fastq: Optional[str] = typer.Option(None, "--out-fastq", help="Output FASTQ file path"),
    out_fasta: Optional[str] = typer.Option(None, "--out-fasta", help="Output FASTA file path"),
    recursive: bool = typer.Option(False, "--recursive", "-r", help="Recursively search directories"),
    plots: bool = typer.Option(False, "--plots", help="Generate quality and trimming plots"),
    ambiguous_calling: bool = typer.Option(False, "--ambiguous-calling", help="Enable ambiguous base calling with IUPAC codes"),
    clonal_context: bool = typer.Option(True, "--clonal-context/--mixed-context", help="Sample context (clonal vs mixed)"),
    spr_noise: float = typer.Option(0.20, "--spr-noise", help="Max SPR for noise threshold"),
    spr_het_low: float = typer.Option(0.33, "--spr-het-low", help="Lower SPR for heterozygous calls"),
    spr_het_high: float = typer.Option(0.67, "--spr-het-high", help="Upper SPR for heterozygous calls"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose logging"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress console output"),
):
    """
    Perform both QC analysis and trimming.
    """
    output_dir = Path(output)
    setup_logging(output_dir, verbose, quiet)

    logger.info("Starting QC and trimming")
    logger.info(f"Parameters: qthreshold={qthreshold}, method={method}, min_length={min_length}")

    if ambiguous_calling:
        logger.info(f"Ambiguous calling enabled: clonal_context={clonal_context}, spr_noise={spr_noise}, spr_het={spr_het_low}-{spr_het_high}")

    # Discover files
    files = discover_files(inputs, recursive)

    if not files:
        logger.error("No valid input files found")
        raise typer.Exit(code=1)

    # Initialize ambiguous base caller if enabled
    caller = None
    if ambiguous_calling:
        caller = create_ambiguous_caller(
            clonal_context=clonal_context,
            spr_noise_max=spr_noise,
            spr_het_low=spr_het_low,
            spr_het_high=spr_het_high
        )

    # Process files
    metrics_list = []
    trimmed_sequences = []
    sequences_data = []  # For plotting
    base_call_annotations = {}  # For ambiguous calling
    base_calls_for_plots = {}  # Store BaseCall objects for plotting
    skipped_count = 0

    for file_path, file_format in tqdm(files, desc="Processing files"):
        result = parse_sequence_file(file_path, file_format)

        if result is None:
            skipped_count += 1
            continue

        seq, quals = result
        sample_id = get_sample_id(file_path)

        # Perform ambiguous base calling if enabled (only for AB1 files)
        recalled_seq = seq
        base_calls_obj = None
        if ambiguous_calling and caller and file_format == "ab1":
            try:
                base_calls_obj = caller.call_bases(file_path, seq, quals)
                recalled_seq = caller.base_calls_to_sequence(base_calls_obj)
                annotations = caller.base_calls_to_annotations(base_calls_obj)
                base_call_annotations[sample_id] = annotations
                # Store for plotting
                if plots:
                    base_calls_for_plots[sample_id] = base_calls_obj
                logger.debug(f"Ambiguous calling completed for {sample_id}: {len(base_calls_obj)} bases")
            except Exception as e:
                logger.warning(f"Ambiguous calling failed for {sample_id}: {e}")
                # Fall back to original sequence
                recalled_seq = seq

        # Apply trimming (use recalled sequence if available)
        trimmed_seq, trimmed_quals, trim_start, trim_end = apply_trim(
            recalled_seq, quals, method, qthreshold
        )

        # Compute QC metrics (use recalled sequence)
        metrics = compute_qc_metrics(
            sample_id=sample_id,
            source_file=str(file_path),
            file_format=file_format,
            seq=recalled_seq,
            quals=quals,
            trim_start=trim_start,
            trim_end=trim_end,
            qthreshold=qthreshold,
            min_length=min_length,
        )
        metrics_list.append(metrics)

        # Store trimmed sequence
        read_id = make_read_id(sample_id, trim_start, trim_end)
        trimmed_sequences.append({
            "read_id": read_id,
            "seq": trimmed_seq,
            "quals": trimmed_quals,
        })

        # Store data for plotting
        if plots:
            sequences_data.append({
                'sample_id': sample_id,
                'quals': quals,
                'trim_start': trim_start,
                'trim_end': trim_end,
                'qthreshold': qthreshold,
            })

    logger.info(f"Processed {len(metrics_list)} reads successfully")
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} files due to errors")

    if not metrics_list:
        logger.error("No valid reads were processed")
        raise typer.Exit(code=1)

    # Write QC outputs
    write_qc_metrics(metrics_list, output_dir)
    summary = compute_summary_stats(metrics_list)
    write_summary_stats(summary, output_dir)

    # Write base call annotations if ambiguous calling was enabled
    if ambiguous_calling and base_call_annotations:
        write_all_base_call_annotations(base_call_annotations, output_dir)
        logger.info(f"Wrote base call annotations for {len(base_call_annotations)} samples")

    # Write trimmed sequences
    if out_fastq:
        fastq_path = Path(out_fastq)
        write_trimmed_fastq(trimmed_sequences, fastq_path)
    else:
        fastq_path = output_dir / "trim" / "trimmed.fastq.gz"
        write_trimmed_fastq(trimmed_sequences, fastq_path)

    if out_fasta:
        fasta_path = Path(out_fasta)
        write_trimmed_fasta(trimmed_sequences, fasta_path)

    # Generate plots if requested
    if plots:
        logger.info("Generating plots...")
        plots_dir = output_dir / "plots"

        # Plot individual sequences (first 10)
        for seq_data in sequences_data[:10]:
            # Static matplotlib plot
            plot_path_png = plots_dir / f"{seq_data['sample_id']}_trim.png"
            plot_sequence_trim(
                sample_id=seq_data['sample_id'],
                quals=seq_data['quals'],
                trim_start=seq_data['trim_start'],
                trim_end=seq_data['trim_end'],
                qthreshold=seq_data['qthreshold'],
                output_path=plot_path_png,
            )

            # Interactive Plotly plot
            plot_path_html = plots_dir / f"{seq_data['sample_id']}_trim_interactive.html"
            plot_sequence_trim_interactive(
                sample_id=seq_data['sample_id'],
                quals=seq_data['quals'],
                trim_start=seq_data['trim_start'],
                trim_end=seq_data['trim_end'],
                qthreshold=seq_data['qthreshold'],
                output_path=plot_path_html,
            )

            # Ambiguous calling plot (if enabled and data available)
            if ambiguous_calling and seq_data['sample_id'] in base_calls_for_plots:
                plot_path_ambig = plots_dir / f"{seq_data['sample_id']}_ambiguous_interactive.html"
                plot_ambiguous_calling_interactive(
                    sample_id=seq_data['sample_id'],
                    quals=seq_data['quals'],
                    base_calls=base_calls_for_plots[seq_data['sample_id']],
                    trim_start=seq_data['trim_start'],
                    trim_end=seq_data['trim_end'],
                    qthreshold=seq_data['qthreshold'],
                    output_path=plot_path_ambig,
                )

        # Plot multi-sequence overview
        if sequences_data:
            plot_multiple_sequences(sequences_data, plots_dir / "sequences_overview.png")

        # Plot summary histograms
        plot_summary_histograms(metrics_list, plots_dir / "summary_histograms.png")

        # Plot length comparison
        plot_length_comparison(metrics_list, plots_dir / "length_comparison.png")

        logger.info(f"Generated plots in: {plots_dir}")

    # Print summary
    logger.info("\n=== Summary Statistics ===")
    logger.info(f"Total reads: {summary['total_reads']}")
    logger.info(f"Mean raw length: {summary['mean_raw_length']:.1f}")
    logger.info(f"Mean trimmed length: {summary['mean_trimmed_length']:.1f}")
    logger.info(f"Mean quality: {summary['mean_mean_q']:.1f}")
    logger.info(f"Passed min length: {summary['reads_passed_minlen']} ({summary['pct_passed']:.1f}%)")

    logger.info("\n=== Output Files ===")
    logger.info(f"Per-read metrics: {output_dir}/qc/per_read_metrics.csv")
    logger.info(f"Per-read metrics: {output_dir}/qc/per_read_metrics.parquet")
    logger.info(f"Summary stats: {output_dir}/qc/summary.json")
    if ambiguous_calling and base_call_annotations:
        logger.info(f"Base call annotations: {output_dir}/base_calls/all_base_calls.csv")
        logger.info(f"Base call annotations: {output_dir}/base_calls/all_base_calls.parquet")
    if out_fastq:
        logger.info(f"Trimmed FASTQ: {out_fastq}")
    else:
        logger.info(f"Trimmed FASTQ: {fastq_path}")
    if out_fasta:
        logger.info(f"Trimmed FASTA: {out_fasta}")
    if plots:
        logger.info(f"Plots: {output_dir}/plots/")
    logger.info(f"Log file: {output_dir}/logs/run.log")


if __name__ == "__main__":
    app()
