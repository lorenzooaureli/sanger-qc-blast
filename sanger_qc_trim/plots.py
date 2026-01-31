"""Plotting functions for visualizing sequence quality and trimming."""

import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

logger = logging.getLogger(__name__)


def plot_sequence_trim(
    sample_id: str,
    quals: List[int],
    trim_start: int,
    trim_end: int,
    qthreshold: int,
    output_path: Path,
) -> None:
    """
    Plot quality scores with trimmed region highlighted.

    Args:
        sample_id: Sample identifier
        quals: List of quality scores
        trim_start: Trim start position
        trim_end: Trim end position
        qthreshold: Quality threshold used
        output_path: Path to save plot
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    positions = np.arange(len(quals))

    # Plot quality scores
    ax.plot(positions, quals, 'k-', linewidth=1, alpha=0.7, label='Quality scores')

    # Highlight trimmed region (kept)
    if trim_end > trim_start:
        ax.axvspan(trim_start, trim_end - 1, alpha=0.3, color='green', label='Kept region')

    # Highlight discarded regions
    if trim_start > 0:
        ax.axvspan(0, trim_start - 1, alpha=0.3, color='red', label='Discarded (5\' end)')
    if trim_end < len(quals):
        ax.axvspan(trim_end, len(quals) - 1, alpha=0.3, color='red')

    # Add threshold line
    ax.axhline(y=qthreshold, color='blue', linestyle='--', linewidth=1.5,
               label=f'Q{qthreshold} threshold')

    # Add Q20 and Q30 reference lines
    ax.axhline(y=20, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax.axhline(y=30, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # Labels and title
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Quality Score (Phred)', fontsize=12)
    ax.set_title(f'{sample_id} - Quality Profile with Trimming', fontsize=14, fontweight='bold')

    # Add text annotations
    raw_length = len(quals)
    trimmed_length = trim_end - trim_start
    mean_q = np.mean(quals) if quals else 0

    textstr = f'Raw length: {raw_length} bp\n'
    textstr += f'Trimmed length: {trimmed_length} bp ({100*trimmed_length/raw_length:.1f}%)\n'
    textstr += f'Trim coordinates: [{trim_start}, {trim_end})\n'
    textstr += f'Mean quality: {mean_q:.1f}'

    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Set y-axis limits
    ax.set_ylim(0, max(50, max(quals) + 5) if quals else 50)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.debug(f"Saved quality plot: {output_path}")


def plot_multiple_sequences(
    sequences_data: List[Dict[str, Any]],
    output_path: Path,
    max_sequences: int = 10,
) -> None:
    """
    Plot multiple sequences in a grid layout.

    Args:
        sequences_data: List of dicts with 'sample_id', 'quals', 'trim_start',
                       'trim_end', 'qthreshold'
        output_path: Path to save plot
        max_sequences: Maximum number of sequences to plot
    """
    n_seqs = min(len(sequences_data), max_sequences)

    if n_seqs == 0:
        logger.warning("No sequences to plot")
        return

    # Calculate grid layout
    n_cols = 2
    n_rows = (n_seqs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4 * n_rows))

    # Handle single subplot case - ensure axes is always 2D array
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    for idx, seq_data in enumerate(sequences_data[:max_sequences]):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        sample_id = seq_data['sample_id']
        quals = seq_data['quals']
        trim_start = seq_data['trim_start']
        trim_end = seq_data['trim_end']
        qthreshold = seq_data['qthreshold']

        positions = np.arange(len(quals))

        # Plot quality scores
        ax.plot(positions, quals, 'k-', linewidth=0.8, alpha=0.7)

        # Highlight regions
        if trim_end > trim_start:
            ax.axvspan(trim_start, trim_end - 1, alpha=0.3, color='green')
        if trim_start > 0:
            ax.axvspan(0, trim_start - 1, alpha=0.2, color='red')
        if trim_end < len(quals):
            ax.axvspan(trim_end, len(quals) - 1, alpha=0.2, color='red')

        # Threshold line
        ax.axhline(y=qthreshold, color='blue', linestyle='--', linewidth=1, alpha=0.7)

        # Labels
        ax.set_xlabel('Position (bp)', fontsize=9)
        ax.set_ylabel('Quality', fontsize=9)
        ax.set_title(f'{sample_id} ({trim_end - trim_start} bp kept)', fontsize=10)
        ax.grid(True, alpha=0.2)

    # Hide unused subplots
    for idx in range(n_seqs, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    # Add legend to the figure
    green_patch = mpatches.Patch(color='green', alpha=0.3, label='Kept region')
    red_patch = mpatches.Patch(color='red', alpha=0.2, label='Discarded')
    blue_line = mpatches.Patch(color='blue', label='Quality threshold')

    fig.legend(handles=[green_patch, red_patch, blue_line],
               loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.98))

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved multi-sequence plot: {output_path}")


def plot_summary_histograms(
    metrics_list: List[Dict[str, Any]],
    output_path: Path,
) -> None:
    """
    Plot summary histograms of key metrics.

    Args:
        metrics_list: List of per-read metrics
        output_path: Path to save plot
    """
    if not metrics_list:
        logger.warning("No metrics to plot")
        return

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    # Extract data
    raw_lengths = [m['raw_length'] for m in metrics_list]
    trimmed_lengths = [m['trimmed_length'] for m in metrics_list]
    mean_qs = [m['mean_q'] for m in metrics_list]
    pct_q20s = [m['pct_q20'] * 100 for m in metrics_list]
    gc_percents = [m['gc_percent'] for m in metrics_list]
    expected_errors = [m['expected_errors'] for m in metrics_list]

    # Plot 1: Raw length distribution
    axes[0].hist(raw_lengths, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Raw Length (bp)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Raw Sequence Length Distribution')
    axes[0].grid(True, alpha=0.3)

    # Plot 2: Trimmed length distribution
    axes[1].hist(trimmed_lengths, bins=20, color='forestgreen', alpha=0.7, edgecolor='black')
    axes[1].set_xlabel('Trimmed Length (bp)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Trimmed Sequence Length Distribution')
    axes[1].grid(True, alpha=0.3)

    # Plot 3: Mean quality distribution
    axes[2].hist(mean_qs, bins=20, color='orange', alpha=0.7, edgecolor='black')
    axes[2].axvline(x=20, color='red', linestyle='--', label='Q20')
    axes[2].axvline(x=30, color='darkred', linestyle='--', label='Q30')
    axes[2].set_xlabel('Mean Quality Score')
    axes[2].set_ylabel('Count')
    axes[2].set_title('Mean Quality Distribution')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    # Plot 4: Q20 percentage
    axes[3].hist(pct_q20s, bins=20, color='purple', alpha=0.7, edgecolor='black')
    axes[3].set_xlabel('% Bases ≥ Q20')
    axes[3].set_ylabel('Count')
    axes[3].set_title('Q20 Percentage Distribution')
    axes[3].grid(True, alpha=0.3)

    # Plot 5: GC content
    axes[4].hist(gc_percents, bins=20, color='teal', alpha=0.7, edgecolor='black')
    axes[4].set_xlabel('GC Content (%)')
    axes[4].set_ylabel('Count')
    axes[4].set_title('GC Content Distribution')
    axes[4].grid(True, alpha=0.3)

    # Plot 6: Expected errors
    axes[5].hist(expected_errors, bins=20, color='crimson', alpha=0.7, edgecolor='black')
    axes[5].set_xlabel('Expected Errors')
    axes[5].set_ylabel('Count')
    axes[5].set_title('Expected Errors Distribution')
    axes[5].grid(True, alpha=0.3)

    plt.suptitle(f'QC Metrics Summary (n={len(metrics_list)} reads)',
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved summary histograms: {output_path}")


def plot_length_comparison(
    metrics_list: List[Dict[str, Any]],
    output_path: Path,
) -> None:
    """
    Plot before/after length comparison.

    Args:
        metrics_list: List of per-read metrics
        output_path: Path to save plot
    """
    if not metrics_list:
        logger.warning("No metrics to plot")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    sample_ids = [m['sample_id'] for m in metrics_list]
    raw_lengths = [m['raw_length'] for m in metrics_list]
    trimmed_lengths = [m['trimmed_length'] for m in metrics_list]

    x = np.arange(len(sample_ids))
    width = 0.35

    # Bar plot comparison
    bars1 = ax1.bar(x - width/2, raw_lengths, width, label='Raw',
                    color='steelblue', alpha=0.7, edgecolor='black')
    bars2 = ax1.bar(x + width/2, trimmed_lengths, width, label='Trimmed',
                    color='forestgreen', alpha=0.7, edgecolor='black')

    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Length (bp)')
    ax1.set_title('Sequence Length: Before vs After Trimming')
    ax1.set_xticks(x)
    ax1.set_xticklabels(sample_ids, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    # Scatter plot: raw vs trimmed
    ax2.scatter(raw_lengths, trimmed_lengths, alpha=0.6, s=100,
                color='purple', edgecolor='black', linewidth=0.5)

    # Add diagonal line
    max_len = max(max(raw_lengths), max(trimmed_lengths))
    ax2.plot([0, max_len], [0, max_len], 'r--', linewidth=2, alpha=0.5, label='No trimming')

    ax2.set_xlabel('Raw Length (bp)')
    ax2.set_ylabel('Trimmed Length (bp)')
    ax2.set_title('Raw vs Trimmed Length Correlation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Add statistics text
    retention = 100 * sum(trimmed_lengths) / sum(raw_lengths) if sum(raw_lengths) > 0 else 0
    textstr = f'Overall retention: {retention:.1f}%\n'
    textstr += f'Mean raw: {np.mean(raw_lengths):.1f} bp\n'
    textstr += f'Mean trimmed: {np.mean(trimmed_lengths):.1f} bp'

    ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved length comparison plot: {output_path}")


def plot_sequence_trim_interactive(
    sample_id: str,
    quals: List[int],
    trim_start: int,
    trim_end: int,
    qthreshold: int,
    output_path: Path,
) -> None:
    """
    Create an interactive Plotly plot of quality scores with trimmed region highlighted.

    Args:
        sample_id: Sample identifier
        quals: List of quality scores
        trim_start: Trim start position
        trim_end: Trim end position
        qthreshold: Quality threshold used
        output_path: Path to save HTML file
    """
    positions = list(range(len(quals)))

    # Create figure
    fig = go.Figure()

    # Add quality score line
    fig.add_trace(go.Scatter(
        x=positions,
        y=quals,
        mode='lines',
        name='Quality scores',
        line=dict(color='black', width=1.5),
        hovertemplate='Position: %{x}<br>Quality: %{y}<extra></extra>'
    ))

    # Add shaded regions for kept/discarded regions
    if trim_end > trim_start:
        # Kept region (green)
        fig.add_vrect(
            x0=trim_start, x1=trim_end - 1,
            fillcolor="green", opacity=0.2,
            layer="below", line_width=0,
            annotation_text="Kept region", annotation_position="top left"
        )

    # Discarded 5' region (red)
    if trim_start > 0:
        fig.add_vrect(
            x0=0, x1=trim_start - 1,
            fillcolor="red", opacity=0.2,
            layer="below", line_width=0,
            annotation_text="Discarded (5')", annotation_position="top left"
        )

    # Discarded 3' region (red)
    if trim_end < len(quals):
        fig.add_vrect(
            x0=trim_end, x1=len(quals) - 1,
            fillcolor="red", opacity=0.2,
            layer="below", line_width=0,
            annotation_text="Discarded (3')", annotation_position="top right"
        )

    # Add threshold line
    fig.add_hline(
        y=qthreshold,
        line_dash="dash",
        line_color="blue",
        annotation_text=f"Q{qthreshold} threshold",
        annotation_position="right"
    )

    # Add Q20 and Q30 reference lines
    fig.add_hline(y=20, line_dash="dot", line_color="gray", opacity=0.5)
    fig.add_hline(y=30, line_dash="dot", line_color="gray", opacity=0.5)

    # Calculate statistics
    raw_length = len(quals)
    trimmed_length = trim_end - trim_start
    mean_q = np.mean(quals) if quals else 0

    # Update layout
    fig.update_layout(
        title=dict(
            text=f'{sample_id} - Quality Profile with Trimming',
            font=dict(size=16, weight='bold')
        ),
        xaxis_title='Position (bp)',
        yaxis_title='Quality Score (Phred)',
        hovermode='x unified',
        showlegend=True,
        height=600,
        width=1200,
        annotations=[
            dict(
                text=f'Raw length: {raw_length} bp<br>'
                     f'Trimmed length: {trimmed_length} bp ({100*trimmed_length/raw_length:.1f}%)<br>'
                     f'Trim coordinates: [{trim_start}, {trim_end})<br>'
                     f'Mean quality: {mean_q:.1f}',
                xref='paper', yref='paper',
                x=0.02, y=0.98,
                xanchor='left', yanchor='top',
                showarrow=False,
                bgcolor='wheat',
                bordercolor='black',
                borderwidth=1,
                borderpad=10,
                font=dict(size=12)
            )
        ]
    )

    # Set y-axis range
    fig.update_yaxes(range=[0, max(50, max(quals) + 5) if quals else 50])

    # Save to HTML
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path))

    logger.debug(f"Saved interactive quality plot: {output_path}")


def plot_ambiguous_calling_interactive(
    sample_id: str,
    quals: List[int],
    base_calls: List[Any],  # List of BaseCall objects
    trim_start: int,
    trim_end: int,
    qthreshold: int,
    output_path: Path,
) -> None:
    """
    Create an interactive Plotly plot showing quality scores, trimming, and ambiguous base calls.

    Args:
        sample_id: Sample identifier
        quals: List of quality scores
        base_calls: List of BaseCall objects from ambiguous calling
        trim_start: Trim start position
        trim_end: Trim end position
        qthreshold: Quality threshold used
        output_path: Path to save HTML file
    """
    positions = list(range(len(quals)))

    # Extract information from base calls
    called_bases = [bc.called_base for bc in base_calls]
    call_modes = [bc.call_mode for bc in base_calls]
    sprs = [bc.spr for bc in base_calls]
    snrs = [bc.snr for bc in base_calls]
    allele_fracs = [bc.allele_fraction for bc in base_calls]
    primary_bases = [bc.primary_base for bc in base_calls]
    secondary_bases = [bc.secondary_base for bc in base_calls]
    flags_list = [', '.join(bc.flags) if bc.flags else '' for bc in base_calls]

    # Identify different call types
    ambiguous_positions = [i for i, bc in enumerate(base_calls) if bc.call_mode == 'ambiguous']
    n_positions = [i for i, bc in enumerate(base_calls) if bc.call_mode == 'N']
    single_positions = [i for i, bc in enumerate(base_calls) if bc.call_mode == 'single']

    # Create figure with subplots
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.7, 0.3],
        subplot_titles=('Quality Profile with Ambiguous Base Calling', 'Base Call Types'),
        vertical_spacing=0.12
    )

    # Top plot: Quality scores
    fig.add_trace(go.Scatter(
        x=positions,
        y=quals,
        mode='lines',
        name='Quality scores',
        line=dict(color='black', width=1.5),
        hovertemplate='Position: %{x}<br>Quality: %{y}<extra></extra>',
        showlegend=True
    ), row=1, col=1)

    # Mark ambiguous positions on quality plot
    if ambiguous_positions:
        fig.add_trace(go.Scatter(
            x=ambiguous_positions,
            y=[quals[i] for i in ambiguous_positions],
            mode='markers',
            name='Ambiguous (IUPAC)',
            marker=dict(
                size=8,
                color='orange',
                symbol='diamond',
                line=dict(width=1, color='darkorange')
            ),
            hovertemplate='<b>Ambiguous Call</b><br>' +
                          'Position: %{x}<br>' +
                          'Quality: %{y}<br>' +
                          'Base: %{customdata[0]}<br>' +
                          'SPR: %{customdata[1]:.3f}<br>' +
                          'SNR: %{customdata[2]:.2f}<br>' +
                          'Allele Frac: %{customdata[3]:.3f}<br>' +
                          'Primary: %{customdata[4]}<br>' +
                          'Secondary: %{customdata[5]}<br>' +
                          'Flags: %{customdata[6]}<extra></extra>',
            customdata=[[called_bases[i], sprs[i], snrs[i], allele_fracs[i],
                        primary_bases[i], secondary_bases[i], flags_list[i]]
                       for i in ambiguous_positions],
            showlegend=True
        ), row=1, col=1)

    # Mark N positions on quality plot
    if n_positions:
        fig.add_trace(go.Scatter(
            x=n_positions,
            y=[quals[i] for i in n_positions],
            mode='markers',
            name='No call (N)',
            marker=dict(
                size=8,
                color='red',
                symbol='x',
                line=dict(width=2)
            ),
            hovertemplate='<b>No Call (N)</b><br>' +
                          'Position: %{x}<br>' +
                          'Quality: %{y}<br>' +
                          'SPR: %{customdata[0]:.3f}<br>' +
                          'SNR: %{customdata[1]:.2f}<br>' +
                          'Flags: %{customdata[2]}<extra></extra>',
            customdata=[[sprs[i], snrs[i], flags_list[i]] for i in n_positions],
            showlegend=True
        ), row=1, col=1)

    # Add shaded regions for kept/discarded regions
    if trim_end > trim_start:
        fig.add_vrect(
            x0=trim_start, x1=trim_end - 1,
            fillcolor="green", opacity=0.15,
            layer="below", line_width=0,
            annotation_text="Kept region", annotation_position="top left",
            row=1, col=1
        )

    if trim_start > 0:
        fig.add_vrect(
            x0=0, x1=trim_start - 1,
            fillcolor="red", opacity=0.1,
            layer="below", line_width=0,
            annotation_text="Discarded (5')", annotation_position="top left",
            row=1, col=1
        )

    if trim_end < len(quals):
        fig.add_vrect(
            x0=trim_end, x1=len(quals) - 1,
            fillcolor="red", opacity=0.1,
            layer="below", line_width=0,
            annotation_text="Discarded (3')", annotation_position="top right",
            row=1, col=1
        )

    # Add threshold lines
    fig.add_hline(
        y=qthreshold,
        line_dash="dash",
        line_color="blue",
        annotation_text=f"Q{qthreshold}",
        annotation_position="right",
        row=1, col=1
    )
    fig.add_hline(y=20, line_dash="dot", line_color="gray", opacity=0.5, row=1, col=1)
    fig.add_hline(y=30, line_dash="dot", line_color="gray", opacity=0.5, row=1, col=1)

    # Bottom plot: Base call type indicator
    # Create a categorical plot showing call types
    call_type_numeric = []
    call_type_colors = []
    for mode in call_modes:
        if mode == 'single':
            call_type_numeric.append(1)
            call_type_colors.append('green')
        elif mode == 'ambiguous':
            call_type_numeric.append(2)
            call_type_colors.append('orange')
        else:  # 'N'
            call_type_numeric.append(0)
            call_type_colors.append('red')

    fig.add_trace(go.Bar(
        x=positions,
        y=call_type_numeric,
        marker=dict(
            color=call_type_colors,
            line=dict(width=0)
        ),
        name='Call types',
        hovertemplate='Position: %{x}<br>' +
                      'Base: %{customdata[0]}<br>' +
                      'Mode: %{customdata[1]}<br>' +
                      'SPR: %{customdata[2]:.3f}<extra></extra>',
        customdata=[[called_bases[i], call_modes[i], sprs[i]] for i in positions],
        showlegend=False
    ), row=2, col=1)

    # Calculate statistics
    raw_length = len(quals)
    trimmed_length = trim_end - trim_start
    mean_q = np.mean(quals) if quals else 0
    n_ambiguous = len(ambiguous_positions)
    n_no_calls = len(n_positions)
    n_single = len(single_positions)

    # Update layout
    fig.update_layout(
        title=dict(
            text=f'{sample_id} - Ambiguous Base Calling Analysis',
            font=dict(size=18)
        ),
        hovermode='x unified',
        showlegend=True,
        height=800,
        width=1400,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    # Update axes
    fig.update_xaxes(title_text="Position (bp)", row=1, col=1)
    fig.update_yaxes(title_text="Quality Score (Phred)", row=1, col=1)
    fig.update_xaxes(title_text="Position (bp)", row=2, col=1)
    fig.update_yaxes(
        title_text="Call Type",
        tickvals=[0, 1, 2],
        ticktext=['N', 'Single', 'IUPAC'],
        row=2, col=1
    )

    # Set y-axis range for quality plot
    fig.update_yaxes(range=[0, max(50, max(quals) + 5) if quals else 50], row=1, col=1)
    fig.update_yaxes(range=[-0.5, 2.5], row=2, col=1)

    # Add statistics annotation
    stats_text = (
        f'<b>Statistics</b><br>'
        f'Raw length: {raw_length} bp<br>'
        f'Trimmed: {trimmed_length} bp ({100*trimmed_length/raw_length:.1f}%)<br>'
        f'Mean quality: {mean_q:.1f}<br>'
        f'<br>'
        f'<b>Base Calls</b><br>'
        f'Single: {n_single} ({100*n_single/raw_length:.1f}%)<br>'
        f'Ambiguous (IUPAC): {n_ambiguous} ({100*n_ambiguous/raw_length:.1f}%)<br>'
        f'No call (N): {n_no_calls} ({100*n_no_calls/raw_length:.1f}%)'
    )

    fig.add_annotation(
        text=stats_text,
        xref='paper', yref='paper',
        x=0.02, y=0.98,
        xanchor='left', yanchor='top',
        showarrow=False,
        bgcolor='wheat',
        bordercolor='black',
        borderwidth=1,
        borderpad=10,
        font=dict(size=11),
        align='left'
    )

    # Save to HTML
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path))

    logger.info(f"Saved ambiguous calling interactive plot: {output_path}")
