# Sanger QC & Trim

A production-ready command-line tool for quality control and trimming of Sanger sequencing reads from `.ab1` (ABI chromatogram) and `.phd.1` (Phred base-caller) files.

## Features

- **Multi-format support**: Automatically detects and processes `.ab1` and `.phd.1` files
- **Comprehensive QC metrics**: Per-read and aggregated statistics including quality scores, GC content, expected errors, and more
- **Smart trimming algorithms**:
  - **Mott/Kadane** (default): Maximum-score subsequence approach for optimal trim point detection
  - **Ends**: Hard clipping from 5' and 3' ends to first/last high-quality base
- **Ambiguous base calling**: Detect and report mixed signals using IUPAC codes (NEW!)
  - Heterozygous SNP detection with allele frequency estimation
  - Peak intensity analysis from AB1 chromatograms
  - Context-aware calling (clonal vs. mixed samples)
  - Per-base annotations with SPR, SNR, and quality metrics
- **Visual plots**: Generate quality plots showing trimmed regions, summary histograms, and length comparisons
- **Multiple output formats**: CSV, Parquet, JSON for metrics; FASTQ and FASTA for trimmed sequences; PNG plots
- **Production-ready**: Handles errors gracefully, provides detailed logging, and works on large datasets

## What are .ab1 and .phd.1 files?

- **`.ab1` files**: Binary chromatogram files from ABI/Applied Biosystems Sanger sequencers containing raw trace data, base calls, and quality scores
- **`.phd.1` files**: Plain-text output from the Phred base-caller containing sequences and per-base quality scores in a human-readable format

## Installation

### Using pip (recommended)

```bash
# Install from source directory
pip install .

# Or install in development mode
pip install -e .

# Install with development dependencies
pip install -e ".[dev]"
```

### Using pipx (isolated environment)

```bash
pipx install .
```

### Requirements

- Python 3.8 or higher
- Dependencies (automatically installed):
  - biopython ≥1.79
  - pandas ≥1.3.0
  - numpy ≥1.21.0
  - pyarrow ≥6.0.0
  - typer ≥0.9.0
  - tqdm ≥4.62.0
  - matplotlib ≥3.5.0

## Quick Start

```bash
# QC + Trimming on a directory (Mott algorithm, Q20 threshold, min length 50 bp)
sangerqc all seq_data/Y19/ -o output --recursive

# QC analysis only with Q30 threshold
sangerqc qc data/ -o results --qthreshold 30 --recursive

# Trim sequences using ends method and output to custom FASTQ
sangerqc trim data/ -o results --method ends --out-fastq trimmed_reads.fastq.gz
```

## Usage

### Command Structure

The tool provides three subcommands:

```bash
sangerqc qc <INPUT...> -o OUTDIR [OPTIONS]    # QC analysis only
sangerqc trim <INPUT...> -o OUTDIR [OPTIONS]  # Trimming only
sangerqc all <INPUT...> -o OUTDIR [OPTIONS]   # Both QC and trimming
```

### Common Options

- `INPUT...`: One or more input files or directories
- `-o, --output OUTDIR`: Output directory (required)
- `--qthreshold INT`: Quality threshold for trimming (default: 20)
- `--method {mott,ends}`: Trimming method (default: mott)
- `--min-length INT`: Minimum acceptable trimmed length (default: 50)
- `-r, --recursive`: Recursively search directories
- `--plots`: Generate quality and trimming visualization plots
- `-v, --verbose`: Verbose logging (DEBUG level)
- `-q, --quiet`: Suppress console output (log to file only)

### Ambiguous Calling Options (NEW!)

- `--ambiguous-calling`: Enable ambiguous base calling with IUPAC codes
- `--clonal-context` / `--mixed-context`: Sample context (clonal: default, mixed: for diploid/environmental)
- `--spr-noise FLOAT`: Max SPR for noise threshold (default: 0.20)
- `--spr-het-low FLOAT`: Lower SPR bound for heterozygous calls (default: 0.33)
- `--spr-het-high FLOAT`: Upper SPR bound for heterozygous calls (default: 0.67)

See [AMBIGUOUS_CALLING.md](AMBIGUOUS_CALLING.md) for detailed documentation on ambiguous base calling features.

### Trim-Specific Options

- `--out-fastq PATH`: Output FASTQ file path (default: OUTDIR/trim/trimmed.fastq.gz)
- `--out-fasta PATH`: Optional FASTA output file path

## Trimming Methods

### Mott/Kadane Algorithm (default)

Uses a maximum-score subsequence approach where each base contributes a weight of `Q - T` (quality minus threshold). Keeps the contiguous region with the highest cumulative score.

**Best for**: Optimal quality-based trimming that balances read length and quality

**Example**: Given qualities `[10, 10, 30, 30, 30, 10]` with threshold 20:
- Weights: `[-10, -10, 10, 10, 10, -10]`
- Best region: indices 2-5 (bases with Q=30)
- Result: `trim_start=2, trim_end=5`

### Ends Algorithm

Hard clips from the 5' and 3' ends to the first and last base meeting the quality threshold.

**Best for**: Conservative trimming that preserves internal low-quality regions if flanked by high-quality bases

**Example**: Given qualities `[15, 25, 25, 15]` with threshold 20:
- First Q≥20: index 1
- Last Q≥20: index 2
- Result: `trim_start=1, trim_end=3`

## QC Metrics

### Per-Read Metrics

Each read generates the following metrics:

| Metric | Description |
|--------|-------------|
| `sample_id` | Sample identifier (filename without extension) |
| `source_file` | Path to source file |
| `format` | File format (ab1 or phd.1) |
| `raw_length` | Original sequence length |
| `mean_q` | Mean quality score |
| `median_q` | Median quality score |
| `pct_q20` | Fraction of bases with Q≥20 |
| `pct_q30` | Fraction of bases with Q≥30 |
| `gc_percent` | GC content (excluding N bases) |
| `n_count` | Number of N (ambiguous) bases |
| `expected_errors` | Sum of error probabilities: Σ10^(-Q/10) |
| `hq_longest_stretch_len` | Longest contiguous high-quality region |
| `trim_start` | Trimming start position (0-based) |
| `trim_end` | Trimming end position (0-based, exclusive) |
| `trimmed_length` | Length after trimming |
| `passed_minlen` | "yes" if trimmed_length ≥ min_length, else "no" |

### Aggregated Summary

Summary statistics across all reads include:

- Total read count
- Mean and median raw/trimmed lengths
- Mean quality metrics (Q20%, Q30%, GC%, expected errors)
- Pass/fail counts and percentages

## Output Files

All outputs are organized under the specified output directory:

```
OUTDIR/
├── qc/
│   ├── per_read_metrics.csv       # Per-read metrics (CSV)
│   ├── per_read_metrics.parquet   # Per-read metrics (Parquet)
│   └── summary.json               # Aggregated summary statistics
├── trim/
│   ├── trimmed.fastq.gz          # Trimmed sequences (FASTQ, gzipped)
│   └── trimmed.fasta.gz          # Trimmed sequences (FASTA, optional)
├── base_calls/                    # Generated with --ambiguous-calling flag
│   ├── all_base_calls.csv        # Per-base annotations (SPR, SNR, IUPAC codes)
│   └── all_base_calls.parquet    # Same data, Parquet format
├── plots/                         # Generated with --plots flag
│   ├── <sample>_trim.png         # Individual sequence quality plots
│   ├── sequences_overview.png    # Multi-sequence grid view
│   ├── summary_histograms.png    # QC metrics distributions
│   └── length_comparison.png     # Before/after trimming comparison
└── logs/
    └── run.log                    # Detailed run log
```

## Visualization Plots

When using the `--plots` flag, the tool generates graphical visualizations to help understand quality and trimming:

### Individual Sequence Plots
- Quality profile with trim regions highlighted
- Green shaded area: kept sequence
- Red shaded areas: discarded ends
- Blue dashed line: quality threshold
- Shows trim coordinates and statistics

### Sequences Overview
- Grid layout showing multiple sequences
- Quick visual comparison of trimming across samples
- Color-coded kept/discarded regions

### Summary Histograms
Six distribution plots showing:
1. Raw sequence lengths
2. Trimmed sequence lengths
3. Mean quality scores
4. Q20 percentages
5. GC content
6. Expected errors

### Length Comparison
- Bar chart: raw vs trimmed lengths per sample
- Scatter plot: correlation between raw and trimmed lengths
- Shows overall retention statistics

**Example with plots:**
```bash
sangerqc all data/ -o results --plots --recursive
# View plots in: results/plots/
```

### Read IDs in Trimmed Outputs

Trimmed sequences are labeled with the format:
```
sample_id/trim:start-end
```

Example: `sample_001/trim:15-542`

## Examples

### Example 1: Full QC and Trimming Pipeline

```bash
sangerqc all seq_data/Y19/ -o results_y19 \
  --qthreshold 20 \
  --min-length 100 \
  --method mott \
  --recursive
```

**Output**:
- QC metrics in `results_y19/qc/`
- Trimmed FASTQ in `results_y19/trim/trimmed.fastq.gz`
- Summary stats in `results_y19/qc/summary.json`

### Example 2: High-Stringency QC Only

```bash
sangerqc qc data/sanger_reads/ -o qc_output \
  --qthreshold 30 \
  --min-length 200 \
  --recursive \
  --verbose
```

### Example 3: Custom Trimmed Output Paths

```bash
sangerqc trim data/ -o results \
  --method ends \
  --out-fastq final_trimmed.fastq.gz \
  --out-fasta final_trimmed.fasta.gz
```

### Example 4: Process Specific Files

```bash
sangerqc all sample1.ab1 sample2.phd.1 sample3.ab1 \
  -o mixed_samples \
  --qthreshold 25
```

### Example 5: Ambiguous Base Calling for Heterozygous SNP Detection

```bash
# Diploid samples with IUPAC ambiguity codes
sangerqc all diploid_samples/ -o snp_results \
  --ambiguous-calling \
  --mixed-context \
  --recursive

# Output includes:
# - FASTQ with IUPAC codes (R, Y, S, W, K, M)
# - base_calls/all_base_calls.csv with SPR, SNR, allele fractions
```

### Example 6: Plasmid Clone Purity Check

```bash
# Check for contamination in clonal samples
sangerqc qc plasmid_clones/ -o purity_check \
  --ambiguous-calling \
  --clonal-context \
  --spr-noise 0.15 \
  --recursive

# Flag clones with mixed signals for re-streaking
```

## Interpreting Results

### Expected Errors

The expected error count is calculated as:

```
EE = Σ 10^(-Q/10) for all bases
```

Lower values indicate higher overall quality. A read with EE < 1 is generally high quality.

### Passed vs Failed Reads

Reads are flagged as:
- **Passed** (`passed_minlen = yes`): Trimmed length ≥ `--min-length`
- **Failed** (`passed_minlen = no`): Trimmed length < `--min-length`

**Note**: Trimming still occurs for failed reads; they're simply flagged for filtering.

### Quality Thresholds

Common quality thresholds:
- **Q20**: 99% base call accuracy (1 error per 100 bases)
- **Q30**: 99.9% base call accuracy (1 error per 1,000 bases)

## Limitations

1. **AB1 files without quality scores**: Files lacking embedded quality annotations are skipped with a warning
2. **No base-calling**: This tool does not perform base-calling; it requires pre-called sequences with quality scores
3. **Single-read files**: Each `.phd.1` file is assumed to contain one read (the first record is used)
4. **Memory usage**: Very large files are processed one-by-one, but aggregated statistics are held in memory

## Development

### Running Tests

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run with coverage
pytest --cov=sanger_qc_trim --cov-report=html
```

### Code Quality

```bash
# Format code
black sanger_qc_trim tests

# Lint code
ruff check sanger_qc_trim tests
```

## Troubleshooting

### "No valid input files found"

- Ensure files have correct extensions: `.ab1` or `.phd.1`
- Use `--recursive` flag for nested directories
- Check file paths are correct

### "No quality scores found in file"

- Some AB1 files may lack embedded quality data
- These files are automatically skipped with warnings logged

### "No valid reads were processed"

- All input files failed to parse (check logs in `OUTDIR/logs/run.log`)
- Verify files are valid Sanger sequence files

## Citation

If you use this tool in your research, please cite:

```
Sanger QC & Trim: A quality control and trimming tool for Sanger sequencing reads
JGI Science, 2025
```

## License

MIT License - see LICENSE file for details

## Support

For bug reports and feature requests, please open an issue at the project repository.

## Acknowledgments

Built with:
- [Biopython](https://biopython.org/) for sequence parsing
- [Typer](https://typer.tiangolo.com/) for CLI
- [Pandas](https://pandas.pydata.org/) for data handling
