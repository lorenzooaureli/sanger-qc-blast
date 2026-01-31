# Quick Start Guide

## Installation

```bash
# Install the tool
pip install .

# Verify installation
sangerqc --help
```

## Basic Usage

### 1. QC + Trimming (Most Common)

Process all `.ab1` and `.phd.1` files in a directory:

```bash
sangerqc all seq_data/Y19/ -o results --recursive
```

### 2. QC Analysis Only

Get quality metrics without trimming:

```bash
sangerqc qc seq_data/Y19/ -o qc_results --recursive
```

### 3. Trimming Only

Trim sequences and output FASTQ:

```bash
sangerqc trim seq_data/Y19/ -o trim_results --out-fastq trimmed.fastq.gz --recursive
```

## Common Options

- `--qthreshold 20` - Quality threshold (default: 20)
- `--method mott` - Trimming method: `mott` or `ends` (default: mott)
- `--min-length 50` - Minimum acceptable length (default: 50)
- `--recursive` - Search directories recursively
- `--verbose` - Show detailed logging
- `--quiet` - Suppress console output (log to file only)

### Ambiguous Calling Options (NEW!)

- `--ambiguous-calling` - Enable IUPAC ambiguity codes
- `--clonal-context` / `--mixed-context` - Sample type (default: clonal)
- `--spr-noise 0.20` - SPR noise threshold
- `--spr-het-low 0.33` - Heterozygous lower bound
- `--spr-het-high 0.67` - Heterozygous upper bound

## Output Files

All results are organized under your output directory:

```
results/
├── qc/
│   ├── per_read_metrics.csv      # Per-read QC metrics
│   ├── per_read_metrics.parquet  # Same data in Parquet format
│   └── summary.json              # Aggregated statistics
├── trim/
│   └── trimmed.fastq.gz         # Trimmed sequences (with IUPAC codes)
├── base_calls/                   # With --ambiguous-calling flag
│   ├── all_base_calls.csv       # Per-base annotations
│   └── all_base_calls.parquet   # Same in Parquet format
└── logs/
    └── run.log                   # Detailed run log
```

## Examples

### High-stringency QC
```bash
sangerqc qc data/ -o high_qc --qthreshold 30 --min-length 100 --recursive
```

### Conservative Trimming
```bash
sangerqc trim data/ -o conservative --method ends --qthreshold 15 --recursive
```

### Custom Output Paths
```bash
sangerqc all data/ -o results \
  --out-fastq my_trimmed.fastq.gz \
  --out-fasta my_trimmed.fasta.gz \
  --recursive
```

### Heterozygous SNP Detection (NEW!)
```bash
sangerqc all diploid_samples/ -o snp_results \
  --ambiguous-calling \
  --mixed-context \
  --recursive
```

### Plasmid Clone Purity Check (NEW!)
```bash
sangerqc qc plasmid_clones/ -o purity_check \
  --ambiguous-calling \
  --clonal-context \
  --spr-noise 0.15 \
  --recursive
```

## Interpreting Results

### Key Metrics

- **mean_q**: Average quality score (higher is better)
- **pct_q20**: Fraction of bases with Q≥20 (0.8 = 80%)
- **expected_errors**: Sum of error probabilities (lower is better)
- **trimmed_length**: Length after trimming
- **passed_minlen**: "yes" if length ≥ min_length

### Quality Thresholds

- **Q20**: 99% accuracy (1 error per 100 bases)
- **Q30**: 99.9% accuracy (1 error per 1,000 bases)

## Troubleshooting

**No files found?**
- Check file extensions: `.ab1` or `.phd.1`
- Use `--recursive` for nested directories

**Files skipped?**
- Check `logs/run.log` for details
- Some AB1 files may lack quality scores

**Need help?**
```bash
sangerqc --help
sangerqc qc --help
sangerqc trim --help
sangerqc all --help
```

## Running Tests

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run with coverage
pytest --cov=sanger_qc_trim
```
