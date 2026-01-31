# File Inventory

## Core Package Files

### Source Code (`sanger_qc_trim/`)
- `__init__.py` - Package initialization, version info
- `cli.py` - Typer-based CLI with qc/trim/all subcommands
- `io_utils.py` - File discovery, format detection, sequence parsing
- `qc.py` - QC metrics computation (per-read and aggregated)
- `trim.py` - Trimming algorithms (Mott/Kadane and Ends)
- `writers.py` - Output writers for CSV, Parquet, JSON, FASTQ, FASTA

### Test Suite (`tests/`)
- `__init__.py` - Test package initialization
- `conftest.py` - Pytest configuration and fixtures
- `test_qc.py` - Unit tests for QC metrics (16 tests)
- `test_trim.py` - Unit tests for trimming algorithms (23 tests)
- `data/` - Directory for test data fixtures

## Configuration Files

- `pyproject.toml` - Python package configuration (PEP 621)
  - Build system setup
  - Package metadata
  - Dependencies
  - Entry points (sangerqc command)
  - Tool configuration (black, ruff)

- `LICENSE` - MIT License

## Documentation Files

- `README.md` - Complete project documentation
  - Installation instructions
  - Usage guide with examples
  - Method descriptions
  - QC metrics explained
  - Output file formats
  - Troubleshooting guide

- `QUICKSTART.md` - Quick start guide
  - Basic usage examples
  - Common options
  - Output interpretation

- `PROJECT_SUMMARY.md` - Project overview
  - Project structure
  - Testing status
  - Acceptance criteria verification
  - Feature list

- `CHANGELOG.md` - Version history and changes

- `FILES.md` - This file (file inventory)

## Utility Scripts

- `verify_algorithms.py` - Verification script
  - Tests algorithms against specification examples
  - Validates Mott, Ends, and QC metrics

- `run_examples.sh` - Example commands script
  - Demonstrates all three subcommands
  - Creates example outputs

## Output Directory Structure

When running the tool, it creates:

```
OUTDIR/
├── qc/
│   ├── per_read_metrics.csv      # Per-read QC metrics
│   ├── per_read_metrics.parquet  # Same data in Parquet format
│   └── summary.json              # Aggregated statistics
├── trim/
│   └── trimmed.fastq.gz         # Trimmed sequences (FASTQ)
│   └── trimmed.fasta.gz         # Trimmed sequences (FASTA, optional)
└── logs/
    └── run.log                   # Detailed run log
```

## File Count Summary

- **Source files**: 6
- **Test files**: 4
- **Configuration**: 2
- **Documentation**: 5
- **Utilities**: 2

**Total**: 19 files (excluding generated files and caches)

## Dependencies Declared

In `pyproject.toml`:
- biopython ≥1.79
- pandas ≥1.3.0
- numpy ≥1.21.0
- pyarrow ≥6.0.0
- typer ≥0.9.0
- tqdm ≥4.62.0

Development dependencies:
- pytest ≥7.0.0
- pytest-cov ≥3.0.0
- black ≥22.0.0
- ruff ≥0.1.0
