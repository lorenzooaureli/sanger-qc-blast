# DEVELOPER_GUIDE.md

This file provides guidance for developers working with code in this repository.

## Project Overview

Sanger QC & Trim is a production-ready bioinformatics tool for quality control and trimming of Sanger sequencing reads from `.ab1` (ABI chromatogram) and `.phd.1` (Phred) files. The tool provides comprehensive QC metrics, intelligent trimming algorithms, and optional visualization plots.

## Development Commands

### Installation & Setup
```bash
# Install in development mode
pip install -e .

# Install with dev dependencies (testing, formatting, linting)
pip install -e ".[dev]"
```

### Testing
```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=sanger_qc_trim --cov-report=html

# Run specific test file
pytest tests/test_trim.py

# Run specific test
pytest tests/test_trim.py::TestTrimMott::test_basic_trimming

# Verify algorithms against spec
python verify_algorithms.py
```

### Code Quality
```bash
# Format code (100 char line length)
black sanger_qc_trim tests

# Lint code
ruff check sanger_qc_trim tests
```

### Running the Tool
```bash
# Full QC + trimming pipeline
sangerqc all seq_data/Y19/ -o results --recursive

# QC only
sangerqc qc seq_data/Y19/ -o results --recursive

# Trimming only
sangerqc trim seq_data/Y19/ -o results --recursive

# With plots (adds --plots flag)
sangerqc all seq_data/Y19/ -o results --plots --recursive

# Test on sample data
sangerqc all seq_data/Y19/ -o test_output --recursive
```

## Architecture

### Module Structure

The codebase follows a clean separation of concerns with 6 core modules:

```
sanger_qc_trim/
├── cli.py          # Typer-based CLI with 3 subcommands (qc, trim, all)
├── io_utils.py     # File discovery, deduplication, format detection, Biopython parsing
├── trim.py         # Trimming algorithms (Mott/Kadane, Ends)
├── qc.py           # QC metrics computation (per-read + aggregated)
├── writers.py      # Output writers (CSV, Parquet, JSON, FASTQ, FASTA)
└── plots.py        # Matplotlib + Plotly visualization (4 plot types, static + interactive)
```

### Data Flow

1. **Input Processing** (`io_utils.py`):
   - `discover_files()` recursively finds .ab1 and .phd.1 files
   - **Automatic deduplication**: Files with same sample_id are deduplicated, preferring .ab1 over .phd.1
   - `detect_format()` auto-detects file type by extension
   - `parse_sequence_file()` uses Biopython SeqIO to extract sequences and quality scores
   - Returns `(sequence_string, quality_list)` tuples

2. **Trimming** (`trim.py`):
   - Two algorithms implemented as pure functions:
     - `trim_mott()`: Kadane-style maximum subsequence (weights = Q - threshold)
     - `trim_ends()`: Conservative 5'/3' clipping to first/last high-quality base
   - Returns 0-based `(start, end)` indices (half-open interval)
   - If all qualities below threshold → `(0, 0)` (empty trim)

3. **QC Metrics** (`qc.py`):
   - `compute_qc_metrics()`: Computes 16 per-read metrics including:
     - Quality stats: mean_q, median_q, pct_q20, pct_q30
     - Content: gc_percent, n_count
     - Error metrics: expected_errors (Σ10^(-Q/10))
     - Trimming: trim_start, trim_end, trimmed_length, passed_minlen
   - `compute_summary_stats()`: Aggregates 13 summary metrics across all reads

4. **Output Generation** (`writers.py`):
   - CSV and Parquet for per-read metrics
   - JSON for summary statistics
   - Gzipped FASTQ/FASTA for trimmed sequences
   - Read IDs formatted as: `sample_id/trim:start-end`

5. **Visualization** (`plots.py`, optional with `--plots`):
   - 4 plot types generated:
     - **Individual quality plots**: Both static (PNG) and interactive (HTML) versions
       - Matplotlib: `<sample>_trim.png` (green=kept, red=discarded)
       - Plotly: `<sample>_trim_interactive.html` (interactive hover, zoom, pan)
     - Multi-sequence overview grid (2x5 layout, PNG)
     - Summary histograms (6 metrics, PNG)
     - Length comparison (bar + scatter plots, PNG)
   - Saves to `OUTDIR/plots/`
   - Handles both single and multiple sequence inputs correctly

### CLI Architecture (`cli.py`)

Three Typer subcommands share common processing logic:
- `qc`: Computes metrics only (no FASTQ output)
- `trim`: Generates trimmed sequences only (no QC CSV/Parquet)
- `all`: Combines both QC and trimming in single pass

All commands:
- Use tqdm progress bars for file processing
- Log to both console and `OUTDIR/logs/run.log`
- Exit with code 1 on fatal errors, 0 on success (even if some files skipped)
- Support `--plots` flag to generate visualizations

### Key Design Patterns

**Error Handling Philosophy**:
- Skip individual bad files with warnings (continue processing)
- Only fail on zero valid reads processed
- Log skipped files to help users identify issues

**Streaming Processing**:
- Files processed one at a time (no batching)
- Metrics accumulated in lists (memory efficient for typical datasets)
- Large files handled via Biopython's streaming parsers

**Quality Score Handling**:
- All quality scores are Phred scores (0-60 typical range)
- Extracted from `record.letter_annotations["phred_quality"]`
- Files without quality annotations are skipped with warnings

## Testing Strategy

### Test Organization
- `tests/test_trim.py`: 23 tests covering both trimming algorithms
- `tests/test_qc.py`: 16 tests covering QC metrics computation
- `tests/conftest.py`: Pytest fixtures for test data

### Critical Test Cases
- Boundary conditions: empty sequences, all low/high quality
- Spec examples: Verify against documented examples in README
- Algorithm correctness: Mott vs Ends produce different expected results

### Verification Script
`verify_algorithms.py` validates core algorithms against specification:
- Mott trimming: `[10,10,30,30,30,10]` @ Q20 → `(2,5)` ✓
- Ends trimming: `[15,25,25,15]` @ Q20 → `(1,3)` ✓
- Q20 percentage: `[20,20,10,30]` → `0.75` ✓

## File Format Notes

### AB1 Files
- Binary format from ABI sequencers
- Parsed with `SeqIO.read(path, "abi")`
- Quality scores in `record.letter_annotations["phred_quality"]`
- Some files may lack quality data → skip with warning

### PHD.1 Files
- Plain-text Phred output
- Parsed with `SeqIO.parse(path, "phd")` (yields iterator)
- Usually one record per file (first record used)
- Quality scores always present in valid files

## Output Structure

```
OUTDIR/
├── qc/
│   ├── per_read_metrics.csv      # Pandas DataFrame → CSV
│   ├── per_read_metrics.parquet  # Same data, Parquet format
│   └── summary.json              # Aggregated stats
├── trim/
│   └── trimmed.fastq.gz          # Gzipped FASTQ (qualities preserved)
├── plots/                         # With --plots flag
│   ├── <sample>_trim.png         # Individual static plots (max 10)
│   ├── <sample>_trim_interactive.html  # Individual interactive plots (max 10)
│   ├── sequences_overview.png    # Multi-sequence grid
│   ├── summary_histograms.png    # 6 metric distributions
│   └── length_comparison.png     # Before/after analysis
└── logs/
    └── run.log                    # Detailed execution log
```

## Common Development Scenarios

### Adding a New QC Metric

1. Add computation to `qc.py::compute_qc_metrics()`
2. Update return dictionary with new metric
3. Add to `compute_summary_stats()` if aggregation needed
4. Update tests in `tests/test_qc.py`
5. Document in README.md metrics table

### Adding a New Trimming Algorithm

1. Implement as pure function in `trim.py`:
   ```python
   def trim_mymethod(quals: list[int], threshold: int) -> Tuple[int, int]:
       # Return (start, end) 0-based indices
   ```
2. Add to `apply_trim()` elif chain
3. Update CLI `--method` choices in `cli.py`
4. Add tests to `tests/test_trim.py`
5. Document in README.md "Trimming Methods" section

### Modifying Plot Appearance

All plot styling in `plots.py`:
- Figure sizes: `figsize=(width, height)`
- DPI: Currently 150 (adjustable in `plt.savefig(dpi=...)`)
- Colors: Green (#008000), Red (#FF0000), Blue (#0000FF)
- Max individual plots: 10 (configurable in `cli.py`)
- Interactive plots use Plotly with hover tooltips and zoom/pan controls
- **Important**: `plot_multiple_sequences()` handles single sequence case with special axes reshaping logic

## Important Implementation Details

### File Deduplication
- `discover_files()` automatically deduplicates files with the same sample_id
- Sample ID extracted from filename (without extension)
- **Preference order**: .ab1 files preferred over .phd.1 files
- Rationale: Both formats contain identical sequence data, .ab1 is native format
- Logs show deduplication statistics: `"Deduplicated N files to M unique samples"`
- Warning logged if duplicate sample_ids found with same format

### Interactive Plotting
- Individual quality plots generated in both formats:
  - Static PNG via matplotlib (`plot_sequence_trim()`)
  - Interactive HTML via Plotly (`plot_sequence_trim_interactive()`)
- Plotly plots support:
  - Hover tooltips showing position and quality score
  - Zoom and pan controls
  - Exportable as PNG from browser
- Both versions show identical visualizations (kept/discarded regions)
- HTML files are ~4.5 MB each (include Plotly.js library)

### Trimming Coordinate Convention
- All coordinates are **0-based**
- Indices are **half-open intervals**: `[start, end)`
- Example: `trim_start=2, trim_end=5` means bases at positions 2, 3, 4
- Empty trim: `(0, 0)` means trimmed sequence is empty string

### Quality Threshold Logic
- Mott: Uses weights `Q - threshold` (positive = good, negative = bad)
- Ends: First/last base with `Q >= threshold`
- Default threshold: Q20 (99% accuracy)
- Typical range: Q15-Q30

### Parallel Processing
Currently NOT implemented (single-process only). If adding:
- Use multiprocessing.Pool for file-level parallelism
- Avoid threading (GIL limitations with Biopython)
- Aggregate results from workers before writing

## Documentation Files

- `README.md`: Complete user documentation
- `QUICKSTART.md`: Quick reference guide
- `INSTALL.md`: Installation details
- `PLOTS_GUIDE.md`: Comprehensive plotting documentation
- `CHANGELOG.md`: Version history
- `PROJECT_SUMMARY.md`: Project overview and metrics

## Dependencies

Core (required):
- biopython: SeqIO parsers for .ab1 and .phd.1
- pandas: DataFrame operations for metrics
- numpy: Numerical computations
- pyarrow: Parquet file format support
- typer: CLI framework
- tqdm: Progress bars
- matplotlib: Static plotting (for --plots feature)
- plotly: Interactive HTML plots (for --plots feature)

Dev (optional):
- pytest: Testing framework
- pytest-cov: Coverage reporting
- black: Code formatting (100 char lines)
- ruff: Linting
