# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2025-10-15

### Added

#### Ambiguous Base Calling (Major Feature)
- **IUPAC ambiguity code support**: Detect and report mixed signals using standard codes (R, Y, S, W, K, M)
- **Peak intensity extraction**: Analyze raw trace data from AB1 chromatogram files
- **SPR calculation**: Secondary-to-Primary Ratio for quantifying mixture strength
- **SNR calculation**: Signal-to-Noise Ratio for assessing signal quality
- **Six-rule calling criteria**: Hierarchical decision tree for context-aware base calling
  - Rule 1: Low quality/SNR → N
  - Rule 2: Clear single base → Primary base
  - Rule 3: Minor mixture → Context-dependent (clonal vs mixed)
  - Rule 4: Heterozygous (balanced) → IUPAC code
  - Rule 5: Unbalanced mixture → Context-dependent
  - Rule 6: Very high SPR → N (unclear primary)
- **Context-aware calling**: Different strategies for clonal vs. mixed samples
- **Allele fraction estimation**: Quantify relative abundance (H1/(H1+H2))
- **Per-base annotations**: Detailed CSV/Parquet output with intensities, SPR, SNR, flags

#### New CLI Options
- `--ambiguous-calling` - Enable ambiguous base calling with IUPAC codes
- `--clonal-context` / `--mixed-context` - Specify sample type (default: clonal)
- `--spr-noise FLOAT` - Max SPR for noise threshold (default: 0.20)
- `--spr-het-low FLOAT` - Lower SPR bound for heterozygous calls (default: 0.33)
- `--spr-het-high FLOAT` - Upper SPR bound for heterozygous calls (default: 0.67)

#### New Output Files
- `base_calls/all_base_calls.csv` - Per-base annotations (when --ambiguous-calling enabled)
- `base_calls/all_base_calls.parquet` - Same data in Parquet format
- FASTQ files now include IUPAC codes when ambiguous calling is enabled

#### New Module
- `sanger_qc_trim/ambiguous_calling.py` (675 lines)
  - `AmbiguousBaseCaller` class with full implementation
  - `PeakIntensityExtractor` class for AB1 trace analysis
  - `AmbiguousCallingConfig` dataclass for configuration
  - `BaseCall` dataclass for per-base results
  - Factory function `create_ambiguous_caller()` for easy setup

#### Testing
- 30 new unit tests for ambiguous calling module
- Total test count: 69 tests (all passing)
- Coverage of IUPAC mapping, calling criteria, peak extraction, fallback calling

#### Documentation
- **AMBIGUOUS_CALLING.md**: 400+ line comprehensive guide
  - Detailed explanation of all features
  - IUPAC code reference table
  - Six-rule decision tree documentation
  - Usage examples for different scenarios
  - Configuration guide with conservative/sensitive presets
  - Output format specifications
  - Use cases (SNP detection, clone purity, contamination detection)
  - Troubleshooting guide
- Updated **README.md** with ambiguous calling features
- Updated **QUICKSTART.md** with new examples
- Updated **GET_STARTED.txt** with new options
- Updated **FINAL_SUMMARY.txt** with new metrics

### Enhanced
- `writers.py` - Added `write_base_call_annotations()` and `write_all_base_call_annotations()`
- `cli.py` - Integrated ambiguous calling into all three commands (qc, trim, all)
- All commands now support ambiguous calling as optional enhancement

### Use Cases
1. **Heterozygous SNP detection** in diploid samples
2. **Clone purity assessment** for plasmid sequencing
3. **Contamination detection** in clonal cultures
4. **Viral quasispecies analysis** with minority variant detection
5. **Mixed infection detection** in clinical samples

### Technical Notes
- Ambiguous calling works only with AB1 files (requires raw trace data)
- PHD.1 files fall back to quality-based calling when --ambiguous-calling is enabled
- Configurable thresholds allow tuning for different applications
- Per-base annotations include 12 metrics per position

### Dependencies
No new dependencies required (uses existing biopython, numpy, pandas stack)

---

## [0.1.0] - 2025-10-01

### Added

#### Core Functionality
- Initial release of Sanger QC & Trim tool
- Support for `.ab1` (ABI chromatogram) and `.phd.1` (Phred) file formats
- Automatic format detection based on file extension
- Recursive directory scanning

#### Trimming Algorithms
- **Mott/Kadane algorithm**: Maximum-score subsequence approach for optimal quality-based trimming
- **Ends algorithm**: Conservative 5' and 3' hard clipping to high-quality bases
- Configurable quality threshold (default: Q20)
- Minimum length filtering with pass/fail flagging

#### QC Metrics (Per-Read)
- Sample ID and source file tracking
- Raw and trimmed sequence lengths
- Mean and median quality scores
- Percentage of bases at Q20 and Q30 thresholds
- GC content calculation (excluding N bases)
- N base counting
- Expected errors calculation (sum of error probabilities)
- Longest high-quality stretch detection
- Trim coordinates (0-based start, end)
- Pass/fail status vs minimum length threshold

#### QC Metrics (Aggregated)
- Total read count
- Mean and median raw/trimmed lengths
- Mean quality metrics across all reads
- Pass/fail counts and percentages

#### Output Formats
- **CSV**: Per-read metrics table
- **Parquet**: Per-read metrics (efficient binary format)
- **JSON**: Aggregated summary statistics
- **FASTQ**: Trimmed sequences with qualities (gzipped)
- **FASTA**: Trimmed sequences without qualities (optional, gzipped)

#### CLI Interface
- `sangerqc qc` - QC analysis only
- `sangerqc trim` - Trimming only
- `sangerqc all` - Combined QC and trimming
- Common options: `--qthreshold`, `--method`, `--min-length`, `--recursive`
- Trim-specific: `--out-fastq`, `--out-fasta`
- Logging options: `--verbose`, `--quiet`

#### User Experience
- Progress bars with tqdm
- Clear console output with summary statistics
- Detailed logging to file
- Graceful error handling (skip bad files, continue processing)
- Non-zero exit codes on fatal errors
- Read IDs with trim coordinate annotations

#### Developer Features
- Type hints throughout codebase
- Comprehensive docstrings
- Modular architecture (separate modules for IO, QC, trimming, writing)
- Clean separation of concerns

#### Testing
- 39 unit tests with pytest
- Test coverage for all core algorithms
- Synthetic test data
- Fixtures for common test scenarios
- Verification script for spec examples

#### Documentation
- Complete README.md with installation, usage, examples
- QUICKSTART.md for rapid onboarding
- PROJECT_SUMMARY.md with project overview
- Inline code documentation
- Troubleshooting guide
- Method descriptions with examples

#### Packaging
- Standard pyproject.toml configuration
- Minimal dependency footprint
- pip-installable
- Entry point: `sangerqc` command
- MIT License

### Tested
- All unit tests passing (39/39)
- Verified against specification examples
- Tested on real Sanger sequencing data (`seq_data/Y19/`)
- Tested on Linux platform (Python 3.10.13)

### Dependencies
- biopython ≥1.79
- pandas ≥1.3.0
- numpy ≥1.21.0
- pyarrow ≥6.0.0
- typer ≥0.9.0
- tqdm ≥4.62.0

### Known Limitations
- AB1 files without embedded quality scores are skipped
- No base-calling functionality (requires pre-called sequences)
- Single-read PHD.1 files (first record used if multiple)
- Aggregated statistics held in memory (not an issue for typical datasets)

### Platform Support
- ✅ Linux (tested)
- ⚠️ macOS (expected to work, not tested)
- ⚠️ Windows (expected to work, not tested)

---

## Future Releases

### Planned Features
- Parallel processing with `--threads` option
- Quality histogram plots with `--plots` option
- Separate output for failed reads with `--keep-failed` option
- Additional visualization options
- Extended platform testing (macOS, Windows)

### Potential Enhancements
- Support for additional Sanger formats
- Batch processing optimizations
- Integration with common bioinformatics pipelines
- Docker containerization
- Conda package
