# Ambiguous Base Calling for Sanger Sequencing

## Overview

The Sanger QC & Trim tool now includes advanced ambiguous base calling capabilities that can detect and report mixed signals in Sanger chromatograms. This feature is particularly useful for:

- **Heterozygous SNP detection** in diploid samples
- **Mixed template detection** (contamination, co-infections)
- **Clone purity assessment** in plasmid sequencing
- **Mutation detection** with allele frequency estimation

## Key Features

- **IUPAC ambiguity codes**: Automatically detects and reports mixed signals using standard IUPAC codes (R, Y, S, W, K, M)
- **Configurable thresholds**: Tune sensitivity for different sample types (clonal vs. mixed)
- **Peak intensity analysis**: Extracts raw trace data from AB1 files for accurate mixture detection
- **Quality metrics**: Calculates SPR (Secondary-to-Primary Ratio), SNR (Signal-to-Noise Ratio)
- **Per-base annotations**: Detailed output with intensity, quality, and mixture information
- **Context-aware calling**: Different strategies for clonal vs. mixed samples

## IUPAC Ambiguity Codes

| Code | Bases | Meaning |
|------|-------|---------|
| R    | A/G   | puRine |
| Y    | C/T   | pYrimidine |
| S    | G/C   | Strong |
| W    | A/T   | Weak |
| K    | G/T   | Keto |
| M    | A/C   | aMino |
| N    | any   | uNknown |

## How It Works

### 1. Peak Intensity Extraction

For each base position in an AB1 file, the tool:
- Extracts raw trace intensities from all four channels (A, C, G, T)
- Integrates peak areas in a configurable window around each base call
- Identifies primary (H1) and secondary (H2) signal intensities

### 2. Metric Calculation

**Secondary-to-Primary Ratio (SPR)**: `SPR = H2 / H1`
- Low SPR (< 0.20): Secondary signal is likely noise
- Medium SPR (0.33 - 0.67): Heterozygous mixture
- High SPR (> 0.67): Unbalanced mixture or unclear primary

**Signal-to-Noise Ratio (SNR)**: `SNR = H1 / baseline_noise`
- Measures signal quality relative to background noise
- Low SNR (< 4): Unreliable signal

**Allele Fraction**: `H1 / (H1 + H2)`
- Estimates relative abundance of primary allele
- Useful for quantifying mixture ratios

### 3. Six-Rule Calling Criteria

The tool applies a hierarchical decision tree:

**Rule 1: Too Noisy → Call N**
- If Q < 12 OR SNR < 4
- Emit N (no call)

**Rule 2: Clear Single Base → Call Primary**
- If Q ≥ 30 AND SPR < 0.20
- Call primary base (secondary is noise)

**Rule 3: Minor Mixture (Borderline)**
- If 0.20 ≤ SPR < 0.33 AND Q ≥ 20
- **Clonal context**: Call primary base with `minor_secondary` flag
- **Mixed context**: Call IUPAC code

**Rule 4: Heterozygous → Call IUPAC**
- If 0.33 ≤ SPR ≤ 0.67 AND Q ≥ 20 AND SNR ≥ 5
- Call IUPAC code with `heterozygous` flag
- Report allele fraction

**Rule 5: Unbalanced Mixture**
- If 0.67 < SPR < 0.85 AND Q ≥ 20
- **Mixed samples**: Call IUPAC with `unbalanced_mixture` flag
- **Clonal samples**: Call primary with `uncertain_mixture` flag

**Rule 6: Unclear Primary → Call N**
- If SPR ≥ 0.85
- Emit N with `unclear_primary` flag

## Usage

### Basic Usage

Enable ambiguous calling with the `--ambiguous-calling` flag:

```bash
# QC + ambiguous calling
sangerqc qc seq_data/Y19/ -o results --ambiguous-calling --recursive

# Trim + ambiguous calling
sangerqc trim seq_data/Y19/ -o results --ambiguous-calling --recursive

# Full pipeline with ambiguous calling
sangerqc all seq_data/Y19/ -o results --ambiguous-calling --recursive
```

### Advanced Configuration

#### Clonal vs. Mixed Context

For **clonal samples** (plasmids, single colonies):
```bash
sangerqc all seq_data/ -o results --ambiguous-calling --clonal-context
```

For **mixed samples** (diploid, environmental, co-infections):
```bash
sangerqc all seq_data/ -o results --ambiguous-calling --mixed-context
```

#### Custom Thresholds

Adjust sensitivity thresholds:

```bash
sangerqc all seq_data/ -o results \
  --ambiguous-calling \
  --spr-noise 0.15 \          # Stricter noise threshold (default: 0.20)
  --spr-het-low 0.30 \        # Lower heterozygous bound (default: 0.33)
  --spr-het-high 0.70         # Upper heterozygous bound (default: 0.67)
```

**Common presets:**

Conservative (fewer ambiguous calls):
```bash
--spr-noise 0.15 --spr-het-low 0.40 --spr-het-high 0.60
```

Sensitive (more ambiguous calls):
```bash
--spr-noise 0.25 --spr-het-low 0.25 --spr-het-high 0.75
```

## Output Files

When ambiguous calling is enabled, additional output files are generated:

```
OUTDIR/
├── base_calls/
│   ├── all_base_calls.csv       # Combined per-base annotations
│   └── all_base_calls.parquet   # Same data, Parquet format
├── qc/
│   ├── per_read_metrics.csv     # Standard QC metrics (updated sequences)
│   └── summary.json
└── trim/
    └── trimmed.fastq.gz         # Trimmed sequences with IUPAC codes
```

### Base Call Annotations Format

The `all_base_calls.csv` file contains detailed per-base information:

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| position | 0-based position in sequence |
| called_base | Final base call (A/C/G/T/IUPAC/N) |
| primary_base | Primary base by intensity |
| secondary_base | Secondary base by intensity |
| primary_intensity | H1 (integrated peak area) |
| secondary_intensity | H2 (integrated peak area) |
| spr | Secondary-to-Primary Ratio |
| snr | Signal-to-Noise Ratio |
| quality | Phred quality score |
| call_mode | single, ambiguous, or N |
| allele_fraction | H1 / (H1 + H2) |
| flags | Comma-separated flags (e.g., heterozygous, minor_secondary) |

### Example Output

```csv
sample_id,position,called_base,primary_base,secondary_base,primary_intensity,secondary_intensity,spr,snr,quality,call_mode,allele_fraction,flags
sample_001,0,A,A,G,1234.5,123.4,0.10,12.5,35,single,0.9091,
sample_001,1,R,A,G,1000.0,500.0,0.50,10.2,28,ambiguous,0.6667,heterozygous
sample_001,2,C,C,T,1500.0,200.0,0.13,11.8,32,single,0.8824,
sample_001,3,N,A,G,200.0,180.0,0.90,3.2,15,N,0.5263,unclear_primary
```

## Interpreting Results

### Flag Meanings

- **heterozygous**: Balanced mixture (SPR 0.33-0.67), likely true heterozygous SNP
- **minor_secondary**: Weak secondary signal (SPR 0.20-0.33), flagged but called as primary
- **unbalanced_mixture**: Unbalanced mixture (SPR 0.67-0.85), one allele dominant
- **uncertain_mixture**: Unclear mixture in clonal context
- **unclear_primary**: Very high SPR (>0.85), can't determine which is primary
- **low_quality**: Quality or SNR below threshold
- **no_trace_data**: Trace data unavailable (PHD.1 files or AB1 parsing failure)

### Quality Assessment

**High-confidence calls:**
- `call_mode=single`, no flags, Q ≥ 30, SNR ≥ 10

**Likely heterozygous:**
- `call_mode=ambiguous`, `heterozygous` flag, SPR 0.4-0.6, Q ≥ 25

**Potential contamination:**
- `call_mode=ambiguous`, `unbalanced_mixture` flag in clonal samples
- Multiple positions with SPR 0.2-0.4

**Low-quality positions:**
- `call_mode=N`, `low_quality` or `unclear_primary` flags
- Should be excluded from downstream analysis

## Use Cases

### 1. Diploid SNP Calling

```bash
# Mixed context for diploid samples
sangerqc all diploid_samples/ -o snp_results \
  --ambiguous-calling \
  --mixed-context \
  --spr-het-low 0.30 \
  --spr-het-high 0.70
```

**Analysis:**
- Look for positions with `heterozygous` flag and SPR 0.4-0.6
- Cross-reference forward and reverse reads
- Filter by Q ≥ 25 and SNR ≥ 5

### 2. Plasmid Clone Purity

```bash
# Clonal context for plasmid sequencing
sangerqc qc plasmid_clones/ -o purity_check \
  --ambiguous-calling \
  --clonal-context \
  --spr-noise 0.15
```

**Analysis:**
- Flag clones with multiple `minor_secondary` or `uncertain_mixture` positions
- Calculate purity score: % positions without mixture flags
- Re-streak clones with < 95% purity

### 3. Viral Quasispecies Analysis

```bash
# Sensitive detection for viral populations
sangerqc all viral_samples/ -o quasispecies \
  --ambiguous-calling \
  --mixed-context \
  --spr-noise 0.25 \
  --spr-het-low 0.25
```

**Analysis:**
- Quantify minority variants using allele_fraction
- Track positions with consistent mixtures across samples
- Estimate within-host diversity

### 4. Contamination Detection

```bash
# Detect cross-contamination in clonal samples
sangerqc qc all_samples/ -o contamination_check \
  --ambiguous-calling \
  --clonal-context
```

**Analysis:**
- Flag samples with > 5 ambiguous positions
- Check for shared IUPAC patterns across samples (barcode swaps)
- Validate with independent sequencing

## Limitations

1. **AB1 files only**: Peak intensity extraction requires raw chromatogram data. PHD.1 files fall back to quality-based calling.

2. **Single-position analysis**: Does not consider linkage or haplotype information.

3. **Homopolymer regions**: May produce false positives due to peak spacing issues.

4. **Trace quality**: Requires clean chromatograms; noisy data may produce unreliable SPR/SNR.

5. **Three-way mixtures**: Current implementation focuses on two-base mixtures (does not call B/D/H/V codes).

## Validation

The ambiguous calling module includes 30 unit tests covering:
- IUPAC code mapping
- All 6 calling criteria rules
- SPR/SNR calculations
- Fallback quality-based calling
- Configuration handling

Run tests:
```bash
pytest tests/test_ambiguous_calling.py -v
```

## Troubleshooting

### Issue: "Ambiguous calling failed" warnings

**Cause**: AB1 file missing trace data or corrupted

**Solution**: Check if files open in chromatogram viewer. Some AB1 files lack raw trace data.

### Issue: Too many N calls

**Cause**: Low SNR or poor quality chromatograms

**Solution**:
- Lower `--snr-min` threshold (default: 4.0)
- Check sequencing quality
- Try different trace analysis settings

### Issue: Expected heterozygous sites called as single base

**Cause**: SPR thresholds too strict

**Solution**:
```bash
--spr-het-low 0.25 --spr-het-high 0.75  # More sensitive
```

### Issue: Too many false positive mixtures

**Cause**: Noisy traces or thresholds too permissive

**Solution**:
```bash
--spr-noise 0.15 --spr-het-low 0.40  # More conservative
```

## Technical Details

### Peak Integration Window

By default, peak intensities are integrated over a ±3 base window around each peak location. This can be adjusted in the configuration:

```python
from sanger_qc_trim.ambiguous_calling import create_ambiguous_caller

caller = create_ambiguous_caller(peak_window=5)  # Larger window
```

### SNR Calculation

Signal-to-noise ratio is estimated using:
1. Identify peak position from trace index
2. Sample baseline regions ±20 positions from peak
3. Calculate median of baseline samples as noise estimate
4. SNR = H1 / noise

### Fallback Calling

When trace data is unavailable (PHD.1 files or parsing errors), the tool falls back to quality-based calling:
- Q < 12: Call N
- Q ≥ 30: Call original base (high confidence)
- 12 ≤ Q < 30: Call original base with `moderate_quality` flag

All fallback calls include `no_trace_data` flag.

## References

- IUPAC Nucleotide Code: https://www.bioinformatics.org/sms/iupac.html
- Sanger Sequencing Quality: Ewing & Green (1998) Genome Research
- AB1 File Format: Applied Biosystems ABI PRISM specifications

## Citation

If you use the ambiguous calling feature in your research, please cite:

```
Sanger QC & Trim: Production-ready bioinformatics tool for quality control
and trimming of Sanger sequencing reads with ambiguous base calling support.
```
