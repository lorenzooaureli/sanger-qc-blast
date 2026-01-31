# Visualization Plots Guide

## Overview

The `--plots` flag generates graphical visualizations showing quality profiles and trimming results. This helps you visually understand where sequences were trimmed and why.

## Enabling Plots

Add `--plots` to any `sangerqc` command:

```bash
# With QC only
sangerqc qc data/ -o results --plots --recursive

# With trimming
sangerqc trim data/ -o results --plots --recursive

# Full pipeline with plots
sangerqc all data/ -o results --plots --recursive
```

## Generated Plots

### 1. Individual Sequence Quality Plots

**Filename**: `<sample_id>_trim.png` (one per sequence, up to 10 sequences)

**What it shows**:
- Black line: Quality scores across the entire sequence
- **Green shaded region**: Kept (trimmed) portion
- **Red shaded regions**: Discarded portions (5' and 3' ends)
- Blue dashed line: Quality threshold used for trimming
- Gray dotted lines: Q20 and Q30 reference markers
- Text box with statistics:
  - Raw length
  - Trimmed length and percentage
  - Trim coordinates [start, end)
  - Mean quality score

**Use this to**:
- See exactly where each sequence was trimmed
- Understand why certain regions were discarded
- Verify trimming makes sense for your data

**Example interpretation**:
- If most of the sequence is green → good quality, minimal trimming needed
- If large red regions → low quality ends were trimmed off
- If mostly red → sequence failed quality thresholds

---

### 2. Sequences Overview (Grid View)

**Filename**: `sequences_overview.png`

**What it shows**:
- Grid layout with multiple sequences (2x5 grid showing up to 10 sequences)
- Each subplot shows one sequence's quality profile
- Same color coding as individual plots:
  - Green: kept regions
  - Red: discarded regions
  - Blue line: quality threshold
- Subplot titles show sample ID and kept length

**Use this to**:
- Compare trimming across multiple samples quickly
- Identify samples with unusual trimming patterns
- Get an overview of dataset quality

**Example interpretation**:
- Consistent green regions → uniform quality across samples
- Variable patterns → quality varies between samples
- Many short sequences → overall low quality or aggressive trimming

---

### 3. Summary Histograms

**Filename**: `summary_histograms.png`

**What it shows**:
Six histogram panels showing distributions of key metrics:

1. **Raw Sequence Length** (blue)
   - Length distribution before trimming
   - Shows original sequence sizes

2. **Trimmed Sequence Length** (green)
   - Length distribution after trimming
   - Compare to raw lengths to see trimming impact

3. **Mean Quality Distribution** (orange)
   - Average quality per sequence
   - Red dashed lines at Q20 and Q30

4. **Q20 Percentage** (purple)
   - Fraction of bases with Q≥20 per sequence
   - Higher is better

5. **GC Content** (teal)
   - GC% distribution across sequences
   - Usually centered around 50% for random sequences

6. **Expected Errors** (crimson)
   - Predicted error count per sequence
   - Lower is better

**Use this to**:
- Understand overall dataset characteristics
- Identify outliers
- Assess data quality at a glance
- Plan downstream filtering thresholds

**Example interpretation**:
- Narrow distributions → consistent quality
- Wide distributions → variable quality
- Bimodal distributions → mixed quality populations

---

### 4. Length Comparison

**Filename**: `length_comparison.png`

**What it shows**:
Two panels comparing raw and trimmed lengths:

**Left Panel - Bar Chart**:
- Blue bars: Raw sequence lengths
- Green bars: Trimmed sequence lengths
- Side-by-side comparison per sample
- X-axis: Sample names
- Y-axis: Length in base pairs

**Right Panel - Scatter Plot**:
- Each point: one sequence
- X-axis: Raw length
- Y-axis: Trimmed length
- Red dashed diagonal: "no trimming" line
- Points below line: sequences that were trimmed
- Text box shows:
  - Overall retention percentage
  - Mean raw length
  - Mean trimmed length

**Use this to**:
- See how much sequence was lost to trimming
- Identify samples with excessive trimming
- Calculate overall data retention
- Assess trimming uniformity

**Example interpretation**:
- Points close to diagonal → minimal trimming
- Points far below diagonal → heavy trimming
- High retention % → good quality data
- Low retention % → poor quality or stringent thresholds

---

## Output Location

All plots are saved in:
```
OUTDIR/plots/
├── Y19_ITS4_trim.png           # Individual quality plot
├── sequences_overview.png      # Multi-sequence grid
├── summary_histograms.png      # Metrics distributions
└── length_comparison.png       # Before/after comparison
```

## Viewing Plots

### Option 1: Local Machine (Recommended)
Copy plots to your local computer for viewing:
```bash
scp -r results/plots/ your_local_machine:/path/to/destination/
```

### Option 2: File Browser
Navigate to the plots directory and open in any image viewer:
```bash
# On Linux with display
eog results/plots/sequences_overview.png

# On macOS
open results/plots/sequences_overview.png
```

### Option 3: Jupyter Notebook
```python
from IPython.display import Image, display
display(Image('results/plots/sequences_overview.png'))
```

### Option 4: Web Browser
If you have a web server running:
```bash
python -m http.server 8000
# Then open: http://localhost:8000/results/plots/
```

## Tips

### For Large Datasets
- Individual sequence plots are limited to first 10 sequences
- All sequences included in overview grid (up to 10)
- Summary histograms include all sequences
- Length comparison includes all sequences

### Customization
To modify plot appearance, edit `sanger_qc_trim/plots.py`:
- Change colors
- Adjust figure sizes
- Modify DPI (currently 150)
- Add additional plots

### Performance
- Plotting adds minimal time (<5 seconds for typical datasets)
- Plot files are ~100-250 KB each
- Matplotlib required (installed automatically)

## Interpreting Trimming Quality

### Good Trimming
- Most sequence retained (green regions large)
- Consistent patterns across samples
- High overall retention percentage (>80%)
- Trim regions make sense (low quality at ends)

### Problematic Trimming
- Very short trimmed sequences
- Inconsistent trimming across samples
- Low retention (<50%)
- High quality regions discarded

### Troubleshooting Bad Trimming

If plots show excessive trimming:

1. **Lower quality threshold**:
   ```bash
   sangerqc all data/ -o results --qthreshold 15 --plots
   ```

2. **Try different method**:
   ```bash
   sangerqc all data/ -o results --method ends --plots
   ```

3. **Check raw data quality**:
   - Look at summary histograms
   - If mean quality is low, data may need resequencing

4. **Adjust minimum length**:
   ```bash
   sangerqc all data/ -o results --min-length 30 --plots
   ```

## Example Workflows

### Workflow 1: Quality Assessment
```bash
# Generate plots to assess data quality
sangerqc qc data/ -o qc_check --plots --recursive

# Review plots in qc_check/plots/
# - Check summary_histograms.png for overall quality
# - Look at sequences_overview.png for patterns
# - Decide on appropriate thresholds
```

### Workflow 2: Optimizing Trimming
```bash
# Try different thresholds and compare
sangerqc all data/ -o trim_q15 --qthreshold 15 --plots --recursive
sangerqc all data/ -o trim_q20 --qthreshold 20 --plots --recursive
sangerqc all data/ -o trim_q25 --qthreshold 25 --plots --recursive

# Compare length_comparison.png across all three
# Choose best balance of quality and retention
```

### Workflow 3: Method Comparison
```bash
# Compare Mott vs Ends trimming
sangerqc all data/ -o mott_trim --method mott --plots --recursive
sangerqc all data/ -o ends_trim --method ends --plots --recursive

# Compare sequences_overview.png
# Mott: optimal quality-based trimming
# Ends: conservative end-clipping
```

## Frequently Asked Questions

**Q: Why are only 10 sequences shown in individual plots?**
A: To keep file sizes reasonable and processing fast. All sequences are included in summary statistics.

**Q: Can I generate plots for specific samples only?**
A: Yes, provide specific files as input:
```bash
sangerqc all sample1.ab1 sample2.ab1 -o results --plots
```

**Q: What if I don't see green regions?**
A: All sequences failed quality thresholds. Try lowering `--qthreshold` or check data quality.

**Q: Can I change plot colors?**
A: Yes, edit `sanger_qc_trim/plots.py` and modify the color parameters.

**Q: Do plots slow down processing?**
A: Minimally (~1-5 seconds for typical datasets). Most time is spent on sequence parsing.

## Summary

The plotting feature provides:
- ✅ Visual validation of trimming decisions
- ✅ Quick quality assessment
- ✅ Dataset-wide statistics
- ✅ Before/after comparisons
- ✅ Easy sharing and reporting

Use `--plots` whenever you want to:
- Understand your data quality
- Optimize trimming parameters
- Generate figures for reports
- Troubleshoot unexpected results
- Share results with collaborators
