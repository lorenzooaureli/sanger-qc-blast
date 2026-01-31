#!/bin/bash
# Example script demonstrating all sangerqc commands

set -e

echo "============================================"
echo "Sanger QC & Trim - Example Commands"
echo "============================================"
echo ""

# Example 1: QC only
echo "1. Running QC analysis only..."
sangerqc qc seq_data/Y19/ -o example_qc --qthreshold 20 --min-length 50
echo "   ✓ QC complete. Results in example_qc/"
echo ""

# Example 2: Trim only
echo "2. Running trimming only..."
sangerqc trim seq_data/Y19/ -o example_trim --method mott --out-fastq example_trim.fastq.gz
echo "   ✓ Trimming complete. Results in example_trim/"
echo ""

# Example 3: Both QC and trimming
echo "3. Running QC and trimming together..."
sangerqc all seq_data/Y19/ -o example_all --qthreshold 20 --method mott --min-length 50
echo "   ✓ QC and trimming complete. Results in example_all/"
echo ""

echo "============================================"
echo "All examples completed successfully!"
echo "============================================"
echo ""
echo "Output directories created:"
echo "  - example_qc/     (QC only)"
echo "  - example_trim/   (Trimming only)"
echo "  - example_all/    (QC + Trimming)"
echo ""
echo "To view results:"
echo "  cat example_all/qc/summary.json"
echo "  cat example_all/qc/per_read_metrics.csv"
echo "  zcat example_all/trim/trimmed.fastq.gz | head -8"
