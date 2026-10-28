#!/usr/bin/env bash
# Step 4 - Quality control of raw reads with FastQC, summarised with MultiQC.
# Usage (egf-ngs environment active):
#   caffeinate -i bash scripts/04_qc.sh 6
# The argument is the number of threads (default: 6).

set -euo pipefail   # stop at the first error, undefined variable or failed pipe

THREADS="${1:-6}"

# Project root, built from the location of this script (same idea as in step 3)
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
FASTQ_DIR="$PROJECT_DIR/data/raw/fastq"
QC_DIR="$PROJECT_DIR/results/qc"

mkdir -p "$QC_DIR/fastqc"

echo "$(date '+%F %T')  FastQC on $(ls "$FASTQ_DIR"/*.fastq.gz | wc -l) files, $THREADS threads"
fastqc --threads "$THREADS" --outdir "$QC_DIR/fastqc" "$FASTQ_DIR"/*.fastq.gz

echo "$(date '+%F %T')  MultiQC"
multiqc "$QC_DIR/fastqc" --outdir "$QC_DIR" --filename multiqc_report.html --force

echo "$(date '+%F %T')  Done. Report: $QC_DIR/multiqc_report.html"
