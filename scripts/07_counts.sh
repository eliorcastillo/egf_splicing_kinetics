#!/usr/bin/env bash
# Step 7 - Gene-level counts with featureCounts.
# Single-end, unstranded library (-s 0, verified in step 6).
# Usage (egf-ngs environment active):
#   bash scripts/07_counts.sh 6

set -euo pipefail

THREADS="${1:-6}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
GTF="$PROJECT_DIR/data/ref/gencode.v50.primary_assembly.annotation.gtf"
BAM_DIR="$PROJECT_DIR/results/bam"
OUT_DIR="$PROJECT_DIR/results/counts"
mkdir -p "$OUT_DIR"

featureCounts -T "$THREADS" -s 0 \
    -a "$GTF" \
    --extraAttributes gene_name,gene_type \
    -o "$OUT_DIR/gene_counts.txt" \
    "$BAM_DIR"/*.bam
