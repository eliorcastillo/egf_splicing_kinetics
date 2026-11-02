#!/usr/bin/env bash
# Step 9 - Splicing events and PSI per sample with rMATS.
# All 18 samples in one group, no statistics (--statoff): PSI per event per sample.
# Single-end, unstranded, 101 bp (checked on SRR2061414: 100% of reads = 101 bp).
# HISAT2 soft-clips ~12% of reads -> --allow-clipping so they are not discarded.
# Run inside the egf-rmats conda environment.
set -euo pipefail

THREADS="${1:-6}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
GTF="$PROJECT_DIR/data/ref/gencode.v50.primary_assembly.annotation.gtf"
BAM_DIR="$PROJECT_DIR/results/bam"
SAMPLES="$PROJECT_DIR/docs/metadata/samples.tsv"
RMATS_DIR="$PROJECT_DIR/results/rmats"
OUT_DIR="$RMATS_DIR/all_samples"
TMP_DIR="$RMATS_DIR/tmp_all"
LOG="$PROJECT_DIR/results/logs/09_rmats.log"
mkdir -p "$OUT_DIR" "$TMP_DIR" "$(dirname "$LOG")"

# BAM list in samples.tsv order (t000_r1 ... t300_r3), comma-separated, one line
B1="$RMATS_DIR/b1_all_samples.txt"
tail -n +2 "$SAMPLES" | cut -f1 | sed "s#.*#$BAM_DIR/&.bam#" | paste -sd, - > "$B1"

# Same order with sample names: column k of IncLevel1 = row k of this file
tail -n +2 "$SAMPLES" | cut -f1,3 > "$RMATS_DIR/sample_order.tsv"

# Stop if any BAM is missing
for BAM in $(tr ',' ' ' < "$B1"); do
    [[ -s "$BAM" ]] || { echo "Missing: $BAM"; exit 1; }
done

echo "$(date '+%F %T')  rMATS start ($(rmats.py --version))" | tee -a "$LOG"
rmats.py \
    --b1 "$B1" \
    --gtf "$GTF" \
    -t single \
    --readLength 101 \
    --allow-clipping \
    --libType fr-unstranded \
    --nthread "$THREADS" \
    --statoff \
    --od "$OUT_DIR" \
    --tmp "$TMP_DIR" \
    >> "$LOG" 2>&1
echo "$(date '+%F %T')  rMATS done" | tee -a "$LOG"
