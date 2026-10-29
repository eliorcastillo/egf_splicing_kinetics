#!/usr/bin/env bash
# Step 6 - Align all FASTQ files with HISAT2 -> sorted, indexed BAM.
# Single-end, unstranded library (see README).
# Usage (egf-ngs environment active):
#   caffeinate -i bash scripts/06_align.sh 6

set -euo pipefail

THREADS="${1:-6}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
FASTQ_DIR="$PROJECT_DIR/data/raw/fastq"
INDEX="$PROJECT_DIR/data/ref/hisat2_index/GRCh38"
SPLICE="$PROJECT_DIR/data/ref/splice_sites.txt"
BAM_DIR="$PROJECT_DIR/results/bam"
LOG_DIR="$PROJECT_DIR/results/logs"
mkdir -p "$BAM_DIR" "$LOG_DIR"

LOG="$LOG_DIR/06_align.log"
log() { echo "$(date '+%F %T')  $*" | tee -a "$LOG"; }

for FQ in "$FASTQ_DIR"/*.fastq.gz; do
    SAMPLE="$(basename "$FQ" .fastq.gz)"
    BAM="$BAM_DIR/$SAMPLE.bam"

    if [[ -s "$BAM.bai" ]]; then
        log "$SAMPLE: already aligned, skipping"
        continue
    fi

    log "$SAMPLE: aligning"
    hisat2 -p "$THREADS" \
        -x "$INDEX" \
        -U "$FQ" \
        --known-splicesite-infile "$SPLICE" \
        --new-summary --summary-file "$LOG_DIR/$SAMPLE.hisat2.txt" \
        2>> "$LOG" \
    | samtools sort -@ 4 -m 1G -o "$BAM" -
    samtools index "$BAM"

    RATE="$(grep "Overall alignment rate" "$LOG_DIR/$SAMPLE.hisat2.txt")"
    log "$SAMPLE: done ($RATE)"
done

log "All samples aligned."
