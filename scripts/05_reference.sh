#!/usr/bin/env bash
# Step 5 - Reference genome, annotation and HISAT2 index (run once).
# Source: GENCODE human release 50 (GRCh38.p14, primary assembly).
# Usage (egf-ngs environment active):
#   caffeinate -i bash scripts/05_reference.sh 6

set -euo pipefail

THREADS="${1:-6}"
RELEASE="50"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}"
GENOME="GRCh38.primary_assembly.genome.fa.gz"
GTF="gencode.v${RELEASE}.primary_assembly.annotation.gtf.gz"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
REF_DIR="$PROJECT_DIR/data/ref"
INDEX_DIR="$REF_DIR/hisat2_index"
mkdir -p "$INDEX_DIR"
cd "$REF_DIR"

log() { echo "$(date '+%F %T')  $*"; }

# 1. Download genome, annotation and the official checksums
log "Downloading GENCODE release $RELEASE"
for FILE in "$GENOME" "$GTF" MD5SUMS; do
    curl --fail --location --continue-at - --remote-name "$BASE_URL/$FILE"
done

# 2. Verify the downloads against the checksums published by GENCODE
log "Checking MD5 checksums"
for FILE in "$GENOME" "$GTF"; do
    EXPECTED=$(grep " $FILE\$" MD5SUMS | cut -d ' ' -f 1)
    ACTUAL=$(md5 -q "$FILE")
    if [ "$EXPECTED" != "$ACTUAL" ]; then
        log "ERROR: checksum mismatch for $FILE"
        exit 1
    fi
    log "  $FILE OK"
done

# 3. Decompress (keep the .gz copies with -k)
log "Decompressing"
gunzip --keep --force "$GENOME" "$GTF"

# 4. Splice sites from the annotation (used later by HISAT2 during alignment)
log "Extracting known splice sites"
hisat2_extract_splice_sites.py "${GTF%.gz}" > splice_sites.txt

# 5. Build the HISAT2 index (genome only; about 8 GB RAM, 1-2 h)
log "Building HISAT2 index with $THREADS threads"
hisat2-build -p "$THREADS" "${GENOME%.gz}" "$INDEX_DIR/GRCh38"

log "Done. Index files:"
ls -lh "$INDEX_DIR"
