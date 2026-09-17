#!/bin/bash
# =============================================================================
# test.sh — Blind-copy samples and annotate with variants (TNLS)
#
# Usage:
#   ./scripts/modules/statistics/test.sh
#
# Prerequisites:
#   data/modules/statistics/sample_mapping.tsv  — local sample mapping file
#   /home/tan/workspaces/mtdna_raw/temp/test/    — wiped before run
#
# Steps:
#   1. Run blind_copy.py to blind-copy samples, generate metrics JSON, and copy
#      Sequencher TXT variant tables into the sequencher/ subfolder.
#   2. Annotate the sample mapping with variants (add_variants_to_samples.py).
#   3. Write samples.txt — the sample LIDs in encrypted (blind) ID form.
# =============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

WORKSPACE=$(pwd)

DATA_DIR="/mnt/nas/bca/mtDNA/science/data"
RESULTS_DIR="/mnt/nas/bca/mtDNA/science/results"
OUTPUT_DIR="/home/tan/workspaces/mtdna_raw/temp/test"
MAPPING_FILE="${WORKSPACE}/data/modules/statistics/sample_mapping.tsv"

# Clean slate
rm -rf "$OUTPUT_DIR"

# Step 1: Run the blind-copy generator
python -m src.modules.statistics.blind_copy \
    --mapping "$MAPPING_FILE" \
    --output "$OUTPUT_DIR" \
    --data-dir "$DATA_DIR" \
    --results-dir "$RESULTS_DIR"

# Step 2: Annotate samples with variants
python -m src.modules.statistics.add_variants_to_samples \
    -i "${MAPPING_FILE}" \
    -j "${OUTPUT_DIR}/JSON/merged.json" \
    -o "${OUTPUT_DIR}/TNLS_with_variants_$(date +%Y%m%d).tsv"

# Step 3: Write samples.txt — sample LIDs in encrypted (blind) ID form.
# One encrypted ID per line, in mapping order; skips blank lines and # comments,
# mirroring read_mapping_file() in src/modules/statistics/blind_copy.py.
awk -F'\t' 'NF==2 && $1 !~ /^[[:space:]]*#/ {print $1}' "$MAPPING_FILE" > "${OUTPUT_DIR}/samples.txt"
