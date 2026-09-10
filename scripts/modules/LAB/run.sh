#!/bin/bash
# =============================================================================
# run.sh — Run the mtDNA LAB HTML report generator for all batches
#
# Usage:
#   ./scripts/modules/lab/run.sh           # Process all batches from fasta/
#   ./scripts/modules/lab/run.sh <batch>    # Process a specific batch only
#
# Environment:
#   Requires .env with RESULTS_DIR, DATA_DIR defined.
#
# Outputs go to: temp/<batch>/
# =============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load environment variables
source .env

WORKSPACE=$(pwd)

# Dynamically list all batches from fasta folder
BATCHES=()
for folder in "${RESULTS_DIR}/fasta"/*; do
    if [ -d "$folder" ]; then
        batch_name=$(basename "$folder")
        BATCHES+=("$batch_name")
    fi
done

# Process each batch
for BATCH in "${BATCHES[@]}"; do
    echo ">>> Processing batch: $BATCH"
    python -m src.modules.LAB.generate_mtdna_html_report \
        -b "$BATCH" \
        -m /mnt/bca/mtDNA/science/data/metadata/$BATCH.xlsx \
        -o temp/$BATCH \
        -r "${WORKSPACE}/src/modules/LAB/rerun.xlsx" \
        -t "${WORKSPACE}/src/modules/LAB/tracking.tsv"
done
