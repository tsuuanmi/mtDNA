#!/bin/bash
# =============================================================================
# run.sh — Run the mtDNA LAB HTML report generator for all FASTA batches
#
# Usage:
#   ./scripts/modules/LAB/run.sh
#
# Environment:
#   Requires .env with RESULTS_DIR and DATA_DIR defined.
#
# Inputs:
#   data/modules/LAB/tracking.xlsx
#   data/modules/LAB/rerun.xlsx
#
# Outputs go to: temp/modules/LAB/<batch>.html and temp/modules/LAB/<batch>.xlsx
# =============================================================================

set -e
set -u

# Load environment variables
source .env

WORKSPACE=$(pwd)
TRACKING_EXCEL="${WORKSPACE}/data/modules/LAB/tracking.xlsx"
TRACKING_TSV="${WORKSPACE}/data/modules/LAB/tracking.tsv"
RERUN_EXCEL="${WORKSPACE}/data/modules/LAB/rerun.xlsx"
FASTA_DIR="${RESULTS_DIR}/fasta"
ARCHIVE_DIR="${RESULTS_DIR}/archive"
METADATA_DIR="${DATA_DIR}/metadata"

# Regenerate tracking TSV from the current workbook before generating reports.
python -m src.modules.LAB.manual \
    --tracking-excel "$TRACKING_EXCEL" \
    --output "$TRACKING_TSV"

# Collect each batch containing FASTA files from active and archived results.
shopt -s nullglob
declare -A SEEN_BATCHES=()
BATCHES=()
for ROOT in "$FASTA_DIR" "$ARCHIVE_DIR"; do
    if [ ! -d "$ROOT" ]; then
        echo ">>> FASTA directory not found, skipping: $ROOT"
        continue
    fi

    for FOLDER in "$ROOT"/*; do
        [ -d "$FOLDER" ] || continue
        compgen -G "$FOLDER/*.fasta" > /dev/null || continue

        BATCH=$(basename "$FOLDER")
        if [ -n "${SEEN_BATCHES[$BATCH]+x}" ]; then
            continue
        fi

        SEEN_BATCHES["$BATCH"]=1
        BATCHES+=("$BATCH")
    done
done

if [ "${#BATCHES[@]}" -eq 0 ]; then
    echo ">>> No batches with FASTA files found in $FASTA_DIR or $ARCHIVE_DIR"
    exit 0
fi

# Process each batch unless both report artifacts already exist.
for BATCH in "${BATCHES[@]}"; do
    OUTPUT_PREFIX="${WORKSPACE}/temp/modules/LAB/${BATCH}"
    HTML_OUTPUT="${OUTPUT_PREFIX}.html"
    EXCEL_OUTPUT="${OUTPUT_PREFIX}.xlsx"

    if [ -f "$HTML_OUTPUT" ] && [ -f "$EXCEL_OUTPUT" ]; then
        echo ">>> Skipping completed batch: $BATCH"
        continue
    fi

    METADATA_EXCEL="${METADATA_DIR}/${BATCH}.xlsx"
    if [ ! -f "$METADATA_EXCEL" ]; then
        echo ">>> Metadata Excel not found; cannot generate report: $METADATA_EXCEL" >&2
        exit 1
    fi

    echo ">>> Processing batch: $BATCH"
    python -m src.modules.LAB.generate_mtdna_html_report \
        -b "$BATCH" \
        -m "$METADATA_EXCEL" \
        -o "$OUTPUT_PREFIX" \
        -r "$RERUN_EXCEL" \
        -t "$TRACKING_TSV"
done
