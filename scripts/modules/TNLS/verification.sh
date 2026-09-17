#!/usr/bin/env bash
# =============================================================================
# verification.sh — End-to-end TNLS/HCLS mtDNA matching pipeline
#
# Usage:
#   bash verification.sh                    Run comparison + statistics (1:1 to 1:8 output)
#   bash verification.sh --skip-filter     Skip LID filter (use full merged set)
#   bash verification.sh --max-ratio 3     Only output 1:1 to 1:3 (default: 15)
#
# Prerequisites:
#   - results/modules/HCLS/<batch>/json/statistic_fullbatch.json  (HCLS samples, produced by Step 0)
#   - data/TNLS/filter.txt                                      (optional LID filter — omit with --skip-filter)
#   - data/TNLS/TNLS_metadata.xlsx                              (optional; enables identity fields)
#
# Outputs:
#   - HCLS samples: results/modules/HCLS/<batch>/
#   - TNLS/HCLS verification: results/modules/TNLS/verification/
#
# =============================================================================

set -euo pipefail

WORKSPACE=$(pwd)

# shellcheck source=/dev/null
source .env

log_message() {
    local message="$1"
    python -c "from loguru import logger; logger.info('${message}')"
}

fail() {
    log_message "ERROR: $1"
    exit 1
}

clean_verification_output() {
    log_message "Cleaning verification output: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
    find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
}

# TNLS data directory (filter list, metadata)
DATA_DIR="$WORKSPACE/data/TNLS"
OUTPUT_DIR="$WORKSPACE/results/modules/TNLS/verification"

HCLS_RAW_ROOT="/mnt/nas/bca/mtDNA/science/HCLS"
HCLS_BATCH="20260317"
HCLS_RAW_DIR="$HCLS_RAW_ROOT/$HCLS_BATCH"
HCLS_OUTPUT_DIR="$WORKSPACE/results/modules/HCLS/$HCLS_BATCH"
HCLS_SAMPLES="/home/tan/workspaces/mtdna_raw/results/tools/sequencher/20260817/json/statistic_fullbatch.json"

# TNLS samples
TNLS_SAMPLES="/mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json"

# LID filter list — optionally filters TNLS_SAMPLES to a subset of LIDs
FILTER_LIST="$DATA_DIR/filter.txt"
TNLS_METADATA="$DATA_DIR/TNLS_metadata.xlsx"

# Matching thresholds — controls Stage 1/2/3 classification logic
MATCHING_RULE="src/modules/TNLS/verification/resources/matching_rule.json"

# ── Options ────────────────────────────────────────────────────────────────────
SKIP_FILTER=1
RUN_STEP0=0
MAX_RATIO=15

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-filter) SKIP_FILTER=1; shift ;;
        --run-step0)   RUN_STEP0=1;  shift ;;
        --max-ratio)   MAX_RATIO="$2"; shift 2 ;;
        *) fail "Unknown option: $1" ;;
    esac
done

# ── Step 0: Generate HCLS samples JSON ────────────────────────────────────────
if (( RUN_STEP0 )); then
    if [[ ! -d "$HCLS_RAW_DIR" ]]; then
        fail "Raw HCLS batch directory not found: $HCLS_RAW_DIR"
    fi
    log_message "Step 0/2: Generate HCLS statistic_fullbatch.json"

    # HCLS data is nested per-LID: <batch>/<LID>/<LID>-<ranges>-<chem>.TXT.
    # Flatten every TXT file into a temp dir, drop the FASTA/MISMATCH exports
    # (not variant tables), then rename to the standard Sequencher format the
    # pipeline QC expects (LID.TXT or LID-a1-b1 a2-b2 ....TXT). Raw names carry a
    # chemistry/run-code suffix (e.g. -22N.045) and hyphen-joined range pairs.
    HCLS_FLAT_DIR=$(mktemp -d -t "hcls_${HCLS_BATCH}_XXXXXX")
    trap 'rm -rf "$HCLS_FLAT_DIR"' EXIT

    # 1) Copy all TXT files flat into the temp dir.
    find "$HCLS_RAW_DIR" -type f -name "*.TXT" -exec cp -t "$HCLS_FLAT_DIR" {} +
    # 2) Filter out FASTA sequence exports and mismatch reports.
    rm -f "$HCLS_FLAT_DIR"/*-FASTA.TXT "$HCLS_FLAT_DIR"/*-MISMATCH.TXT
    # 3) Rename to standard Sequencher filenames so QC accepts them.
    python - "$HCLS_FLAT_DIR" <<'PYSTEP0'
import sys
from pathlib import Path
from src.tools.sequencher.quality_control import parse_filename_ranges

d = Path(sys.argv[1])
for f in sorted(d.glob("*.TXT")):
    p = parse_filename_ranges(f.name)
    if not p.sample_id:
        continue
    ranges = "" if not p.ranges or p.ranges == "FULL REGION" \
        else "-" + " ".join(f"{a}-{b}" for a, b in p.ranges)
    name = f"{p.sample_id}{ranges}.TXT"
    if name != f.name:
        f.rename(d / name)
PYSTEP0

    python -m src.tools.sequencher.pipeline \
        --input-dir "$HCLS_FLAT_DIR" \
        --output-dir "$HCLS_OUTPUT_DIR" \
        --ref-path "$REFERENCE" \
        --batch-id "$HCLS_BATCH"
fi

# ── Clean results ──────────────────────────────────────────────────────────────
# Wipe only TNLS/HCLS verification outputs. HCLS samples live outside this tree
# (results/modules/HCLS/<batch>/), so they survive re-runs when Step 0 is skipped.
clean_verification_output

# ── Validate required inputs ───────────────────────────────────────────────────
log_message "Validating required inputs"
for label in \
    "MATCHING_RULE:$MATCHING_RULE" \
    "HCLS_SAMPLES:$HCLS_SAMPLES" \
    "TNLS_SAMPLES:$TNLS_SAMPLES"
do
    IFS=: read -r name path <<< "$label"
    if [[ ! -f "$path" ]]; then
        fail "$name not found: $path"
    fi
done
log_message "All required files present"

# ── Step 2: Compare HCLS vs TNLS + compute all stats (1:1/1:N ratios, summary TXT, clean output) ────
log_message "Step 2/2: Compare HCLS vs TNLS + statistics"
COMPARE_ARGS=(
    --results-dir "$OUTPUT_DIR"
    --base-samples "$HCLS_SAMPLES"
    --target-samples "$TNLS_SAMPLES"
    --matching-rule "$MATCHING_RULE"
    --max-ratio "$MAX_RATIO"
    --tnls-metadata "$TNLS_METADATA"
)

if (( ! SKIP_FILTER )) && [[ -f "$FILTER_LIST" ]]; then
    COMPARE_ARGS+=(--filter-list "$FILTER_LIST")
fi

python -m src.modules.TNLS.verification.compare_TNLS_HCLS "${COMPARE_ARGS[@]}"
log_message "Done. Results in: $OUTPUT_DIR/"
