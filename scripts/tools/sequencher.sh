#!/bin/bash

# =============================================================================
# Sequencher Pipeline — Step-Based Shell Script
#
# Runs Sequencher ETL (TXT → Sample JSON) with per-region intermediate JSON
# output. Produces whole-sample JSONs and region JSONs for the unify step.
#
# Supports step selection via -s (numbers or names).
# =============================================================================

# =============================================================================
# Usage
# =============================================================================

usage() {
    echo "Usage: $0 [options] <BATCH_ID> <IN_DIR> <OUT_DIR>"
    echo ""
    echo "Options:"
    echo "  -h, --help              Display this help message"
    echo "  -s, --steps STEPS       Steps to run (comma-separated, numbers or names)"
    echo "  -j, --json-dir DIR      Directory for per-region JSON output (default: \$SEQUENCHER_DIR/\$BATCH_ID)"
    echo ""
    echo "Steps:"
    echo "  1 | extract  - Extract Sequencher ZIP files to temp directory"
    echo "  2 | etl      - Run Sequencher ETL + per-region JSON output"
    echo ""
    echo "Examples:"
    echo "  $0 MS_210426_003 /in /out                                         # Run all"
    echo "  $0 -s 2 MS_210426_003 /in /out                                    # Skip extract"
    echo "  $0 -j /path/to/json MS_210426_003 /in /out                        # Custom JSON dir"
    exit 0
}

# =============================================================================
# Logging
# =============================================================================

log_message() {
    local message="$1"
    python -c "from loguru import logger; logger.info('$message')" | tee -a "$LOG_FILE"
}

# =============================================================================
# Step flags — all default to true
# =============================================================================

RUN_EXTRACT=true
RUN_ETL=true
JSON_DIR=""

BATCH_ID=""
IN_DIR=""
OUT_DIR=""

# =============================================================================
# Argument parsing
# =============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -s|--steps)
            RUN_EXTRACT=false
            RUN_ETL=false

            OLD_IFS="$IFS"
            IFS=',' read -ra STEPS <<< "$2"
            for step in "${STEPS[@]}"; do
                case "$step" in
                    1|extract)  RUN_EXTRACT=true ;;
                    2|etl)      RUN_ETL=true ;;
                    *) echo "Invalid step: $step"; usage ;;
                esac
            done
            IFS="$OLD_IFS"
            shift 2
            ;;
        -j|--json-dir)
            JSON_DIR="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            if [ -z "$BATCH_ID" ]; then
                BATCH_ID="$1"
            elif [ -z "$IN_DIR" ]; then
                IN_DIR="$1"
            elif [ -z "$OUT_DIR" ]; then
                OUT_DIR="$1"
            else
                echo "Error: Too many positional arguments"
                usage
            fi
            shift
            ;;
    esac
done

# =============================================================================
# Validate required arguments
# =============================================================================

if [ -z "$BATCH_ID" ]; then
    echo "Error: BATCH_ID is required"
    usage
fi
if [ -z "$IN_DIR" ]; then
    echo "Error: IN_DIR is required"
    usage
fi
if [ -z "$OUT_DIR" ]; then
    echo "Error: OUT_DIR is required"
    usage
fi

# Source environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ENV_FILE="$PROJECT_DIR/.env"

if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo "Warning: .env not found at $ENV_FILE"
fi

BASE_DIR="${BASE_DIR:-$PROJECT_DIR}"
REFERENCE="${REFERENCE:-$BASE_DIR/ref/rCRS.fasta}"
SEQUENCHER_DIR="${SEQUENCHER_DIR:-/mnt/nas/bca/mtDNA/Sequencher}"

# Default JSON_DIR to SEQUENCHER_DIR/BATCH_ID (region JSONs stay in analysis dir)
if [ -z "$JSON_DIR" ]; then
    JSON_DIR="${SEQUENCHER_DIR}/${BATCH_ID}"
fi

LOG_FILE="${LOGS_DIR}/${BATCH_ID}.log"
mkdir -p "$(dirname "$LOG_FILE")"

log_message "[sequencher] Starting Sequencher pipeline for ${BATCH_ID}"

# =============================================================================
# Step 1: Extract Sequencher ZIP files
# =============================================================================

if [ "$RUN_EXTRACT" == "true" ]; then
    log_message "[Step 1] Extract Sequencher ZIP files"

    mkdir -p "$IN_DIR"

    # Clean previous TXT and JSON output subdirectories only (not the whole batch directory)
    if [ -d "$IN_DIR" ]; then
        log_message "Cleaning TXT files from ${IN_DIR}"
        find "$IN_DIR" -maxdepth 1 -name "*.TXT" -type f -delete
    fi
    if [ -d "$JSON_DIR" ]; then
        log_message "Cleaning region JSON directories in ${JSON_DIR}"
        find "$JSON_DIR" -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} +
    fi

    ZIP_FILES=$(find "$SEQUENCHER_DIR" -type f -name "*${BATCH_ID}*.zip" 2>/dev/null)

    if [ -n "$ZIP_FILES" ]; then
        echo "$ZIP_FILES" | while read -r zip_file; do
            log_message "Extracting $(basename "$zip_file")"
            unzip -q -o "$zip_file" -d "$IN_DIR" >> "$LOG_FILE" 2>&1
        done
    else
        log_message "No ZIP files found for ${BATCH_ID} in ${SEQUENCHER_DIR}"
    fi

    TXT_COUNT=$(find "$IN_DIR" -name "*.TXT" | wc -l)
    log_message "Found ${TXT_COUNT} TXT files in ${IN_DIR}"

    if [ "$TXT_COUNT" -eq 0 ]; then
        echo "Error: No TXT files found in ${IN_DIR} after extraction"
        exit 1
    fi
else
    log_message "Skipping Step 1"
fi

# =============================================================================
# Step 2: Run Sequencher ETL + per-region JSON output
# =============================================================================

if [ "$RUN_ETL" == "true" ]; then
    log_message "[Step 2] Sequencher ETL"

    if [ ! -d "$IN_DIR" ]; then
        echo "Error: Input directory not found: ${IN_DIR} (run step 1 first)"
        exit 1
    fi

    python -m src.tools.sequencher.pipeline \
        --input-dir "$IN_DIR" \
        --output-dir "$OUT_DIR" \
        --ref-path "$REFERENCE" \
        --json-dir "$JSON_DIR" \
        --batch-id "$BATCH_ID" >> "$LOG_FILE" 2>&1

    if [ $? -ne 0 ]; then
        echo "Error: Sequencher ETL failed (see $LOG_FILE)"
        exit 1
    fi

    log_message "Sequencher ETL complete — region JSONs in ${JSON_DIR}/"
else
    log_message "Skipping Step 2"
fi

# =============================================================================
# Done
# =============================================================================

log_message "Sequencher pipeline completed for ${BATCH_ID}"

# =============================================================================
# Cleanup: Remove extracted TXT directory and batch directory after successful run
# =============================================================================
# The txt/ directory is a temp working area holding files extracted from ZIPs.
# On re-run, Step 1 re-extracts from the source ZIPs, so this is safe to delete.
# The batch directory under SEQUENCHER_DIR is also temp (only holds txt/ here),
# so remove it too to avoid leaving empty stubs like .../Sequencher_temp/MS_210426_003.

if [ "$RUN_ETL" == "true" ] && [ -d "$IN_DIR" ]; then
    log_message "Cleaning up extracted TXT directory: ${IN_DIR}"
    rm -rf "${IN_DIR}"

    BATCH_DIR="${SEQUENCHER_DIR}/${BATCH_ID}"
    if [ -d "$BATCH_DIR" ]; then
        # Remove the batch directory if empty (only txt/ was in it)
        rmdir "$BATCH_DIR" 2>/dev/null && \
            log_message "Removed empty batch directory: ${BATCH_DIR}" || \
            log_message "Batch directory not empty, keeping: ${BATCH_DIR}"
    fi
fi
