#!/bin/bash

# =============================================================================
# Tracy Tool Script — Single-command wrapper around src/tools/tracy/pipeline.py
#
# Runs Tracy batch processing (AB1 decompose → JSON ETL → Sample JSON).
# No multi-step orchestration; simply invokes the Python pipeline.
# =============================================================================

set -eo pipefail

# =============================================================================
# Usage
# =============================================================================

usage() {
    echo "Usage: $0 [options] <BATCH_ID>"
    echo ""
    echo "Run Tracy batch processing for a single batch."
    echo ""
    echo "Options:"
    echo "  -h, --help                Display this help message"
    echo "  -s, --step STEP           (ignored; kept for pipeline.sh compatibility)"
    echo "  --data-dir DIR            Input directory containing AB1 files"
    echo "  --sample-list FILE        Path to TXT file with sample IDs (one per line)"
    echo "  --output-dir DIR          Output directory for results"
    echo "  --ref-path PATH           Path to rCRS reference FASTA (default: \$REFERENCE or ref/rCRS.fasta)"
    echo "  --json-dir DIR            Directory for per-region JSON output"
    echo ""
    echo "Examples:"
    echo "  $0 --data-dir /data/raw/MS_100426_004 \\"
    echo "       --sample-list /data/raw/MS_100426_004.txt \\"
    echo "       --output-dir /results \\"
    echo "       MS_100426_004"
    exit 0
}

# =============================================================================
# Load environment (before arg parsing so CLI args can override)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ENV_FILE="$PROJECT_DIR/.env"

if [ -f "$ENV_FILE" ]; then
    # shellcheck source=/dev/null
    source "$ENV_FILE"
fi

BASE_DIR="${BASE_DIR:-$PROJECT_DIR}"
REFERENCE="${REFERENCE:-$BASE_DIR/ref/rCRS.fasta}"

# Set defaults from env — CLI args will override these
DATA_DIR="${DATA_DIR:-}"
SAMPLE_LIST="${SAMPLE_LIST:-}"
OUTPUT_DIR="${RESULTS_DIR:-}"
REF_PATH="${REF_PATH:-$REFERENCE}"
JSON_DIR="${JSON_DIR:-}"

# =============================================================================
# Argument parsing (CLI args override env defaults)
# =============================================================================

BATCH_ID=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -s|--step)
            # Ignored — kept for pipeline.sh compatibility
            shift 2
            ;;
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --sample-list)
            SAMPLE_LIST="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --ref-path)
            REF_PATH="$2"
            shift 2
            ;;
        --json-dir)
            JSON_DIR="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            ;;
        *)
            BATCH_ID="$1"
            shift
            ;;
    esac
done

# =============================================================================
# Validate required arguments
# =============================================================================

if [ -z "$BATCH_ID" ]; then
    echo "Error: BATCH_ID is required" >&2
    usage
fi
if [ -z "$DATA_DIR" ]; then
    echo "Error: --data-dir is required" >&2
    usage
fi
if [ -z "$SAMPLE_LIST" ]; then
    echo "Error: --sample-list is required" >&2
    usage
fi
if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: --output-dir is required" >&2
    usage
fi

# =============================================================================
# Build and run
# =============================================================================

mkdir -p "$OUTPUT_DIR"

CMD=(
    uv run python -m src.tools.tracy.pipeline
    --input-dir "$DATA_DIR"
    --output-dir "$OUTPUT_DIR"
    --samples "$SAMPLE_LIST"
    --batch-id "$BATCH_ID"
    --ref-path "$REF_PATH"
)

if [ -n "$JSON_DIR" ]; then
    CMD+=(--json-dir "$JSON_DIR")
fi

echo "[tracy] Running Tracy pipeline for ${BATCH_ID}"
echo "[tracy]   data-dir:    ${DATA_DIR}"
echo "[tracy]   sample-list: ${SAMPLE_LIST}"
echo "[tracy]   output-dir:  ${OUTPUT_DIR}"
echo "[tracy]   ref-path:    ${REF_PATH}"
if [ -n "$JSON_DIR" ]; then
    echo "[tracy]   json-dir:    ${JSON_DIR}"
fi

cd "$BASE_DIR"
exec "${CMD[@]}"
