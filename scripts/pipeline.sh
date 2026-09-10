#!/bin/bash

# Exit on error
set -eo pipefail

# =============================================================================
# Unified mtDNA Pipeline
# Supports two analysis engines: 'tracy' (default) and 'blastn'
# =============================================================================

# Define usage and help function
usage() {
    echo "Usage: $0 [options] <batch_id>"
    echo ""
    echo "Options:"
    echo "  -h, --help              Display this help message"
    echo "  -p, --pipeline MODE     Analysis engine: 'tracy' (default) or 'blastn'"
    echo "  -s, --steps STEPS       Steps to run (comma-separated, numbers or names)"
    echo "  -d, --data-dir DIR      Specify custom DATA_DIR path (isolated run, requires -l)"
    echo "  -m, --manual-dir DIR    Specify custom MANUAL_DIR path"
    echo "  -o, --output-dir DIR    Specify custom RESULTS_DIR path"
    echo "  -l, --sample-list FILE  Custom sample list file (skips Excel extraction)"
    echo "      --sequencher-dir DIR  Override SEQUENCHER_DIR from .env (Sequencher zip source)"
    echo "  --force, -f             Force regenerate all FASTA files (wipe output dir first)"
    echo ""
    echo "Steps (use numbers or names with -s):"
    echo "  1 | validate   - Positive Control validation"
    echo "  2 | prepare    - Organize AB1 files and extract sample list"
    echo "  3 | automate   - Run analysis engine (varies by --pipeline mode)"
    echo "  4 | manual     - Extract manual data + Sequencher comparison"
    echo "  5 | fasta      - Generate FASTA files"
    echo "  6 | report     - Generate PDF reports"
    echo ""
    echo "Pipeline modes (for -p):"
    echo "  tracy   - Integrated Tracy analysis with variant calling thresholds"
    echo "  blastn  - BLASTN analysis → per-sample Tracy decompose → automate pipeline"
    echo ""
    echo "Examples:"
    echo "  $0 20241218_mtDNA_11                              # Run all steps with tracy (default)"
    echo "  $0 -p blastn 20241218_mtDNA_11                    # Run all steps with blastn engine"
    echo "  $0 -s 2,3,4 20241218_mtDNA_11                     # Run prepare, automate, manual only"
    echo "  $0 -s prepare,automate,manual 20241218_mtDNA_11   # Same as above, using names"
    echo "  $0 -l /path/to/samples.txt -s 3,4 20241218_mtDNA_11  # Rerun specific samples"
    echo "  $0 -d /custom/data -l /path/to/samples.txt 20241218_mtDNA_11  # Isolated run (requires -l)"
    exit 0
}

# =============================================================================
# Default configuration
# =============================================================================

PIPELINE_MODE="tracy"

# Step flags — all default to true
RUN_VALIDATE=true
RUN_PREPARE=true
RUN_AUTOMATE=true
RUN_MANUAL=true
RUN_FASTA=true
RUN_REPORT=true
FORCE_MERGE=""

BATCH_ID=""
CUSTOM_DATA_DIR=""
CUSTOM_MANUAL_DIR=""
CUSTOM_RESULTS_DIR=""
CUSTOM_SAMPLE_LIST=""
CUSTOM_SEQUENCHER_DIR=""

# =============================================================================
# Argument parsing
# =============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -p|--pipeline)
            PIPELINE_MODE="$2"
            if [[ "${PIPELINE_MODE}" != "blastn" && "${PIPELINE_MODE}" != "tracy" ]]; then
                echo "Error: Invalid pipeline mode '${PIPELINE_MODE}'. Use 'blastn' or 'tracy'."
                exit 1
            fi
            shift 2
            ;;
        -d|--data-dir)
            CUSTOM_DATA_DIR="$2"
            shift 2
            ;;
        -m|--manual-dir)
            CUSTOM_MANUAL_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            CUSTOM_RESULTS_DIR="$2"
            shift 2
            ;;
        -l|--sample-list)
            CUSTOM_SAMPLE_LIST="$2"
            shift 2
            ;;
        --sequencher-dir)
            CUSTOM_SEQUENCHER_DIR="$2"
            shift 2
            ;;
        --force|-f)
            FORCE_MERGE=true
            shift
            ;;
        -s|--steps)
            # Reset all steps to false first
            RUN_VALIDATE=false
            RUN_PREPARE=false
            RUN_AUTOMATE=false
            RUN_MANUAL=false
            RUN_FASTA=false
            RUN_REPORT=false

            # Parse comma-separated list of steps (numbers or names)
            OLD_IFS="${IFS}"
            IFS=',' read -ra STEPS <<< "$2"
            for step in "${STEPS[@]}"; do
                case "${step}" in
                    1|validate) RUN_VALIDATE=true ;;
                    2|prepare)  RUN_PREPARE=true ;;
                    3|automate) RUN_AUTOMATE=true ;;
                    4|manual)   RUN_MANUAL=true ;;
                    5|fasta)    RUN_FASTA=true ;;
                    6|report)   RUN_REPORT=true ;;
                    *) echo "Invalid step: ${step}"; usage ;;
                esac
            done
            IFS="${OLD_IFS}"
            shift 2
            ;;
        *)
            if [[ "$1" != -* ]]; then
                BATCH_ID="$1"
                shift
            else
                echo "Unknown option: $1"
                usage
            fi
            ;;
    esac
done

# =============================================================================
# Load configuration
# =============================================================================

# shellcheck source=/dev/null
source .env

# Override DATA_DIR if custom path is provided
if [ -n "${CUSTOM_DATA_DIR}" ]; then
    DATA_DIR="${CUSTOM_DATA_DIR}"
fi

# Override RESULTS_DIR if custom path is provided
if [ -n "${CUSTOM_RESULTS_DIR}" ]; then
    RESULTS_DIR="${CUSTOM_RESULTS_DIR}"
fi

# Override SEQUENCHER_DIR if custom path is provided.
# Must come after `source .env`, which unconditionally reassigns it.
if [ -n "${CUSTOM_SEQUENCHER_DIR}" ]; then
    SEQUENCHER_DIR="${CUSTOM_SEQUENCHER_DIR}"
fi

# Require batch ID as argument
if [ -z "${BATCH_ID}" ]; then
    echo "Error: No batch ID provided as argument"
    echo "Usage: $0 [options] <batch_id>"
    echo "Run '$0 --help' for more information."
    exit 1
fi

# Validate custom directories exist
if [ -n "${CUSTOM_DATA_DIR}" ] && [ ! -d "${CUSTOM_DATA_DIR}" ]; then
    echo "Error: Custom DATA_DIR '${CUSTOM_DATA_DIR}' does not exist"
    exit 1
fi

if [ -n "${CUSTOM_MANUAL_DIR}" ] && [ ! -d "${CUSTOM_MANUAL_DIR}" ]; then
    echo "Error: Custom MANUAL_DIR '${CUSTOM_MANUAL_DIR}' does not exist"
    exit 1
fi

if [ -n "${CUSTOM_RESULTS_DIR}" ] && [ ! -d "${CUSTOM_RESULTS_DIR}" ]; then
    echo "Error: Custom RESULTS_DIR '${CUSTOM_RESULTS_DIR}' does not exist"
    exit 1
fi

if [ -n "${CUSTOM_SAMPLE_LIST}" ] && [ ! -f "${CUSTOM_SAMPLE_LIST}" ]; then
    echo "Error: Custom sample list file '${CUSTOM_SAMPLE_LIST}' does not exist"
    exit 1
fi

# Validate -d requires -l for isolated run
if [ -n "${CUSTOM_DATA_DIR}" ] && [ -z "${CUSTOM_SAMPLE_LIST}" ]; then
    echo "Error: -d requires -l for isolated run. Use: -d /custom/data -l /path/to/samples.txt"
    exit 1
fi

# Extract date from BATCH_ID (format: YYYYMMDD_mtDNA_XX)
YEAR=${BATCH_ID:0:4}
MONTH=${BATCH_ID:4:2}
DAY=${BATCH_ID:6:2}
DATE="${DAY}-${MONTH}-${YEAR}"

# =============================================================================
# Logging
# =============================================================================

LOG_FILE="${LOGS_DIR}/${BATCH_ID}.log"

log_message() {
    local message="$1"
    python -c "from loguru import logger; logger.info('${message}')" | tee -a "${LOG_FILE}"
}

log_message "=========================================="
log_message "Pipeline mode: ${PIPELINE_MODE}"
log_message "Batch ID: ${BATCH_ID}"
log_message "=========================================="

# =============================================================================
# Directory structure
# =============================================================================

RAW_DIR="${DATA_DIR}/raw/${BATCH_ID}"
METADATA_DIR="${DATA_DIR}/metadata"

# Sample list the analysis step actually reads. A custom list (-l) is used
# where it lies instead of being copied over "${RAW_DIR}.txt": that file is
# the batch's canonical list on shared storage, while -l carries only a
# subset (the rerun samples). Resolving it here rather than inside step 2 is
# what makes `-s 3,4 -l list.txt` work -- the SOP rerun skips step 2, so a
# list consumed only there never reached the analysis at all.
if [ -n "${CUSTOM_SAMPLE_LIST}" ]; then
    SAMPLE_LIST_FILE="${CUSTOM_SAMPLE_LIST}"
else
    SAMPLE_LIST_FILE="${RAW_DIR}.txt"
fi

AUTOMATE_DIR="${RESULTS_DIR}/automate_pipeline"
TRACY_RESULTS_DIR="${RESULTS_DIR}/tools/tracy/${BATCH_ID}"
BLASTN_RESULTS_DIR="${RESULTS_DIR}/tools/blastn/${BATCH_ID}"
AUTOMATE_REGENERATE_DIR="${RESULTS_DIR}/automate_pipeline/regenerate/${BATCH_ID}"

MANUAL_DIR="${RESULTS_DIR}/manual_pipeline/manual/${BATCH_ID}"
SEQUENCHER_RESULTS_DIR="${RESULTS_DIR}/tools/sequencher"
MANUAL_REGENERATE_DIR="${RESULTS_DIR}/manual_pipeline/regenerate/${BATCH_ID}"

# Expected PC control profiles (Sequencher TXT); defaults to project ./data
PC_REFERENCE_DIR="${PC_REFERENCE_DIR:-data}"

COMPARISON_DIR="${RESULTS_DIR}/modules/comparison/${BATCH_ID}"

FASTA_DIR="${RESULTS_DIR}/fasta/${BATCH_ID}"
REPORTS_DIR="${RESULTS_DIR}/reports/${BATCH_ID}"

if [ -n "${CUSTOM_MANUAL_DIR}" ]; then
    MANUAL_DIR="${CUSTOM_MANUAL_DIR}"
fi

log_message "RAW_DIR: ${RAW_DIR}"
log_message "MANUAL_DIR: ${MANUAL_DIR}"
log_message "RESULTS_DIR: ${RESULTS_DIR}"

# Create required directories based on which steps are running
mkdir_if_needed() {
    if [ "$1" == "true" ]; then
        for dir in "${@:2}"; do
            mkdir -p "${dir}"
        done
    fi
}

# Filter a tool's raw JSON output into the canonical regenerate (final) directory.
# Variants are scoped to the sequenced intervals; the raw tool JSON is retained
# unchanged under results/tools/<tool>/json for the comparison step, so the
# reviewer still sees out-of-range/additional variants there.
run_regenerate() {
    local src="${1}"
    local dest="${2}"
    if [ ! -d "${src}" ]; then
        log_message "Warning: source ${src} not found; ${dest} left unchanged"
        return 0
    fi
    rm -rf "${dest:?}"
    mkdir -p "${dest}"
    python -m src.core.regenerate \
        -i "${src}" \
        -o "${dest}" \
        -r "${REFERENCE}" \
        --batch-id "${BATCH_ID}" >> "${LOG_FILE}" 2>&1
    log_message "Regenerated (filtered)"
}

mkdir_if_needed "${RUN_VALIDATE}" "${PC_DATA_DIR}/${BATCH_ID}" "${PC_RESULTS_DIR}/${BATCH_ID}"
mkdir_if_needed "${RUN_PREPARE}" "${RAW_DIR}" "${METADATA_DIR}"
mkdir_if_needed "${RUN_AUTOMATE}" "${BLASTN_RESULTS_DIR}" "${TRACY_RESULTS_DIR}" "${AUTOMATE_REGENERATE_DIR}"
mkdir_if_needed "${RUN_MANUAL}" "${MANUAL_DIR}" "${MANUAL_REGENERATE_DIR}"
mkdir_if_needed "${RUN_FASTA}" "${FASTA_DIR}"
mkdir_if_needed "${RUN_REPORT}" "${REPORTS_DIR}"

# Chỉ để hiển thị trong log. Đừng dùng để khớp file: nó cắt mất ngày,
# và khớp chuỗi con thì mtDNA_43 trúng luôn mtDNA_430.
BATCH="${BATCH_ID#*_}"

# Batch-scoped match patterns. The old `*${BATCH}*` form stripped the leading
# field and matched as a substring, so e.g. batch 20241222_mtDNA_12 also pulled
# in mtDNA_120, mtDNA_121 and a same-numbered batch from another date --
# silently mixing other batches' AB1/TXT into this run. Anchor on the full
# BATCH_ID instead and require a separator (_ . or space) after it.
BATCH_DIR_GLOB="${BATCH_ID}[_. ]*"
BATCH_ZIP_GLOB="${BATCH_ID}[_. ]*.zip"

# =============================================================================
# Step 1: Validate — Positive Control processing
# =============================================================================

if [ "${RUN_VALIDATE}" == "true" ]; then
    log_message "[Step 1/6 - validate] Running Positive Control validation"

    # pc_ntc sources its AB1 files from --data-dir. For an isolated/local run
    # (-d/--data-dir) that is the user-provided batch dir; for a normal NAS run
    # it is the lab batch directory LAB_DATA_DIR/<batch_id>.
    if [ -n "${CUSTOM_DATA_DIR}" ]; then
        PC_NT_DATA_DIR="${DATA_DIR}"
    else
        PC_NT_DATA_DIR="${LAB_DATA_DIR}/${BATCH_ID}"
    fi

    if python -m src.modules.quality_control.pc_ntc \
        --batch-id "${BATCH_ID}" \
        --data-dir "${PC_NT_DATA_DIR}" \
        --results-dir "${RESULTS_DIR}" \
        --ref-path "${REFERENCE}" \
        --pc-reference-dir "${PC_REFERENCE_DIR}" >> "${LOG_FILE}" 2>&1; then
        log_message "Positive Control validation completed successfully"
    else
        log_message "Positive Control validation failed, but continuing with main pipeline"
    fi
else
    log_message "Skipping Step 1: Positive Control validation"
fi

# =============================================================================
# Step 2: Prepare — Organize AB1 files and extract sample list
# =============================================================================

if [ "${RUN_PREPARE}" == "true" ]; then
    if [ -n "${CUSTOM_DATA_DIR}" ]; then

        # Organize AB1 files from the custom DATA_DIR into RAW_DIR (mirror normal run)
        log_message "[Step 2/6 - prepare] Organizing AB1 files from custom DATA_DIR: ${DATA_DIR}"
        # Thư mục này có thể đang giữ read chạy lại do rerun đưa vào. Bước 2
        # nạp lại từ LAB_DATA_DIR -- nơi chỉ có read gốc -- nên sẽ mang đúng
        # read đã hỏng quay về mà không ai biết.
        if [ -f "${RAW_DIR}/.rerun_applied.json" ]; then
            log_message "WARNING: ${RAW_DIR} dang giu read chay lai (.rerun_applied.json)"
            log_message "WARNING: Buoc 2 se nap lai tu nguon goc va dua read cu quay ve"
        fi
        [ -d "${RAW_DIR}" ] && rm -rf "${RAW_DIR}"

        mkdir -p "${RAW_DIR}"

        if [ -f "${CUSTOM_SAMPLE_LIST}" ]; then
            log_message "Using custom sample list: ${CUSTOM_SAMPLE_LIST}"
            cp "${CUSTOM_SAMPLE_LIST}" "${RAW_DIR}.txt"
        fi

        find "${DATA_DIR}" -type f -name "*.ab1" ! -path "${RAW_DIR}/*" -exec cp {} "${RAW_DIR}" \;
        find "${RAW_DIR}" -type f -name "*.ab1" \( -name "*NTC*" -o -name "*PCR*" -o -name "*PC*" \) -delete
        find "${RAW_DIR}" -mindepth 1 -type d -exec rm -rf {} + 2>/dev/null || true
    else
        log_message "[Step 2/6 - prepare] Organizing AB1 files"

        # Thư mục này có thể đang giữ read chạy lại do rerun đưa vào. Bước 2
        # nạp lại từ LAB_DATA_DIR -- nơi chỉ có read gốc -- nên sẽ mang đúng
        # read đã hỏng quay về mà không ai biết.
        if [ -f "${RAW_DIR}/.rerun_applied.json" ]; then
            log_message "WARNING: ${RAW_DIR} dang giu read chay lai (.rerun_applied.json)"
            log_message "WARNING: Buoc 2 se nap lai tu nguon goc va dua read cu quay ve"
        fi

        # Remove existing RAW_DIR if it exists
        if [ -d "${RAW_DIR}" ]; then
            rm -rf "${RAW_DIR}"
        fi

        mkdir -p "${RAW_DIR}"

        # Copy only .ab1 files from the source directory.
        # Neo vào trọn BATCH_ID như phần copy metadata ngay bên dưới. Dạng cũ
        # `-path "*${BATCH}*"` cắt mất ngày rồi khớp chuỗi con, nên batch
        # 20250624_mtDNA_43 kéo về cả mtDNA_430/431/…  (3850 file thay vì 2),
        # batch 20250702_mtDNA_136 nhận trọn 377 read của 20250705_mtDNA_136 --
        # hai batch đó dùng chung cả 90 mã mẫu nên tracy nhận nhầm mà không
        # cách nào biết -- và cả thư mục rerun/ cũng bị lôi vào.
        log_message "Copying .ab1 files from source directory"
        find "${LAB_DATA_DIR}" -type f -name "*.ab1" \
            \( -path "*/${BATCH_ID}/*" -o -path "*/${BATCH_DIR_GLOB}/*" \) \
            -exec cp {} "${RAW_DIR}" \; > /dev/null 2>&1

        # Remove control samples and organize files
        find "${RAW_DIR}" -type f -name "*.ab1" \( -name "*NTC*" -o -name "*PCR*" -o -name "*PC*" \) -delete
        find "${RAW_DIR}" -type f -name "*.ab1" ! -path "${RAW_DIR}/*" -exec mv {} "${RAW_DIR}" \; 2>/dev/null || true
        find "${RAW_DIR}" -mindepth 1 -type d -exec rm -rf {} + 2>/dev/null || true
        find "${RAW_DIR}" -type f ! -name "*.ab1" -delete

        # Copy Excel metadata file
        log_message "Copying Excel metadata file"
        find "${LAB_DATA_DIR}" -type f -name "*.xlsx" \
            \( -path "*/${BATCH_ID}/*" -o -path "*/${BATCH_DIR_GLOB}/*" \) \
            -exec cp -f "{}" "${METADATA_DIR}/${BATCH_ID}.xlsx" \;

        # Use custom sample list if provided, otherwise extract from Excel
        if [ -n "${CUSTOM_SAMPLE_LIST}" ] && [ -f "${CUSTOM_SAMPLE_LIST}" ]; then
            # Read from where it lies; see SAMPLE_LIST_FILE above for why it is
            # not copied over the canonical "${RAW_DIR}.txt".
            log_message "Using custom sample list: ${CUSTOM_SAMPLE_LIST}"
        else
            # Generate sample list from Excel file LID column (excluding control samples)
            log_message "Extracting sample IDs from Excel file (excluding PC and NTC controls)"
            python -c "import pandas as pd; df = pd.read_excel('${METADATA_DIR}/${BATCH_ID}.xlsx'); lid_series = df['LID'].dropna().drop_duplicates().astype(str); filtered = lid_series[~lid_series.str.contains('PC|NTC', case=False, na=False)]; filtered.to_csv('${RAW_DIR}.txt', index=False, header=False)"
        fi
        # Ensure file ends with newline
        [ -s "${RAW_DIR}.txt" ] && [ "$(tail -c1 "${RAW_DIR}.txt" | wc -l)" -eq 0 ] && echo >> "${RAW_DIR}.txt"
    fi
else
    log_message "Skipping Step 2: Prepare input data"
fi

# =============================================================================
# Step 3: Automate — Routes by --pipeline mode
# =============================================================================

if [ "${RUN_AUTOMATE}" == "true" ]; then

    if [ "${PIPELINE_MODE}" == "blastn" ]; then
        # -----------------------------------------------------------------
        # BLASTN pipeline (direct src.tools.blastn module, same pattern as Tracy)
        # -----------------------------------------------------------------
        log_message "[Step 3/6 - automate] Running BLASTN pipeline"

        # Wipe the tool batch results folder to ensure a fresh run
        if [ -d "${BLASTN_RESULTS_DIR}" ]; then
            log_message "Cleaning up existing BLASTN output"
            rm -rf "${BLASTN_RESULTS_DIR:?}"/*
        fi
        mkdir -p "${BLASTN_RESULTS_DIR}"

        python -m src.tools.blastn.pipeline \
            --samples "${SAMPLE_LIST_FILE}" \
            --batch-id "${BATCH_ID}" \
            --input-dir "${DATA_DIR}/raw/${BATCH_ID}" \
            --output-dir "${BLASTN_RESULTS_DIR}" \
            --ref-path "${REFERENCE}" >> "${LOG_FILE}" 2>&1

        AUTOMATE_RAW_JSON_DIR="${BLASTN_RESULTS_DIR}/json"
        run_regenerate "${AUTOMATE_RAW_JSON_DIR}" "${AUTOMATE_REGENERATE_DIR}"

    elif [ "${PIPELINE_MODE}" == "tracy" ]; then
        # -----------------------------------------------------------------
        # Tracy-only pipeline (standard src.tools.tracy module)
        # -----------------------------------------------------------------
        log_message "[Step 3/6 - automate] Running Tracy-only pipeline"

        # Wipe the tool batch results folder to ensure a fresh run
        if [ -d "${TRACY_RESULTS_DIR}" ]; then
            log_message "Cleaning up existing Tracy output files"
            rm -rf "${TRACY_RESULTS_DIR:?}"/*
        fi
        mkdir -p "${TRACY_RESULTS_DIR}"

        python -m src.tools.tracy.pipeline \
            --samples "${SAMPLE_LIST_FILE}" \
            --batch-id "${BATCH_ID}" \
            --input-dir "${DATA_DIR}/raw/${BATCH_ID}" \
            --output-dir "${TRACY_RESULTS_DIR}" \
            --ref-path "${REFERENCE}" >> "${LOG_FILE}" 2>&1

        AUTOMATE_RAW_JSON_DIR="${TRACY_RESULTS_DIR}/json"
        run_regenerate "${AUTOMATE_RAW_JSON_DIR}" "${AUTOMATE_REGENERATE_DIR}"
    fi
else
    log_message "Skipping Step 3: Automate analysis"
fi

# =============================================================================
# Step 4: Manual — Manual pipeline + Sequencher comparison
# =============================================================================

if [ "${RUN_MANUAL}" == "true" ]; then
    log_message "[Step 4/6 - manual] Extracting manual data and running Sequencher comparison"

    # Clean up and create manual directory (only if not using custom directory)
    if [ -n "${CUSTOM_MANUAL_DIR}" ]; then
        MANUAL_DIR="${CUSTOM_MANUAL_DIR}"
    else
        if [ -d "${MANUAL_DIR}" ]; then
            rm -rf "${MANUAL_DIR}"
        fi
        mkdir -p "${MANUAL_DIR}"

        # Find all zip files containing the batch number
        MANUAL_ZIP_FILES=$(find "${SEQUENCHER_DIR}" -type f \
            \( -name "${BATCH_ID}.zip" -o -name "${BATCH_ZIP_GLOB}" \))

        if [ -n "${MANUAL_ZIP_FILES}" ]; then
            log_message "Found the following zip files:"
            echo "${MANUAL_ZIP_FILES}" | while read -r zip_file; do
                log_message "- $(basename "${zip_file}")"
                unzip -q -o "${zip_file}" -d "${MANUAL_DIR}" > /dev/null 2>&1
            done
        else
            log_message "Warning: No zip files found containing batch pattern ${BATCH} in ${SEQUENCHER_DIR}"
        fi
    fi

    # Run Sequencher ETL pipeline
    if [ -d "${MANUAL_DIR}" ] && [ "$(find "${MANUAL_DIR}" -name "*.TXT" | wc -l)" -gt 0 ]; then
        SEQUENCHER_OUTPUT_DIR="${SEQUENCHER_RESULTS_DIR}/${BATCH_ID}"

        if [ -d "${SEQUENCHER_OUTPUT_DIR}" ]; then
            log_message "Cleaning up existing Sequencher output"
            rm -rf "${SEQUENCHER_OUTPUT_DIR:?}"/*
        fi
        mkdir -p "${SEQUENCHER_OUTPUT_DIR}"

        python -m src.tools.sequencher.pipeline \
            --batch-id "${BATCH_ID}" \
            --input-dir "${MANUAL_DIR}" \
            --output-dir "${SEQUENCHER_OUTPUT_DIR}" \
            --ref-path "${REFERENCE}" >> "${LOG_FILE}" 2>&1

        run_regenerate "${SEQUENCHER_OUTPUT_DIR}/json" "${MANUAL_REGENERATE_DIR}"
    else
        log_message "Skipping Sequencher ETL - no TXT files found in ${MANUAL_DIR}"
    fi

    # Comparison reads the regenerated automated JSON so reviewers compare the
    # same filtered final calls used by downstream FASTA/report generation.
    python -m src.core.comparison \
        -a "${AUTOMATE_REGENERATE_DIR}/statistic_fullbatch.json" \
        -b "${SEQUENCHER_OUTPUT_DIR}/json/statistic_fullbatch.json" \
        -o "${COMPARISON_DIR}" \
        --batch-id "${BATCH_ID}" \
        --legacy-format >> "${LOG_FILE}" 2>&1
else
    log_message "Skipping Step 4: Manual pipeline comparison and Sequencher analysis"
fi

# =============================================================================
# Step 5: FASTA — Generate FASTA files
# =============================================================================

if [ "${RUN_FASTA}" == "true" ]; then
    log_message "[Step 5/6 - fasta] Generating FASTA files"

    # Check if required files exist
    if [ ! -d "${MANUAL_REGENERATE_DIR}" ] || ! ls "${MANUAL_REGENERATE_DIR}/"*.json 1> /dev/null 2>&1; then
        log_message "Warning: No Sequencher JSON results found in ${MANUAL_REGENERATE_DIR}"
        log_message "Make sure to run manual step first or provide the required files manually."
    fi

    python -m src.generate_fasta \
        -i "${MANUAL_REGENERATE_DIR}" \
        -o "${FASTA_DIR}" \
        -r "${REFERENCE}" \
        --results-dir "${RESULTS_DIR}" \
        ${FORCE_MERGE:+--force} >> "${LOG_FILE}" 2>&1

    # Check how many FASTA files were generated
    FASTA_COUNT=$(find "${FASTA_DIR}" -name "*.fasta" | wc -l)
    log_message "Generated ${FASTA_COUNT} FASTA files in ${FASTA_DIR}"
else
    log_message "Skipping Step 5: Generate FASTA files"
fi

# =============================================================================
# Step 6: Report — Generate PDF reports
# =============================================================================

if [ "${RUN_REPORT}" == "true" ]; then
    log_message "[Step 6/6 - report] Generating PDF reports"

    # Check if required files exist
    if [ ! -d "${MANUAL_REGENERATE_DIR}" ] || ! ls "${MANUAL_REGENERATE_DIR}/"*.json 1> /dev/null 2>&1; then
        log_message "Warning: No Sequencher JSON results found in ${MANUAL_REGENERATE_DIR}"
        log_message "Make sure to run manual step first or provide the required files manually."
    fi

    python -m src.generate_reports \
        -b "${BATCH_ID}" \
        -d "${DATA_DIR}" \
        -r "${MANUAL_REGENERATE_DIR}" \
        --report_dir "${REPORTS_DIR}" >> "${LOG_FILE}" 2>&1
else
    log_message "Skipping Step 6: Generate PDF reports"
fi

log_message "Pipeline completed successfully (mode: ${PIPELINE_MODE})"
