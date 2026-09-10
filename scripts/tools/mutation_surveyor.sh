#!/bin/bash

# =============================================================================
# Mutation Surveyor Pipeline — Step-Based Shell Script
#
# Runs 10 steps (including TXT→XLSX conversion and unify) from src/tools/mutation_surveyor/ in dependency order.
# Input files are discovered by content-based globbing in IN_DIR.
#
# Supports step selection via -s (numbers or names), like the main pipeline.sh.
# =============================================================================

# Define usage and help function
usage() {
    echo "Usage: $0 [options] <BATCH_ID> <IN_DIR> <OUT_DIR> <RAW_DIR>"
    echo ""
    echo "Options:"
    echo "  -h, --help              Display this help message"
    echo "  -s, --steps STEPS       Steps to run (comma-separated, numbers or names)"
    echo "      --truthset PATH    Override truthset file (default: temp/tools/sequencher/merged/Truthset.xlsx)"
    echo ""
    echo "Steps (use numbers or names with -s):"
    echo "  1 | convert    - Convert TXT custom report to XLSX (if needed)"
    echo "  2 | trace      - Trace processing (inventory QC + trace process)"
    echo "  3 | control    - Control QC (positive/NTC control validation)"
    echo "  4 | hv23       - HV2-3 merge"
    echo "  5 | hv1        - HV1 merge"
    echo "  6 | profiles   - Final profiles merge (HV1 + HV2-3)"
    echo "  7 | truth      - Final profiles vs truthset"
    echo "  8 | etl        - MS ETL + Sequencher ETL + region filter"
    echo "  9 | unify      - Merge MS + Sequencher into unified profiles"
    echo "  10 | compare   - Compare regenerate vs unified profiles"
    echo ""
    echo "Examples:"
    echo "  $0 MS_210426_003 /in /out /raw                          # Run all steps"
    echo "  $0 -s 1,2,3,4 MS_210426_003 /in /out /raw              # Run steps 1-4 only"
    echo "  $0 -s trace,control,hv23 MS_210426_003 /in /out /raw   # Same, using names"
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

RUN_CONVERT=true
RUN_TRACE=true
RUN_CONTROL=true
RUN_HV23=true
RUN_HV1=true
RUN_PROFILES=true
RUN_TRUTH=true
RUN_ETL=true
RUN_UNIFY=true
RUN_COMPARE=true

BATCH_ID=""
IN_DIR=""
OUT_DIR=""
RAW_DIR=""

# Optional overrides (empty = use defaults computed later)
OPT_TRUTHSET=""

# =============================================================================
# Argument parsing
# =============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -s|--steps)
            # Reset all steps to false first
            RUN_CONVERT=false
            RUN_TRACE=false
            RUN_CONTROL=false
            RUN_HV23=false
            RUN_HV1=false
            RUN_PROFILES=false
            RUN_TRUTH=false
            RUN_ETL=false
            RUN_UNIFY=false
            RUN_COMPARE=false

            # Parse comma-separated list of steps (numbers or names)
            OLD_IFS="$IFS"
            IFS=',' read -ra STEPS <<< "$2"
            for step in "${STEPS[@]}"; do
                case "$step" in
                    1|convert)   RUN_CONVERT=true ;;
                    2|trace)     RUN_TRACE=true ;;
                    3|control)   RUN_CONTROL=true ;;
                    4|hv23)      RUN_HV23=true ;;
                    5|hv1)       RUN_HV1=true ;;
                    6|profiles)  RUN_PROFILES=true ;;
                    7|truth)     RUN_TRUTH=true ;;
                    8|etl)       RUN_ETL=true ;;
                    9|unify)     RUN_UNIFY=true ;;
                    10|compare)  RUN_COMPARE=true ;;
                    *)           echo "Invalid step: $step"; usage ;;
                esac
            done
            IFS="$OLD_IFS"
            shift 2
            ;;
        --truthset)
            OPT_TRUTHSET="$2"
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
            # Positional arguments: BATCH_ID IN_DIR OUT_DIR RAW_DIR
            if [ -z "$BATCH_ID" ]; then
                BATCH_ID="$1"
            elif [ -z "$IN_DIR" ]; then
                IN_DIR="$1"
            elif [ -z "$OUT_DIR" ]; then
                OUT_DIR="$1"
            elif [ -z "$RAW_DIR" ]; then
                RAW_DIR="$1"
            else
                echo "Error: Too many positional arguments"
                usage
            fi
            shift
            ;;
    esac
done

# =============================================================================
# Validate required positional arguments
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
if [ -z "$RAW_DIR" ]; then
    echo "Error: RAW_DIR is required"
    usage
fi

# --- Validate input directories ---
if [ ! -d "${IN_DIR}" ]; then
    echo "Error: IN_DIR does not exist: ${IN_DIR}"
    exit 1
fi
if [ ! -d "${RAW_DIR}" ]; then
    echo "Error: RAW_DIR does not exist: ${RAW_DIR}"
    exit 1
fi

mkdir -p "${OUT_DIR}"

# JSON temp directory for per-region intermediate files (MS + Sequencher)
# Derived from the analysis base directory to match the shared json_temp structure:
#   analysis/json_temp/${BATCH_ID}/mutation_surveyor/HV1/HV1_{LID}.json (and HV2-3/)
#   analysis/json_temp/${BATCH_ID}/sequencher/HV1/HV1_{LID}.json (and HV2-3/)
ANALYSIS_DIR="$(dirname "$(dirname "${OUT_DIR}")")"
JSON_TEMP_DIR="${ANALYSIS_DIR}/json_temp/${BATCH_ID}"
mkdir -p "${JSON_TEMP_DIR}/mutation_surveyor"
mkdir -p "${JSON_TEMP_DIR}/sequencher"

# =============================================================================
# Load configuration and setup logging
# =============================================================================

# shellcheck source=/dev/null
source .env

LOG_FILE="${LOGS_DIR}/${BATCH_ID}.log"
mkdir -p "$(dirname "$LOG_FILE")"

log_message "=========================================="
log_message "Mutation Surveyor pipeline"
log_message "Batch ID: $BATCH_ID"
log_message "=========================================="

# =============================================================================
# Resolve input files by content-based globbing
# =============================================================================

LID_MANIFEST="/mnt/nas/bca/mtDNA/science/data/metadata/${BATCH_ID}.xlsx"
if [ ! -f "${LID_MANIFEST}" ]; then
    LID_MANIFEST=$(ls "${IN_DIR}"/*Danh*.xlsx 2>/dev/null | head -1)
fi

# PC reference: prefer batch-local copy, fall back to module-bundled reference
CONTROL_REF="${IN_DIR}/mtDNA_PC_reference.xlsx"
if [ ! -f "${CONTROL_REF}" ]; then
    CONTROL_REF="src/tools/mutation_surveyor/01_02/mtDNA_PC_reference.xlsx"
fi

TRUTHSET="temp/tools/sequencher/merged/Truthset.xlsx"
# Apply optional override
if [ -n "${OPT_TRUTHSET}" ]; then
    TRUTHSET="${OPT_TRUTHSET}"
fi

# =============================================================================
# File check helper
# =============================================================================

check_file() {
    local label="$1"
    local path="$2"
    if [ ! -f "$path" ]; then
        echo "Error: Missing ${label}: ${path}"
        exit 1
    fi
}

# =============================================================================
# Step execution
# =============================================================================

# --- Step 1: Convert TXT custom report to XLSX (if needed) ---

# Priority-ordered glob: df > default > CR_mutationSurveyor > bare
# When batch is provided, only match files starting with {batch}_ to avoid
# picking up wrong-batch custom reports from a shared IN_DIR.
find_custom_report() {
    local dir="$1" ext="$2" batch="${3:-}"
    for stem in "df_custom_report" "default_custom_report" "CR_mutationSurveyor_CRSetting" "_custom_report"; do
        local hit
        if [ -n "${batch}" ]; then
            hit=$(ls "${dir}"/${batch}_*${stem}*.${ext} 2>/dev/null | head -1)
        else
            hit=$(ls "${dir}"/*${stem}*.${ext} 2>/dev/null | head -1)
        fi
        if [ -n "${hit}" ]; then
            echo "${hit}"
            return
        fi
    done
}

if [ "$RUN_CONVERT" == "true" ]; then
    log_message "[Step 1] Convert TXT to XLSX"

    MS_FILE=$(find_custom_report "${IN_DIR}" xlsx "${BATCH_ID}")

    if [ -n "${MS_FILE}" ]; then
        log_message "XLSX custom report found: ${MS_FILE}"
    else
        TXT_FILE=$(find_custom_report "${IN_DIR}" txt "${BATCH_ID}")

        if [ -n "${TXT_FILE}" ]; then
            log_message "No XLSX custom report found. Converting TXT: ${TXT_FILE}"

            python -m src.tools.mutation_surveyor.txt_to_excel \
                --input "${TXT_FILE}" >> "$LOG_FILE" 2>&1

            MS_FILE=$(ls "${IN_DIR}"/${BATCH_ID}_*df_custom_report*.xlsx 2>/dev/null | head -1)
            if [ ! -f "${MS_FILE}" ]; then
                echo "Error: TXT→XLSX conversion failed — expected output not found in ${IN_DIR}"
                exit 1
            fi
            log_message "Converted XLSX: ${MS_FILE}"
        else
            echo "Error: No XLSX or TXT custom report found in ${IN_DIR}"
            exit 1
        fi
    fi
else
    log_message "Skipping Step 1"
    MS_FILE=$(find_custom_report "${IN_DIR}" xlsx "${BATCH_ID}")
fi

# --- Step 2: Trace processing (inventory QC + trace process) ---

if [ "$RUN_TRACE" == "true" ]; then
    log_message "[Step 2] Trace processing"

    check_file "MS custom report" "${MS_FILE}"
    check_file "LID manifest" "${LID_MANIFEST}"

    python -m src.tools.mutation_surveyor.01_02.ms_trace_process \
        --ms "${MS_FILE}" \
        --lid-manifest "${LID_MANIFEST}" \
        --inventory-out "${OUT_DIR}/${BATCH_ID}_MS_trace_processed_inventory_QC.xlsx" \
        --out "${OUT_DIR}/${BATCH_ID}_MS_trace_processed.xlsx" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 2"
fi

# --- Step 3: Control QC (positive/NTC control validation) ---

if [ "$RUN_CONTROL" == "true" ]; then
    log_message "[Step 3] Control QC"

    TRACE_PROCESSED="${OUT_DIR}/${BATCH_ID}_MS_trace_processed.xlsx"
    if [ ! -f "${TRACE_PROCESSED}" ]; then
        echo "Error: Missing trace-processed file: ${TRACE_PROCESSED} (run step 2 first)"
        exit 1
    fi
    check_file "PC reference" "${CONTROL_REF}"

    python -m src.tools.mutation_surveyor.01_02.ms_control_QC \
        --trace-processed "${TRACE_PROCESSED}" \
        --control-ref "${CONTROL_REF}" \
        --out "${OUT_DIR}/control_QC/${BATCH_ID}_control_QC.xlsx" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 3"
fi

# --- Step 4: HV2-3 merge ---

if [ "$RUN_HV23" == "true" ]; then
    log_message "[Step 4] HV2-3 merge"

    TRACE_PROCESSED="${OUT_DIR}/${BATCH_ID}_MS_trace_processed.xlsx"
    if [ ! -f "${TRACE_PROCESSED}" ]; then
        echo "Error: Missing trace-processed file: ${TRACE_PROCESSED} (run step 2 first)"
        exit 1
    fi

    python -m src.tools.mutation_surveyor.03_04.ms_hv23_merge \
        --input "${TRACE_PROCESSED}" \
        --output "${OUT_DIR}/${BATCH_ID}_HV23_merged.xlsx" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 4"
fi

# --- Step 5: HV1 merge ---

if [ "$RUN_HV1" == "true" ]; then
    log_message "[Step 5] HV1 merge"

    TRACE_PROCESSED="${OUT_DIR}/${BATCH_ID}_MS_trace_processed.xlsx"
    if [ ! -f "${TRACE_PROCESSED}" ]; then
        echo "Error: Missing trace-processed file: ${TRACE_PROCESSED} (run step 2 first)"
        exit 1
    fi

    python -m src.tools.mutation_surveyor.03_04.ms_hv1_merge \
        --input "${TRACE_PROCESSED}" \
        --output "${OUT_DIR}/${BATCH_ID}_HV1_merged.xlsx" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 5"
fi

# --- Step 6: Final profiles merge (HV1 + HV2-3) ---

if [ "$RUN_PROFILES" == "true" ]; then
    log_message "[Step 6] Final profiles merge"

    HV23="${OUT_DIR}/${BATCH_ID}_HV23_merged.xlsx"
    HV1="${OUT_DIR}/${BATCH_ID}_HV1_merged.xlsx"
    if [ ! -f "${HV23}" ]; then
        echo "Error: Missing HV2-3 merged file: ${HV23} (run step 4 first)"
        exit 1
    fi
    if [ ! -f "${HV1}" ]; then
        echo "Error: Missing HV1 merged file: ${HV1} (run step 5 first)"
        exit 1
    fi

    python -m src.tools.mutation_surveyor.05_06.ms_final_profiles_merge \
        --hv23 "${HV23}" \
        --hv1 "${HV1}" \
        --output "${OUT_DIR}/${BATCH_ID}_final_profiles.xlsx" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 6"
fi

# --- Step 7: Final profiles vs truth ---

if [ "$RUN_TRUTH" == "true" ]; then
    log_message "[Step 7] Final profiles vs truthset"

    FINAL_PROFILES="${OUT_DIR}/${BATCH_ID}_final_profiles.xlsx"
    if [ ! -f "${FINAL_PROFILES}" ]; then
        echo "Error: Missing final profiles file: ${FINAL_PROFILES} (run step 6 first)"
        exit 1
    fi
    check_file "Truthset" "${TRUTHSET}"

    python -m src.tools.mutation_surveyor.05_06.ms_final_profiles_vs_truth \
        --merged "${FINAL_PROFILES}" \
        --truth "${TRUTHSET}" \
        --output "${OUT_DIR}/${BATCH_ID}_final_profiles_vs_truth.xlsx" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 7"
fi

# =============================================================================
# Step 8: MS ETL (per-region JSON output) + Sequencher ETL (via sequencher.sh)
# =============================================================================
# MS ETL produces per-region intermediate JSONs (HV1_{LID}.json, HV2-3_{LID}.json)
# for auto-pass regions. Sequencher ETL produces per-region JSONs from TXT files
# (1 or 2 per sample, depending on whether TXTs are full-region or region-specific).
# Both are consumed by Step 9 (unify) for region-based merge with region_sources
# provenance tracking.

if [ "$RUN_ETL" == "true" ]; then
    log_message "[Step 8] MS ETL + Sequencher ETL"

    FINAL_PROFILES="${OUT_DIR}/${BATCH_ID}_final_profiles.xlsx"
    if [ ! -f "${FINAL_PROFILES}" ]; then
        echo "Error: Missing final profiles file: ${FINAL_PROFILES} (run step 6 first)"
        exit 1
    fi

    python -m src.tools.mutation_surveyor.pipeline \
        --input-dir "${FINAL_PROFILES}" \
        --output-dir "${RESULTS_DIR}/tools/mutation_surveyor" \
        --json-dir "${JSON_TEMP_DIR}/mutation_surveyor" \
        --batch-id "${BATCH_ID}" >> "$LOG_FILE" 2>&1

    SEQUENCHER_IN_DIR="${SEQUENCHER_DIR}/${BATCH_ID}/txt"
    SEQUENCHER_JSON_DIR="${JSON_TEMP_DIR}/sequencher"

    # --- Sequencher ETL: delegate to sequencher.sh ---
    bash "scripts/tools/sequencher.sh" \
    "${BATCH_ID}" \
    "${SEQUENCHER_IN_DIR}" \
    "${RESULTS_DIR}/tools/sequencher" \
    --json-dir "${JSON_TEMP_DIR}/sequencher"
else
    log_message "Skipping Step 8"
fi

# =============================================================================
# Step 9: Unify MS + Sequencher per-region JSONs into unified profiles
# =============================================================================
# Discovers per-region JSONs (HV1_{LID}.json, HV2-3_{LID}.json) from both MS and
# Sequencher directories. Merges by region via merge_by_regions() with
# region_sources provenance tracking. Falls back to whole-sample merge for
# legacy combined JSONs. Produces 1 unified {LID}.json per sample with
# source_tool="unified" and region_sources in information.

if [ "$RUN_UNIFY" == "true" ]; then
    log_message "[Step 9] Unify MS + Sequencher profiles"

    MS_JSON_DIR="${JSON_TEMP_DIR}/mutation_surveyor"
    SEQ_JSON_DIR="${JSON_TEMP_DIR}/sequencher"

    if [ ! -d "${MS_JSON_DIR}" ]; then
        echo "Error: Missing MS JSON directory: ${MS_JSON_DIR} (run step 8 first)"
        exit 1
    fi

    if [ ! -d "${SEQ_JSON_DIR}" ]; then
        log_message "Skipping unify: Sequencher JSON directory not found (${SEQ_JSON_DIR})"
    else
        python -m src.core.unify \
            --ms-dir "${MS_JSON_DIR}" \
            --seq-dir "${SEQ_JSON_DIR}" \
            --output-dir "${JSON_TEMP_DIR}/unified" \
            --batch-id "${BATCH_ID}" >> "$LOG_FILE" 2>&1
    fi
else
    log_message "Skipping Step 9"
fi

# =============================================================================
# Step 10: Compare regenerate vs unified profiles
# =============================================================================

if [ "$RUN_COMPARE" == "true" ]; then
    log_message "[Step 10] Compare regenerate vs unified profiles"

    REGEN_BATCH="results/automate_pipeline/regenerate/${BATCH_ID}/statistic_fullbatch.json"
    UNIFIED_BATCH="${JSON_TEMP_DIR}/unified/${BATCH_ID}/statistic_fullbatch.json"

    if [ ! -f "${REGEN_BATCH}" ]; then
        echo "Error: Missing regenerate batch JSON: ${REGEN_BATCH} (run automate pipeline first)"
        exit 1
    fi
    if [ ! -f "${UNIFIED_BATCH}" ]; then
        echo "Error: Missing unified batch JSON: ${UNIFIED_BATCH} (run step 9 first)"
        exit 1
    fi

    python -m src.core.comparison \
        --batch-a "${REGEN_BATCH}" \
        --batch-b "${UNIFIED_BATCH}" \
        --output-dir "results/modules/comparison/${BATCH_ID}" \
        --batch-id "${BATCH_ID}" >> "$LOG_FILE" 2>&1
else
    log_message "Skipping Step 10"
fi

# =============================================================================
# Done
# =============================================================================

log_message "Mutation Surveyor pipeline completed for ${BATCH_ID}"
