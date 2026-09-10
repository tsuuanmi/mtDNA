#!/bin/bash

set -euo pipefail

# =============================================================================
# Batch Pipeline Script
# Supports two modes: normal and rerun
# =============================================================================

# Load environment variables
# shellcheck source=/dev/null
source .env

# Create logs directory at the very beginning
mkdir -p "${LOGS_DIR}"

# =============================================================================
# Logging
# =============================================================================

log_message() {
    local message="$1"
    local batch_id="${2:-}"

    if [ -n "$batch_id" ]; then
        local log_file="${LOGS_DIR}/${batch_id}.log"
        python -c "from loguru import logger; logger.info(\"$message\")" | tee -a "$log_file"
    else
        python -c "from loguru import logger; logger.info(\"$message\")"
    fi
}

# =============================================================================
# Setup directories
# =============================================================================

# Drop the empty placeholder dirs a relative-path run leaves in the repo root.
# Kept to one level deep on purpose: the recursive form walked the whole tree
# and would delete an output directory another run had just created but not
# yet written into.
find . -maxdepth 1 -type d -empty -delete

# # Create results directory if it doesn't exist
# rm -rf "${RESULTS_DIR}/tools/sequencher"
# mkdir -p "${RESULTS_DIR}/tools/sequencher"

# =============================================================================
# Batch configuration — defaults, all overridable from the command line
# =============================================================================

# Analysis engine: tracy | blastn
PIPELINE="tracy"

# Pipeline steps forwarded to scripts/pipeline.sh -s
STEPS="3,4"

# Batches to process when --batches is not supplied.
BATCHES=()

# =============================================================================
# Parse arguments
# =============================================================================

usage() {
    cat <<'USAGE'
Usage: bash scripts/batch_pipeline.sh [OPTIONS]

Options:
  -h, --help                Show this help message
  -r, --rerun               RERUN mode: pick samples from Sequencher ZIPs
  -s, --steps STEPS         Pipeline steps, comma-separated       (default: 3,4)
  -p, --pipeline MODE       Analysis engine: tracy | blastn       (default: tracy)
  -b, --batches LIST        Batch IDs, comma-separated
  -j, --jobs N              Batches to process at once (default: 1)
      --max-cores N         Cores this run may use, split between its batches
      --sequencher-dir DIR  Override SEQUENCHER_DIR from .env
      --results-dir DIR     Override RESULTS_DIR from .env (created if missing)
      --sheet-id ID         Override GOOGLE_SHEET_ID from .env
      --upload              Force Google Sheets upload on
      --no-upload           Force Google Sheets upload off
                            (default: on for normal mode, OFF for --rerun)

Examples:
  # Automate + compare
  bash scripts/batch_pipeline.sh -s 2,3,4 -p tracy -b MS_110826_012

  # Re-compare only, after manual results changed
  bash scripts/batch_pipeline.sh -s 4 -b MS_110826_012

  # Rerun, writing to a private results dir
  bash scripts/batch_pipeline.sh --rerun -s 3,4 -p blastn \
      --sequencher-dir /mnt/nas/bca/mtDNA/Sequencher_temp/rerun/rerun_260425_X \
      --results-dir /home/minhtq/rerun_results
USAGE
    exit 0
}

RERUN_MODE=false
JOBS=1
OPT_MAX_CORES=""
OPT_SEQUENCHER_DIR=""
OPT_RESULTS_DIR=""
OPT_SHEET_ID=""
UPLOAD=""   # empty = auto: on for normal mode, off for rerun

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)         usage ;;
        -r|--rerun)        RERUN_MODE=true; shift ;;
        -s|--steps)        STEPS="$2"; shift 2 ;;
        -p|--pipeline)     PIPELINE="$2"; shift 2 ;;
        -b|--batches)      IFS=',' read -ra BATCHES <<< "$2"; shift 2 ;;
        -j|--jobs)         JOBS="$2"; shift 2 ;;
        --max-cores)       OPT_MAX_CORES="$2"; shift 2 ;;
        --sequencher-dir)  OPT_SEQUENCHER_DIR="$2"; shift 2 ;;
        --results-dir)     OPT_RESULTS_DIR="$2"; shift 2 ;;
        --sheet-id)        OPT_SHEET_ID="$2"; shift 2 ;;
        --upload)          UPLOAD=true; shift ;;
        --no-upload)       UPLOAD=false; shift ;;
        *) echo "Unknown option: $1" >&2; echo "Run with --help for usage." >&2; exit 1 ;;
    esac
done

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [ "$JOBS" -lt 1 ]; then
    echo "Error: --jobs must be a positive integer, got '$JOBS'" >&2
    exit 1
fi

if [[ "$PIPELINE" != "tracy" && "$PIPELINE" != "blastn" ]]; then
    echo "Error: --pipeline must be 'tracy' or 'blastn', got '$PIPELINE'" >&2
    exit 1
fi

if [ ${#BATCHES[@]} -eq 0 ]; then
    echo "Error: no batches to process. Pass -b/--batches BATCH1,BATCH2" >&2
    echo "Run with --help for usage." >&2
    exit 1
fi

# Overrides must be applied here: `source .env` above reassigns these
# unconditionally, so an exported environment variable would be discarded.
if [ -n "$OPT_SEQUENCHER_DIR" ]; then
    SEQUENCHER_DIR="$OPT_SEQUENCHER_DIR"
fi
if [ -n "$OPT_RESULTS_DIR" ]; then
    RESULTS_DIR="$OPT_RESULTS_DIR"
fi
if [ -n "$OPT_SHEET_ID" ]; then
    GOOGLE_SHEET_ID="$OPT_SHEET_ID"
fi

# Rerun writes to a private RESULTS_DIR and must never replace the shared
# production worksheet, so upload stays off there unless asked for explicitly.
if [ -z "$UPLOAD" ]; then
    if [[ "$RERUN_MODE" == "true" ]]; then
        UPLOAD=false
    else
        UPLOAD=true
    fi
fi

METADATA_DIR="${DATA_DIR}/metadata"
mkdir -p "$RESULTS_DIR"

log_message "Config: pipeline=${PIPELINE} steps=${STEPS} rerun=${RERUN_MODE} upload=${UPLOAD}"
log_message "Config: batches=${BATCHES[*]}"
log_message "Config: SEQUENCHER_DIR=${SEQUENCHER_DIR}"
log_message "Config: RESULTS_DIR=${RESULTS_DIR}"

# =============================================================================
# Phase 1: Extract ZIPs and identify rerun samples (rerun mode only)
# =============================================================================

if [[ "$RERUN_MODE" == "true" ]]; then
    log_message "Running in RERUN mode"

    RERUN_LISTS_DIR="${RESULTS_DIR}/rerun_lists"
    mkdir -p "$RERUN_LISTS_DIR"

    # Refresh only the raw sample TXT needed by rerun matching. This avoids
    # running the full prepare step while ensuring the TXT contains all metadata
    # LID values except PC/NTC controls.
    refresh_raw_sample_txt() {
        local batch="$1"
        local metadata_file="${METADATA_DIR}/${batch}.xlsx"
        local raw_list="${DATA_DIR}/raw/${batch}.txt"

        if [ ! -f "$metadata_file" ]; then
            log_message "ERROR: Metadata Excel not found for batch $batch: $metadata_file" "$batch"
            return 1
        fi

        mkdir -p "$(dirname "$raw_list")"
        log_message "Extracting sample IDs from Excel file (excluding PC and NTC controls)" "$batch"
        # `if ! python` rather than a bare call: under `set -e` a bad workbook
        # would abort the whole run and take the remaining batches with it.
        if ! python - "$metadata_file" "$raw_list" <<'PY'
import sys
from pathlib import Path

import pandas as pd

metadata_path = Path(sys.argv[1])
raw_list_path = Path(sys.argv[2])
df = pd.read_excel(metadata_path)

if "LID" not in df.columns:
    print(f"no LID column in {metadata_path}", file=sys.stderr)
    raise SystemExit(1)

lid_series = df["LID"].dropna().drop_duplicates().astype(str)
filtered = lid_series[~lid_series.str.contains("PC|NTC", case=False, na=False)]

if filtered.empty:
    print(f"no sample IDs left after filtering {metadata_path}", file=sys.stderr)
    raise SystemExit(1)

# Write through a temp file. raw_list is the batch's canonical sample list on
# shared storage; truncating it and then failing would leave the batch with an
# empty list and nothing to say why.
tmp_path = raw_list_path.with_name(raw_list_path.name + ".tmp")
filtered.to_csv(tmp_path, index=False, header=False)
tmp_path.replace(raw_list_path)
PY
        then
            log_message "ERROR: Could not extract sample IDs for batch $batch; leaving $raw_list as it was" "$batch"
            return 1
        fi
        return 0
    }

    # A rerun's unit of work is a (batch, source folder) pair, not a batch. The
    # same batch often sits in two rerun folders as two different Sequencher
    # analyses of the same samples, and pipeline.sh unzips every ZIP matching
    # the batch into one directory: identical TXT names meant the second
    # analysis silently overwrote the first, and differing names tripped the
    # duplicate-LID guard and dropped the sample -- or the whole batch. Running
    # each folder on its own, into its own results tree, keeps the two apart.
    PAIRS_FILE="${RERUN_LISTS_DIR}/pairs.tsv"
    : > "$PAIRS_FILE"

    # The folder basename names both a results subtree and a worksheet, and the
    # latter reaches the Sheets API: allowed characters only, capped at the
    # 100-character worksheet title limit. Parameter expansion rather than
    # basename|tr, which would turn the trailing newline into an underscore.
    sanitise_tab() {
        local name="${1%/}"
        name="${name##*/}"
        name="${name//[^A-Za-z0-9_.-]/_}"
        printf '%s' "${name:0:100}"
    }

    extract_rerun_samples() {
        local batch="$1"
        local folders_file="${RERUN_LISTS_DIR}/${batch}_folders.txt"
        local zip zip_parent folder tab

        # Anchor on the full batch id and require a separator after it. The old
        # "*${batch#*_}*" form matched as a substring, so batch 20241222_mtDNA_12
        # also picked up mtDNA_120 and mtDNA_121 zips.
        mapfile -t zip_files < <(find "$SEQUENCHER_DIR" -type f \
            \( -name "${batch}.zip" -o -name "${batch}[_. ]*.zip" \))

        if [ ${#zip_files[@]} -eq 0 ]; then
            log_message "WARNING: No ZIP files found for batch $batch" "$batch"
            return
        fi

        for zip in "${zip_files[@]}"; do
            zip_parent="${zip%/*}"
            printf '%s\n' "${zip_parent%/}"
        done | sort -u > "$folders_file"

        log_message "Batch $batch has ZIPs in $(wc -l < "$folders_file") folder(s)" "$batch"
        while IFS= read -r folder || [ -n "$folder" ]; do
            [ -n "$folder" ] || continue
            tab=$(sanitise_tab "$folder")
            extract_rerun_samples_for_folder "$batch" "$folder" "$tab"
        done < "$folders_file"
        rm -f "$folders_file"
    }

    # Same selection as before -- ZIP member names to sample ids, filtered
    # against the batch's raw TXT -- but restricted to one folder's ZIPs.
    extract_rerun_samples_for_folder() {
        local batch="$1"
        local folder="$2"
        local tab="$3"
        local dir="${RERUN_LISTS_DIR}/${tab}"
        mkdir -p "$dir"
        local rerun_list="${dir}/${batch}_rerun.txt"
        local missing_txt_list="${dir}/${batch}_zip_not_in_txt.txt"
        local zip_entries="${dir}/${batch}_zip_entries.txt"
        local zip_sample_entries="${dir}/${batch}_zip_sample_entries.tsv"
        local zip_samples="${dir}/${batch}_zip_samples.txt"
        local raw_samples="${dir}/${batch}_raw_samples.txt"

        # -maxdepth 1: only this folder's own ZIPs, never a nested folder's.
        mapfile -t zip_files < <(find "$folder" -maxdepth 1 -type f \
            \( -name "${batch}.zip" -o -name "${batch}[_. ]*.zip" \))

        # Extract sample IDs from ZIP member names without unpacking files.
        # The pipeline only needs the rerun sample list here; extracted TXT
        # contents are not used and some ZIPs preserve read-only directories.
        : > "$zip_entries"
        for zip in "${zip_files[@]}"; do
            if ! zipinfo -1 "$zip" | while IFS= read -r zip_entry; do
                [[ "$zip_entry" == *.TXT ]] || continue
                basename "$zip_entry"
            done >> "$zip_entries"; then
                log_message "ERROR: Failed to inspect $zip for batch $batch" "$batch"
                continue
            fi
        done
        sort -u -o "$zip_entries" "$zip_entries"

        python - "$zip_entries" "$zip_sample_entries" <<'PY'
import sys
from pathlib import Path

from src.tools.sequencher.quality_control import parse_filename_ranges

entries_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
with entries_path.open() as entries, output_path.open("w") as output:
    for line in entries:
        zip_entry = line.rstrip("\n")
        if not zip_entry:
            continue
        sample_id = parse_filename_ranges(zip_entry).sample_id
        output.write(f"{sample_id}\t{zip_entry}\n")
PY
        cut -f1 "$zip_sample_entries" | sort -u > "$zip_samples"

        # Filter samples
        local raw_list="${DATA_DIR}/raw/${batch}.txt"
        if [ -f "$raw_list" ]; then
            awk '{ sub(/\r$/, ""); if (NF) print }' "$raw_list" | sort -u > "$raw_samples"
            : > "$rerun_list"
            while IFS= read -r sample || [ -n "$sample" ]; do
                sample="${sample%$'\r'}"
                [[ -z "$sample" ]] && continue
                if grep -qxF "$sample" "$zip_samples" 2>/dev/null; then
                    echo "$sample" >> "$rerun_list"
                fi
            done < "$raw_list"
        else
            log_message "WARNING: Raw sample TXT not found for batch $batch: $raw_list" "$batch"
            : > "$rerun_list"
            : > "$raw_samples"
        fi

        local zip_sample_count
        local found_sample_count
        zip_sample_count=$(awk 'NF { count++ } END { print count + 0 }' "$zip_samples")
        found_sample_count=$(awk 'NF { count++ } END { print count + 0 }' "$rerun_list")
        log_message "Rerun ZIP sample summary for batch $batch in ${tab}: $zip_sample_count unique sample(s) in ZIP, $found_sample_count found in raw TXT" "$batch"

        : > "$missing_txt_list"
        if [ -f "$raw_list" ]; then
            while IFS=$'\t' read -r sample zip_entry || [ -n "$sample" ]; do
                [[ -z "$sample" ]] && continue
                if ! grep -qxF "$sample" "$raw_samples" 2>/dev/null; then
                    printf '%s\t%s\n' "$sample" "$zip_entry" >> "$missing_txt_list"
                    log_message "Sequencher ZIP TXT '$zip_entry' has sample '$sample' not found in raw TXT for batch $batch" "$batch"
                fi
            done < "$zip_sample_entries"
        fi

        # Cleanup
        rm -f "$zip_entries" "$zip_sample_entries" "$zip_samples" "$raw_samples"

        # Only a pair with something to run becomes a unit of work. Recording
        # the folder path as well as the tab name means Phase 2 never has to
        # rediscover it, and a rerun pointed straight at one folder produces
        # exactly one pair -- the old behaviour.
        if [ -s "$rerun_list" ]; then
            printf '%s\t%s\t%s\n' "$batch" "$tab" "$folder" >> "$PAIRS_FILE"
        else
            log_message "[$batch] No rerun samples found in ${tab}, skipping" "$batch"
        fi
    }

    # Read do lab giải trình tự lại nằm ở Lab/Raw_data/rerun và chưa bao giờ
    # được nạp vào thư mục làm việc của batch, nên tới giờ một lượt rerun vẫn
    # phân tích lại đúng read cũ -- chính cái đã hỏng nên mới phải chạy lại.
    RERUN_AB1_ROOT="${RERUN_AB1_ROOT:-${LAB_DATA_DIR%/}/rerun}"

    apply_rerun_reads() {
        local batch="$1"
        local samples="${RERUN_LISTS_DIR}/${batch}_reads_wanted.txt"

        # "Bước 2" của rerun là thay read lab đã chạy lại -- không phải bước 2
        # thường (gom AB1 từ đầu). Bỏ chọn thì rerun phân tích lại read gốc,
        # đó là hành vi cũ và đôi khi vẫn cần để đối chiếu.
        if [[ ",${STEPS}," != *",2,"* ]]; then
            log_message "Skipping rerun read swap (step 2 not selected)" "$batch"
            return 0
        fi

        # Gộp mẫu của mọi thư mục nguồn rồi thay một lần: các thư mục nguồn
        # khác nhau vẫn đọc chung một thư mục AB1 của batch.
        cat "${RERUN_LISTS_DIR}"/*/"${batch}_rerun.txt" 2>/dev/null \
            | sort -u > "$samples" || true
        if [ ! -s "$samples" ]; then
            rm -f "$samples"
            return 0
        fi
        if [ ! -d "$RERUN_AB1_ROOT" ]; then
            log_message "WARNING: No rerun AB1 directory at $RERUN_AB1_ROOT for batch $batch; using original reads" "$batch"
            rm -f "$samples"
            return 0
        fi
        # `if ! python`: không thay được read thì lượt chạy vẫn tiếp tục với
        # read gốc, nhưng phải nói rõ chứ không im lặng.
        if ! python -m src.tools.rerun.apply_reads \
            --batch-id "$batch" \
            --samples "$samples" \
            --rerun-root "$RERUN_AB1_ROOT" \
            --raw-dir "${DATA_DIR}/raw/${batch}" \
            --report "${RERUN_LISTS_DIR}/${batch}_reads_applied.json"; then
            log_message "WARNING: Could not apply rerun reads for batch $batch; using original reads" "$batch"
        fi
        rm -f "$samples"
    }

    # Process each batch
    for BATCH in "${BATCHES[@]}"; do
        if ! refresh_raw_sample_txt "$BATCH"; then
            log_message "ERROR: Cannot prepare rerun sample list for batch $BATCH" "$BATCH"
            exit 1
        fi
        extract_rerun_samples "$BATCH"
        apply_rerun_reads "$BATCH"
    done

else
    log_message "Running in NORMAL mode"
fi

# =============================================================================
# Phase 2: Run pipeline for each batch
# =============================================================================

# Forwarded to every pipeline.sh invocation so the child run uses the same
# overrides as this script, instead of re-reading the unmodified .env.
# Rerun không bao giờ chuyển bước 2 xuống pipeline.sh: ở đó bước 2 là
# "rm -rf RAW_DIR rồi nạp lại từ LAB_DATA_DIR", tức xoá sạch read vừa thay và
# mang chính read hỏng quay về. Bước 2 của rerun đã làm xong ở giai đoạn 1.
PIPELINE_STEPS="$STEPS"
if [[ "$RERUN_MODE" == "true" ]]; then
    PIPELINE_STEPS=$(printf '%s' "$STEPS" | tr ',' '\n' | grep -vx '2' | paste -sd, -)
fi

PIPELINE_ARGS=(
    -p "$PIPELINE"
    -s "$PIPELINE_STEPS"
)

# Each batch writes only to paths carrying its own id -- data/raw/<batch>,
# results/tools/*/<batch>, modules/comparison/<batch> -- so batches do not
# collide. The merge below is the shared part, and it runs after all of them.
#
# Every line of a batch's output is prefixed with its id. Without that the
# streams interleave and there is no way, in the log or in the web UI, to tell
# which batch a "[Step 3/6]" line belongs to.
STATUS_DIR=$(mktemp -d)
trap 'rm -rf "$STATUS_DIR"' EXIT

run_batch() {
    local batch="$1"
    local rc=0

    log_message "[$batch] Running pipeline (normal mode)" "$batch"
    # `sed -u`: without it sed fills a 4KB block before writing, and the web UI
    # captures this stream through a pipe. A single batch never prints 4KB, so
    # the run looked frozen from 08:58 to 09:03 while it was in fact working.
    bash scripts/pipeline.sh "${PIPELINE_ARGS[@]}" \
        -o "$RESULTS_DIR" --sequencher-dir "$SEQUENCHER_DIR" "$batch" 2>&1 \
        | sed -u "s|^|[$batch] |" || rc=$?

    echo "$rc" > "${STATUS_DIR}/${batch}"
    return 0
}

# One (batch, folder) pair. --sequencher-dir is the single folder, so step 4
# unzips only that folder's ZIP; -o is a results tree of its own, so the same
# batch coming from a second folder cannot overwrite this one's tool output,
# manual TXT or comparison.
run_pair() {
    local batch="$1"
    local tab="$2"
    local folder="$3"
    local rc=0
    local rerun_list="${RERUN_LISTS_DIR}/${tab}/${batch}_rerun.txt"
    local out="${RESULTS_DIR}/by_folder/${tab}"
    local sample_count
    sample_count=$(wc -l < "$rerun_list")
    # pipeline.sh rejects a -o that does not exist yet.
    mkdir -p "$out"

    if [ -z "$PIPELINE_STEPS" ]; then
        log_message "[$batch] Only the read swap was requested, nothing to analyse" "$batch"
        echo 0 > "${STATUS_DIR}/${tab}__${batch}"
        return 0
    fi

    log_message "[$batch] Running pipeline with $sample_count samples from ${tab} (rerun mode)" "$batch"
    bash scripts/pipeline.sh "${PIPELINE_ARGS[@]}" \
        -o "$out" --sequencher-dir "$folder" \
        -l "$rerun_list" "$batch" < /dev/null 2>&1 \
        | sed -u "s|^|[$batch] |" || rc=$?

    echo "$rc" > "${STATUS_DIR}/${tab}__${batch}"
    return 0
}

# Left to itself every tool takes os.cpu_count() workers -- all 256 here -- so
# one run saturates the machine and N concurrent batches merely contend. Cap the
# whole run instead and split that cap between the jobs. The cap applies at -j 1
# too: a single run has no business taking the entire machine when other people
# are working on it. Raise or lower with MTDNA_MAX_CORES in .env.
TOTAL_CORES=$(nproc)
# --max-cores wins over the .env default: the web UI hands out shares of a
# machine-wide budget, so it knows better than a fixed number what this run
# should take while other people are running.
MAX_CORES="${OPT_MAX_CORES:-${MTDNA_MAX_CORES:-128}}"
if ! [[ "$MAX_CORES" =~ ^[0-9]+$ ]] || [ "$MAX_CORES" -lt 1 ]; then
    echo "Error: MTDNA_MAX_CORES must be a positive integer, got '$MAX_CORES'" >&2
    exit 1
fi
[ "$MAX_CORES" -gt "$TOTAL_CORES" ] && MAX_CORES="$TOTAL_CORES"

WORKERS_PER_JOB=$(( MAX_CORES / JOBS ))
[ "$WORKERS_PER_JOB" -lt 1 ] && WORKERS_PER_JOB=1
export MTDNA_TRACY_MAX_WORKERS="$WORKERS_PER_JOB"
export MTDNA_BLASTN_MAX_WORKERS="$WORKERS_PER_JOB"
export MTDNA_SEQUENCHER_MAX_WORKERS="$WORKERS_PER_JOB"
log_message "Up to $JOBS batch(es) at once, $WORKERS_PER_JOB workers each (cap ${MAX_CORES}/${TOTAL_CORES} cores)"

# Read the pairs into an array first: a `while read < file` loop would hand
# its stdin to every backgrounded job, which then eats the remaining pairs.
PAIR_LINES=()
if [[ "$RERUN_MODE" == "true" ]]; then
    mapfile -t PAIR_LINES < "$PAIRS_FILE"
    log_message "Rerun units of work: ${#PAIR_LINES[@]} (batch, folder) pair(s)"
fi

start_job() {
    if [ "$JOBS" -le 1 ]; then
        "$@"
    else
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do
            sleep 0.5
        done
        "$@" &
    fi
}

if [[ "$RERUN_MODE" == "true" ]]; then
    for PAIR_LINE in ${PAIR_LINES[@]+"${PAIR_LINES[@]}"}; do
        IFS=$'\t' read -r PAIR_BATCH PAIR_TAB PAIR_FOLDER <<< "$PAIR_LINE"
        [ -n "$PAIR_BATCH" ] || continue
        start_job run_pair "$PAIR_BATCH" "$PAIR_TAB" "$PAIR_FOLDER"
    done
else
    for BATCH in "${BATCHES[@]}"; do
        start_job run_batch "$BATCH"
    done
fi
wait

# Name the folder too: the same batch can fail in one folder and succeed in
# another, and "20260616_mtDNA_389" alone would not say which.
FAILED_BATCHES=()
if [[ "$RERUN_MODE" == "true" ]]; then
    for PAIR_LINE in ${PAIR_LINES[@]+"${PAIR_LINES[@]}"}; do
        IFS=$'\t' read -r PAIR_BATCH PAIR_TAB PAIR_FOLDER <<< "$PAIR_LINE"
        [ -n "$PAIR_BATCH" ] || continue
        STATUS_FILE="${STATUS_DIR}/${PAIR_TAB}__${PAIR_BATCH}"
        if [ ! -f "$STATUS_FILE" ] || [ "$(cat "$STATUS_FILE")" != "0" ]; then
            FAILED_BATCHES+=("${PAIR_BATCH} (${PAIR_TAB})")
        fi
    done
else
    for BATCH in "${BATCHES[@]}"; do
        STATUS_FILE="${STATUS_DIR}/${BATCH}"
        if [ ! -f "$STATUS_FILE" ] || [ "$(cat "$STATUS_FILE")" != "0" ]; then
            FAILED_BATCHES+=("$BATCH")
        fi
    done
fi

if [ ${#FAILED_BATCHES[@]} -gt 0 ]; then
    # Carry on to the merge: the batches that did succeed still have results
    # worth publishing. The failure is reported at the end.
    log_message "WARNING: pipeline failed for: ${FAILED_BATCHES[*]}"
fi

# =============================================================================
# Google Sheets Upload Logic
# =============================================================================

HAS_CREDENTIALS=false
if [[ "$UPLOAD" != "true" ]]; then
    log_message "Google Sheets upload disabled (--no-upload, or rerun-mode default)."
elif [ -f "$CREDENTIALS_FILE" ]; then
    HAS_CREDENTIALS=true
else
    log_message "Google Sheets credentials file not found at $CREDENTIALS_FILE"
    log_message "Will generate comparison files but skip Google Sheets upload."
fi

# Per-batch worksheets, normal mode only. A rerun pulls a handful of samples
# from each of many batches, so one worksheet per batch would bury the sheet in
# near-empty tabs; the merged All_Batches sheet is the whole point there.
# Gated on $HAS_CREDENTIALS as well so --no-upload is honoured.
if $HAS_CREDENTIALS && [[ "$RERUN_MODE" != "true" ]]; then
    for BATCH in "${BATCHES[@]}"; do
        MANUAL_DIR="${RESULTS_DIR}/manual_pipeline/manual/${BATCH}"

        # Guard the directory test: under `set -euo pipefail` a find on a
        # missing path exits 1 and would abort the whole script here.
        if [ -d "$MANUAL_DIR" ]; then
            TXT_COUNT=$(find "$MANUAL_DIR" -name "*.TXT" | wc -l)
        else
            TXT_COUNT=0
        fi
        if [ "$TXT_COUNT" -eq 0 ]; then
            log_message "Skipping $BATCH - no TXT files in manual directory" "$BATCH"
            continue
        fi

        if [ -f "${RESULTS_DIR}/modules/comparison/${BATCH}/${BATCH}.xlsx" ]; then
            WORKSHEET_NAME="Batch_${BATCH}"
            log_message "Uploading ${WORKSHEET_NAME} to the configured sheet" "$BATCH"
            python src/modules/sheets/uploader.py \
                -i "${RESULTS_DIR}/modules/comparison/${BATCH}/${BATCH}.xlsx" \
                -s "$GOOGLE_SHEET_ID" \
                -c "$CREDENTIALS_FILE" \
                -w "$WORKSHEET_NAME" \
                -r
        fi
    done
fi

# Merge the comparison results for the batches processed in this run.
# Only the BATCHES listed above are merged (via -b), so historical batch
# subfolders already present under results/modules/comparison are not pulled
# into the consolidated "all_batches_comparison" output.
BATCHES_STR=$(IFS=,; echo "${BATCHES[*]}")

# A rerun keeps one results tree per source folder, so its comparison files
# live a level deeper than a normal run's.
if [[ "$RERUN_MODE" == "true" ]]; then
    COMPARISON_DIR="${RESULTS_DIR}/by_folder"
else
    COMPARISON_DIR="${RESULTS_DIR}/modules/comparison"
fi
if [ -d "$COMPARISON_DIR" ]; then
    COMPARISON_COUNT=$(find "$COMPARISON_DIR" -path "*/*.xlsx" -not -name "all_batches*" | wc -l)
else
    COMPARISON_COUNT=0
fi
if [ "$COMPARISON_COUNT" -eq 0 ]; then
    log_message "No comparison files were generated. Skipping merge step."
    # Not a bare `exit 0`: a run where every batch failed produced no
    # comparison file either, and reporting that as success hid the failure
    # from both the caller and the web UI.
    if [ ${#FAILED_BATCHES[@]} -gt 0 ]; then
        log_message "ERROR: batch pipeline finished with failures: ${FAILED_BATCHES[*]}"
        exit 1
    fi
    exit 0
fi

# Serialise from here to the end of the upload. Per-batch worksheets are named
# Batch_<id> and never collide, but the merge writes one shared
# all_batches_comparison.xlsx and replaces the single All_Batches worksheet.
# Two runs overlapping here can upload a half-written file. Last writer still
# wins, which is the intended "latest run" behaviour -- this only stops the
# two from interleaving mid-write.
#
# The lock lives on local disk, not beside the results: RESULTS_DIR is usually
# NFS, where flock is unreliable. Key it on the results path so runs writing to
# different directories do not block each other.
MERGE_LOCK_DIR="${TMPDIR:-/tmp}/mtdna-locks"
mkdir -p "$MERGE_LOCK_DIR"
MERGE_LOCK="${MERGE_LOCK_DIR}/merge-$(printf '%s' "$RESULTS_DIR" | md5sum | cut -c1-16).lock"

exec {MERGE_LOCK_FD}>"$MERGE_LOCK"
if ! flock -w 600 "$MERGE_LOCK_FD"; then
    log_message "ERROR: timed out waiting for the merge lock at $MERGE_LOCK"
    exit 1
fi

# Everything merged lands in All_Batches. A rerun is additionally split by the
# folder its ZIPs came from: a rerun is normally pointed at the whole
# .../Sequencher_temp/rerun tree, where each session has its own subfolder, so
# naming the tab after the top of the tree would pile every session into one.
# Each folder gets a worksheet holding only the batches whose ZIPs are in it,
# and a batch found in two folders appears in both.
# Merged output stays at the top level, outside the per-folder trees the
# merger scans.
MERGED_DIR="${RESULTS_DIR}/modules/comparison"
mkdir -p "$MERGED_DIR"

MERGED_GROUPS=()

if [[ "$RERUN_MODE" == "true" ]]; then
    # A rerun gets one merged table per source folder and no All_Batches. The
    # folders are separate rerun sessions -- often the same batch and the same
    # samples analysed twice -- so a combined table would carry two rows per
    # sample with nothing to tell them apart.
    declare -A GROUP_BATCHES=()
    for PAIR_LINE in ${PAIR_LINES[@]+"${PAIR_LINES[@]}"}; do
        IFS=$'\t' read -r PAIR_BATCH PAIR_TAB PAIR_FOLDER <<< "$PAIR_LINE"
        [ -n "$PAIR_BATCH" ] || continue
        GROUP_BATCHES["$PAIR_TAB"]+="${PAIR_BATCH},"
    done

    for GROUP in $(printf '%s\n' "${!GROUP_BATCHES[@]}" | sort); do
        log_message "Merging comparison results for folder ${GROUP}: ${GROUP_BATCHES[$GROUP]%,}"
        # `if ! python` rather than a bare call: a folder whose batches all
        # failed has no comparison file to merge, and under `set -e` that would
        # abort the run and take the remaining folders with it.
        if ! python src/modules/sheets/merger.py \
            -i "${RESULTS_DIR}/by_folder/${GROUP}/modules/comparison" \
            -o "${MERGED_DIR}/comparison_${GROUP}" \
            -b "${GROUP_BATCHES[$GROUP]%,}"; then
            log_message "WARNING: no comparison result for folder ${GROUP}, skipping its worksheet"
            continue
        fi
        MERGED_GROUPS+=("$GROUP")
    done
else
    log_message "Merging comparison results for batches: ${BATCHES_STR}"
    python src/modules/sheets/merger.py \
        -i "$COMPARISON_DIR" \
        -o "${MERGED_DIR}/all_batches_comparison" \
        -b "${BATCHES_STR}"
fi

# Upload. HAS_CREDENTIALS is already false when $UPLOAD is not "true", so a run
# that did not ask to upload cannot touch any worksheet.
if $HAS_CREDENTIALS; then
    if [[ "$RERUN_MODE" == "true" ]]; then
        for GROUP in ${MERGED_GROUPS[@]+"${MERGED_GROUPS[@]}"}; do
            log_message "Uploading ${GROUP} to the configured sheet"
            python src/modules/sheets/uploader.py \
                -i "${MERGED_DIR}/comparison_${GROUP}.xlsx" \
                -s "$GOOGLE_SHEET_ID" \
                -c "$CREDENTIALS_FILE" \
                -w "$GROUP" \
                -r
        done
    elif [ -f "${MERGED_DIR}/all_batches_comparison.xlsx" ]; then
        log_message "Uploading All_Batches to the configured sheet"
        python src/modules/sheets/uploader.py \
            -i "${MERGED_DIR}/all_batches_comparison.xlsx" \
            -s "$GOOGLE_SHEET_ID" \
            -c "$CREDENTIALS_FILE" \
            -w "All_Batches" \
            -r
    fi
fi

flock -u "$MERGE_LOCK_FD"
exec {MERGE_LOCK_FD}>&-

if [ ${#FAILED_BATCHES[@]} -gt 0 ]; then
    log_message "ERROR: batch pipeline finished with failures: ${FAILED_BATCHES[*]}"
    exit 1
fi

log_message "Batch pipeline completed successfully!"
