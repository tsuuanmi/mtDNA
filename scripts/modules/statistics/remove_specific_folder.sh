#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR="/mnt/nas/bca/mtDNA/science/results"

# Folder names / batches to find and remove under RESULTS_DIR
BATCHES=(
    "20260626_mtDNA_484SSS"
    "20260529_mtDNA_377s"
    "202060529_mtDNA_372"
    "20200610_mtDNA_430"
)

declare -A batch_names=()
declare -A found_batches=()
matched_paths=()

for batch in "${BATCHES[@]}"; do
    batch_names["$batch"]=1
    found_batches["$batch"]=false
done

while IFS= read -r -d '' path; do
    folder_name="${path##*/}"

    if [[ -n "${batch_names[$folder_name]:-}" ]]; then
        matched_paths+=("$path")
        found_batches["$folder_name"]=true
    fi
done < <(find "$RESULTS_DIR" -type d -print0)

for path in "${matched_paths[@]}"; do
    echo "Removing: $path"
    rm -rf -- "$path"
done

for batch in "${BATCHES[@]}"; do
    if [[ "${found_batches[$batch]}" == false ]]; then
        echo "Not found, skipping: $batch"
    fi
done

echo "Done."
