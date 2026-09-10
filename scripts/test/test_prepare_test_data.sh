#!/bin/bash
set -euo pipefail

# =============================================================================
# Prepare 72-sample test data folder for pytest compatibility tests
# Extracts from the zip and reorganizes into the pipeline's expected structure.
# =============================================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)

# Configuration
ZIP_FILE="${ZIP_FILE:-/mnt/nas/bca/mtDNA/science/temp/test_set_72_samples.zip}"
TEST_DIR="${TEST_DIR:-/mnt/nas/bca/mtDNA/science/tests}"
BATCH_ID="test"
TMP_EXTRACT="/tmp/test_set_72_extract"

echo "=== Preparing test data ==="
echo "ZIP:  $ZIP_FILE"
echo "DIR:  $TEST_DIR"
echo "BATCH: $BATCH_ID"

# Validate zip exists
if [ ! -f "$ZIP_FILE" ]; then
    echo "ERROR: Zip file not found: $ZIP_FILE"
    exit 1
fi

# Step 1: Extract zip to temp
echo "[1/5] Extracting zip..."
rm -rf "$TMP_EXTRACT"
mkdir -p "$TMP_EXTRACT"
unzip -q -o "$ZIP_FILE" -d "$TMP_EXTRACT"

SOURCE="$TMP_EXTRACT/test_set_72_samples"

if [ ! -d "$SOURCE" ]; then
    echo "ERROR: Expected directory not found after extraction: $SOURCE"
    ls "$TMP_EXTRACT"
    exit 1
fi

# Step 2: Create target directory structure
echo "[2/5] Creating directory structure..."
mkdir -p "$TEST_DIR/data/raw/$BATCH_ID"
mkdir -p "$TEST_DIR/results/manual_pipeline/manual/$BATCH_ID"
mkdir -p "$TEST_DIR/results/manual_pipeline/regenerate/$BATCH_ID"
mkdir -p "$TEST_DIR/results/fasta/$BATCH_ID"

# Step 3: Copy and reorganize files
echo "[3/5] Copying files..."

# Sample list → data/raw/test.txt
cp "$SOURCE/List_72.txt" "$TEST_DIR/data/raw/$BATCH_ID.txt"
echo "  Sample list: $(wc -l < "$TEST_DIR/data/raw/$BATCH_ID.txt") samples"

# AB1 files (flatten subdirectories) → data/raw/test/
find "$SOURCE/raw_ab1" -name '*.ab1' -exec cp -f {} "$TEST_DIR/data/raw/$BATCH_ID/" \;
AB1_COUNT=$(find "$TEST_DIR/data/raw/$BATCH_ID" -name '*.ab1' | wc -l)
echo "  AB1 files: $AB1_COUNT"

# Sequencher TXT files → results/manual_pipeline/manual/test/
find "$SOURCE/sequencher_txt" -name '*.TXT' -exec cp -f {} "$TEST_DIR/results/manual_pipeline/manual/$BATCH_ID/" \;
TXT_COUNT=$(find "$TEST_DIR/results/manual_pipeline/manual/$BATCH_ID" -name '*.TXT' | wc -l)
echo "  TXT files: $TXT_COUNT"

# FASTA files → results/fasta/test/
find "$SOURCE/fasta" -name '*.fasta' -exec cp -f {} "$TEST_DIR/results/fasta/$BATCH_ID/" \;
FASTA_COUNT=$(find "$TEST_DIR/results/fasta/$BATCH_ID" -name '*.fasta' | wc -l)
echo "  FASTA files: $FASTA_COUNT"

# JSON files → results/manual_pipeline/regenerate/test/{sample_id}/{sample_id}.json
JSON_COUNT=0
for json_file in "$SOURCE/json"/*.json; do
    sample_id=$(basename "$json_file" .json)
    mkdir -p "$TEST_DIR/results/manual_pipeline/regenerate/$BATCH_ID/$sample_id"
    cp -f "$json_file" "$TEST_DIR/results/manual_pipeline/regenerate/$BATCH_ID/$sample_id/$sample_id.json"
    JSON_COUNT=$((JSON_COUNT + 1))
done
echo "  JSON files: $JSON_COUNT"

# Step 4: Summary
echo "[4/5] Verifying structure..."
echo ""
echo "$TEST_DIR/"
find "$TEST_DIR" -type d | sort | sed 's|'"$TEST_DIR"'||' | while read -r dir; do
    count=$(find "$TEST_DIR$dir" -maxdepth 1 -type f | wc -l)
    echo "  $dir/ ($count files)"
done

echo ""
echo "[5/5] Done! Test data prepared at: $TEST_DIR"
echo ""
echo "Expected sample count: 72"
echo "Samples in list:      $(wc -l < "$TEST_DIR/data/raw/$BATCH_ID.txt")"
echo "AB1 files:            $AB1_COUNT"
echo "TXT files:            $TXT_COUNT"
echo "FASTA files:          $FASTA_COUNT"
echo "JSON files:           $JSON_COUNT"

# Cleanup
rm -rf "$TMP_EXTRACT"
echo ""
echo "Cleaned up temp extraction directory."
