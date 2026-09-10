#!/usr/bin/env bash
set -euo pipefail

source .env

# # Merge JSON files for statistics of Sanger
# python src/modules/NGS/merge_json_files.py \
#     results/20250710_mtDNA_85_validation_NGS.json \
#     resuyeah \
#     -o results/merged_statistics.json

# # Merge JSON files for statistics of Sanger
# python src/modules/NGS/merge_json_files.py \
#     results/statistic_fullbatch.json \
#     results/merged_statistics.json \
#     -o results/merged_statistics.json

####################################################################################################

# # Merge variants from FIS
# python src/modules/NGS/FIS_merge_variants.py \
#     -o results/FIS_merged_variants_20250820.tsv \
#     -e "data/20250820/*MT.xlsx"

# Remove generated outputs from the previous run while preserving unrelated files.
BATCH_ID="20260820"
OUTPUT_DIR="results/modules/NGS/$BATCH_ID"
COMPARISON_DIR="$OUTPUT_DIR/comparison"
CORRECTION_FILE="results/modules/NGS/correct_mapping_sanger.tsv"

test -f "$CORRECTION_FILE" || {
    echo "Missing authoritative FIS-to-Sanger correction file: $CORRECTION_FILE" >&2
    exit 1
}

rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Stage 1: create canonical FIS Sample JSON and established NGS artifacts.
python -m src.tools.fis.pipeline \
    --r1-file temp/20260817.result_3040_20260817.R1.DataProduction.xlsx \
    --r2-file temp/20260817.result_3040_20260817.R2.MT.xlsx \
    --output-dir "$OUTPUT_DIR" \
    --batch-id "$BATCH_ID" \
    --write-ngs-artifacts

# python src/modules/NGS/compare_fis_sanger.py \
#     -f results/modules/NGS/20260820/FIS_transformed.tsv \
#     -s /mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json \
#     -c data/NGS/mapping.csv \
#     -t data/NGS/STR_mtDNA_with_ids_family_with_variants.tsv \
#     -o results/fis_sanger_comparison.tsv

echo "####################################################################################################"

# Stage 2: map and normalize FIS samples for shared comparison.
python -m src.modules.NGS.prepare_fis_comparison \
    -f "$OUTPUT_DIR/json/statistic_fullbatch.json" \
    -s /mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json \
    -t results/modules/NGS/mapping.tsv \
    -u "$CORRECTION_FILE" \
    -o "$COMPARISON_DIR/FIS.json" \
    -m "$COMPARISON_DIR/sample_mapping.tsv" \
    --batch-id "$BATCH_ID"

# Stage 3: select and normalize only the mapped Sanger samples.
python -m src.modules.NGS.prepare_sanger_comparison \
    -s /mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json \
    -m "$COMPARISON_DIR/sample_mapping.tsv" \
    -o "$COMPARISON_DIR/Sanger.json" \
    --batch-id "$BATCH_ID"

# Stage 4: run the centralized comparison logic.
python -m src.core.comparison \
    -a "$COMPARISON_DIR/FIS.json" \
    -b "$COMPARISON_DIR/Sanger.json" \
    -o "$COMPARISON_DIR" \
    --batch-id "$BATCH_ID"

# Stage 5: add FIS metadata and unmatched samples to the final report.
python -m src.modules.NGS.enrich_fis_sanger_comparison \
    -c "$COMPARISON_DIR/$BATCH_ID.json" \
    -f "$OUTPUT_DIR/json/statistic_fullbatch.json" \
    -t "$OUTPUT_DIR/FIS_transformed.tsv" \
    -s "$COMPARISON_DIR/Sanger.json" \
    -m "$COMPARISON_DIR/sample_mapping.tsv" \
    -o "$COMPARISON_DIR/FIS_sanger_comparison_direct.tsv"

# Stage 6: split the enriched report into the three defined sample sets.
python -m src.modules.NGS.split_fis_sanger_comparison_sets \
    -i "$COMPARISON_DIR/FIS_sanger_comparison_direct.tsv" \
    -s results/modules/NGS/sets.tsv \
    -o "$COMPARISON_DIR"

# Stage 7: write the batch-local called-heteroplasmy report.
python -m src.modules.NGS.generate_heteroplasmy_report \
    -f "$OUTPUT_DIR/FIS.tsv" \
    -c "$COMPARISON_DIR/${BATCH_ID}_variant_level.tsv" \
    -m "$COMPARISON_DIR/sample_mapping.tsv" \
    -s results/modules/NGS/sets.tsv \
    -o "$OUTPUT_DIR/Heteroplasmy.tsv"

####################################################################################################

# # Step 1: Consensus FASTA merging script
# python src/modules/NGS/FIS_merge_fasta.py \
#     -o results/FIS_FASTA \
#     -e "/mnt/nas/earth/Validation_NGS_384/Result_384_FIS/V350346430_L01_10-10_13428/FinalResult/Sub_web/*/2-MT/*-HyperVar-MT.xlsx"

# # Step 2: Convert FASTA to JSON (FASTA-input BLASTN pipeline; writes json/ + regions/ + preprocess/)
# python -m src.tools.blastn.fasta_pipeline \
#     -i results/20250916.R2.MT.FASTA \
#     -o results/20250916.R2.MT.JSON \
#     -r "$REFERENCE" \
#     --batch-id 20250916.R2.MT

# python src/modules/NGS/convert_variants_to_tsv.py \
#     results/20250916.R2.MT.JSON/json/statistic_fullbatch.json \
#     results/FIS_merged_variants_transformed_FASTA.tsv

# # Basic usage
# python -m src.modules.NGS.compare_fis_sanger_direct \
#     -f results/FIS_merged_variants_transformed_FASTA.tsv \
#     -s /mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json \
#     -t results/modules/NGS/mapping.tsv \
#     -o results/FIS_sanger_comparison_direct_FASTA.tsv

# python src/modules/NGS/compare_fis_fasta_direct.py \
#     -f results/FIS_merged_variants_transformed.tsv \
#     -s results/FIS_FASTA_JSON/statistic_fullbatch.json \
#     -o results/FIS_FASTA_comparison.tsv

####################################################################################################

# # Merge comparison results from FIS-Sanger and FIS-FASTA comparisons
# python src/modules/NGS/merge_comparison.py \
#     -s results/fis_sanger_comparison.tsv \
#     -f results/FIS_FASTA_comparison.tsv \
#     -o results/merged_comparison_results.tsv
