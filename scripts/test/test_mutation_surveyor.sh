source .env

BATCHES=(
    "20260615_mtDNA_444"
)

IN_DIR="temp/9_batch_260623"

for BATCH in "${BATCHES[@]}"; do
    # Resolve per-batch truthset file from IN_DIR
    TRUTHSET=$(ls "${IN_DIR}"/${BATCH}_*_Truthset.xlsx 2>/dev/null | head -1)

    # bash scripts/tools/mutation_surveyor.sh \
    #     ${BATCH} \
    #     -s 1,2,3,4,5,6,7,8,9 \
    #     "${IN_DIR}" \
    #     /mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/mutation_surveyor/${BATCH} \
    #     /mnt/nas/bca/mtDNA/science/data/raw/${BATCH} \
    #     ${TRUTHSET:+--truthset "${TRUTHSET}"}

    # bash scripts/tools/tracy.sh \
    #     -s "1" \
    #     --data-dir "/mnt/nas/bca/mtDNA/science/data/raw/${BATCH}" \
    #     --sample-list "/mnt/nas/bca/mtDNA/science/data/raw/${BATCH}.txt" \
    #     --output-dir "results" \
    #     "${BATCH}"

    # bash scripts/tools/blastn.sh \
    #     -s "1" \
    #     --data-dir "/mnt/nas/bca/mtDNA/science/data/raw/${BATCH}" \
    #     --sample-list "/mnt/nas/bca/mtDNA/science/data/raw/${BATCH}.txt" \
    #     --output-dir "results/tools/blastn/${BATCH}" \
    #     "${BATCH}"

    # python -m src.tools.sequencher.pipeline \
    #     --input-dir "/mnt/nas/bca/mtDNA/science/results/manual_pipeline/manual/${BATCH}" \
    #     --output-dir "results/tools/sequencher" \
    #     --ref-path "$REFERENCE" \
    #     --batch-id "${BATCH}"

    # python -m src.core.comparison \
    #     --batch-a "results/tools/sequencher/${BATCH}/statistic_fullbatch.json" \
    #     --batch-b "/mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/json_temp/${BATCH}/unified/statistic_fullbatch.json" \
    #     --output-dir "/mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/comparison/${BATCH}" \
    #     --batch-id "${BATCH}"

    python -m src.core.comparison \
        --batch-a "/mnt/nas/bca/mtDNA/science/results/automate_pipeline/regenerate/${BATCH}/statistic_fullbatch.json" \
        --batch-b "/mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/json_temp/${BATCH}/unified/statistic_fullbatch.json" \
        --output-dir "/mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/comparison/${BATCH}" \
        --batch-id "${BATCH}"
done
