BATCHES=(
    "20260610_mtDNA_427"
    "20260615_mtDNA_438"
    "20260615_mtDNA_439"
    "20260615_mtDNA_440"
    "20260615_mtDNA_441"
    "20260615_mtDNA_442"
    "20260615_mtDNA_443"
    "20260615_mtDNA_444"
    "20260615_mtDNA_446"
)

IN_DIR="temp/9_batch_260623"

for BATCH in "${BATCHES[@]}"; do
    # Resolve per-batch truthset file from IN_DIR
    TRUTHSET=$(ls "${IN_DIR}"/${BATCH}_*_Truthset.xlsx 2>/dev/null | head -1)

    bash scripts/tools/mutation_surveyor.sh \
        ${BATCH} \
        -s 1,2,3,4,5,6,7,8 \
        "${IN_DIR}" \
        /mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/mutation_surveyor/${BATCH} \
        /mnt/nas/bca/mtDNA/science/data/raw/${BATCH} \
        ${TRUTHSET:+--truthset "${TRUTHSET}"}

    # rm -rf /mnt/nas/bca/mtDNA/science/mutation_surveyvor/results/${BATCH}/*.xlsx

    # cp /mnt/nas/bca/mtDNA/science/mutation_surveyvor/results/${BATCH}/*.xlsx /mnt/nas/bca/mtDNA/science/mutation_surveyvor/excel

    # mv /mnt/nas/bca/mtDNA/science/mutation_surveyvor/excel/${BATCH}_df_custom_report.xlsx /mnt/nas/bca/mtDNA/science/mutation_surveyvor/excel/${BATCH}_Analysis_custom_settings_260609.xlsx

    # cp /mnt/nas/bca/mtDNA/science/mutation_surveyvor/results/results_120626/${BATCH}/*.xlsx /mnt/nas/bca/mtDNA/science/mutation_surveyvor/excel
done
