source .env

BATCHES=(
    "20260615_mtDNA_444"
)

for BATCH in "${BATCHES[@]}"; do
    # python -m src.tools.sequencher.pipeline \
    #     --input-dir "/mnt/nas/bca/mtDNA/science/results/manual_pipeline/manual/${BATCH}" \
    #     --output-dir "results/tools/sequencher/${BATCH}" \
    #     --ref-path "$REFERENCE" \
    #     --batch-id "${BATCH}"

    python -m src.tools.sequencher.pipeline \
        --input-dir "temp/20260817" \
        --output-dir "results/tools/sequencher/20260817" \
        --ref-path "$REFERENCE" \
        --batch-id "20260817"

    python -m src.generate_fasta \
        -i "results/tools/sequencher/20260817/json" \
        -o "results/fasta/20260817" \
        -r "ref/rCRS.fasta" \
        --results-dir "results" \

done
