rm -rf logs/test.log

# bash scripts/pipeline.sh \
#     --pipeline "tracy" \
#     --steps "2,3,4,5" \
#     --data-dir "temp/test" \
#     --manual-dir "temp/test/sequencher" \
#     --sample-list "temp/test/samples.txt" \
#     --output-dir "results" \
#     "test"

bash scripts/pipeline.sh \
    --pipeline "tracy" \
    --steps "2,3,4" \
    --output-dir "results" \
    "MS_110926_008"
