rm -rf logs/test.log

# bash scripts/modules/statistics/test.sh

bash scripts/pipeline.sh \
    --pipeline "tracy" \
    --steps "2,3,4,5" \
    --data-dir "temp/test" \
    --manual-dir "temp/MS_130826_006_V565_X" \
    --sample-list "temp/test/samples.txt" \
    --output-dir "results" \
    "test"
