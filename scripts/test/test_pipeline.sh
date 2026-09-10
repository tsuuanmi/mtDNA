rm -rf logs/test.log

# bash scripts/modules/statistics/test.sh

bash scripts/pipeline.sh \
    --pipeline "tracy" \
    --steps "2,3,4,5" \
    --data-dir "/home/tan/workspaces/mtdna_raw/temp/test" \
    --manual-dir "/home/tan/workspaces/mtdna_raw/temp/MS_130826_006_V565_X" \
    --sample-list "/home/tan/workspaces/mtdna_raw/temp/test/samples.txt" \
    --output-dir "results" \
    "test"
