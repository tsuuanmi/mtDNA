# LAB

LAB utilities generate HTML and Excel mtDNA tracking reports.

Run all available FASTA batches from the repository root:

```bash
bash scripts/modules/LAB/run.sh
```

The runner reads the local, ignored files in `data/modules/LAB/`, requires each batch metadata workbook in `${DATA_DIR}/metadata/`, maps date-based and MS batch IDs to tracking data, and writes completed reports to `temp/modules/LAB/`.
