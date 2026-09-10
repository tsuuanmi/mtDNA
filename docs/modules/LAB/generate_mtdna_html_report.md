# LAB HTML Report Generation

> Generate HTML and Excel reports for mtDNA sample tracking and result review.

## Purpose

`src/modules/LAB/generate_mtdna_html_report.py` generates HTML and Excel reports from tracking TSV and metadata Excel files for lab review of mtDNA results.

## Public API

| Function | Description |
|----------|-------------|
| `format_date(date_str)` | Format a date string for report display |
| `load_tracking_tsv(path)` | Load tracking data from TSV file |
| `load_metadata_excel(path)` | Load sample metadata from Excel file |
| `load_rerun_excel(path)` | Load rerun data from Excel file |
| `extract_date_and_batch(name)` | Extract date and batch ID from a name string |
| `generate_html_report(tracking_data, metadata_data, rerun_data, output_dir, batch_id)` | Generate HTML report with tracking and metadata |
| `generate_excel_report(tracking_data, metadata_data, rerun_data, output_dir, batch_id)` | Generate Excel report with tracking and metadata |

### CLI Entry Point

```bash
python -m src.modules.LAB.generate_mtdna_html_report [options]
```

## Configuration

- Input files: tracking TSV, metadata Excel, rerun Excel (optional)
- Output: HTML report and Excel report in specified output directory

## Cross-References

- [pipeline.md](../../pipeline.md) — Pipeline orchestration
