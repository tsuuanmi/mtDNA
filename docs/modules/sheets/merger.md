# Excel File Merger

> Merge multiple Sequencher comparison Excel files into a single consolidated file.

## Purpose

`src/modules/sheets/merger.py` provides a CLI tool for merging per-batch comparison Excel files (recursively found under a directory) into a single consolidated Excel and TSV file. It normalizes variant separators and sorts by batch and sample ID.

## Public API

### `merge_excel_files(input_dir, output_file, batches=None)`

```python
merge_excel_files(input_dir: str | Path, output_file: str | Path,
                   batches: list[str] | None = None) → None
```

Merge per-batch `*.xlsx` files found recursively under `input_dir` (excluding `all_batches*` merged outputs and temporary `~$` files). Writes both `.xlsx` and `.tsv` to the output path (with appropriate suffix).

When `batches` is a list of batch IDs, only the per-batch subfolders whose name matches one of the listed IDs are merged. Pass `None` (the default) to merge every per-batch file found under `input_dir`. `scripts/batch_pipeline.sh` passes the run's `BATCHES` array so the consolidated output reflects only the batches processed in that run instead of every historical batch subfolder.

**Steps:**
1. Load all matching Excel files
2. Normalize tool-parameterized or `(A)`/`(B)` columns to legacy `Pipeline`/`Sequencher` names
3. Concatenate into a single DataFrame
4. Normalize sample IDs to strings
5. Sort by batch and sample ID
6. Normalize variant separators (space-delimited)
7. Write merged `.xlsx` and `.tsv` files

### Module-Level Functions

| Function | Description |
|----------|-------------|
| `arg_parser()` | Parse CLI arguments for `--input_dir` (`-i`), `--output_file` (`-o`), and optional `--batches` (`-b`, comma-separated batch IDs) |
| `_normalize_columns(df)` | Map tool-parameterized or `(A)`/`(B)` column names (e.g. `Variants (tracy)`) to legacy `Pipeline`/`Sequencher` labels so batches align |
| `_load(input_dir, batches=None)` | Recursively load per-batch `*.xlsx` files, normalizing columns before concat (excludes merged/temp files); when `batches` is set, only subfolders matching a listed batch ID are included |
| `_normalize_separator(df, separator)` | Normalize variant column separators to space or comma |
| `_write_excel(df, path)` | Write DataFrame to Excel with auto-sized columns |
| `_write_tsv(df, path)` | Write DataFrame to TSV |

### CLI Entry Point

```bash
python -m src.modules.sheets.merger \
    -i <input_dir> \
    -o <output_file> \
    [-b <batch_id,batch_id,...>]
```

`-b/--batches` (optional) restricts the merge to the listed comma-separated batch IDs (matched against the per-batch subfolder name under `input_dir`). Omit it to merge every per-batch file found. `scripts/batch_pipeline.sh` builds this list from its `BATCHES` array so the consolidated output covers only the batches processed in that run.

## Key Concepts

### Variant Columns

The merger normalizes separators in these legacy summary columns:
- `Variants (Pipeline)`
- `Variants (Sequencher)`
- `Variants Unique (Pipeline)`
- `Variants Unique (Sequencher)`

### Column Normalization

Comparison outputs may name columns by the source tool (e.g. `Variants (tracy)` / `Variants (sequencher)`) or by the generic `A`/`B` labels. Different batches may use different automate tools (tracy vs blastn), so raw column names vary. The merger maps the first tool in column order to `Pipeline` (the automate/pipeline side) and the second to `Sequencher` (the manual side) before concatenating, so all batches align while preserving the legacy summary headers.

### Output Format

Both `.xlsx` and `.tsv` files are written with the same merged data. The Excel file includes auto-sized columns (capped at 50 characters width).

## Cross-References

- [core/comparison.md](../../core/comparison.md) — Comparison results that feed into the merger
- [pipeline.md](../../pipeline.md) — Pipeline orchestration
