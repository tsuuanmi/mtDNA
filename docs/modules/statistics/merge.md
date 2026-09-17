# mtDNA data merger

Consolidate per-batch metadata and results into the merged exports that other
modules (population analysis, variant frequencies, TNLS/NGS runners) consume:

```bash
uv run python -m src.modules.statistics.merge                 # production NAS defaults
uv run python -m src.modules.statistics.merge BATCH_A BATCH_B  # replace default exclusions
uv run python -m src.modules.statistics.merge \
  --data-dir /mnt/nas/bca/mtDNA/science/data \
  --results-dir /home/bca/mtDNA/science/results \
  --output-dir /home/bca/mtDNA/science/results/merged
uv run python -m src.modules.statistics.merge --force-refresh  # bypass all caches
```

To use `.env` values explicitly:

```bash
source .env
uv run python -m src.modules.statistics.merge --data-dir "${DATA_DIR}" --results-dir "${RESULTS_DIR}"
```

## Configuration and CLI

Inputs and outputs are CLI arguments following the population module's
convention: the defaults point at the production NAS paths and every location
can be overridden per run.

| Path | Role |
|---|---|
| `<data-dir>/metadata/` | Per-batch Excel workbooks and TSVs (sample_id, barcode) |
| `<results-dir>/manual_pipeline/regenerate/<batch>/` | Regenerate results; `statistic_fullbatch.json` keys are sample IDs |
| `<results-dir>/fasta/<batch>/`, `<results-dir>/archive/<batch>/` | FASTA sample sets |
| `<output-dir>/` | All outputs below |

| Argument | Meaning |
|---|---|
| `excluded_batches...` | Batch names to exclude from every scan and merge. Positional; **replaces** the default list (`20250429_mtDNA_100`) entirely when given |
| `--data-dir` | Root data directory containing `metadata/` (default: `/mnt/nas/bca/mtDNA/science/data`) |
| `--results-dir` | Results root containing `manual_pipeline/regenerate/`, `fasta/` and `archive/` (default: `/mnt/nas/bca/mtDNA/science/results`) |
| `--output-dir` | Directory for merged outputs and the scan cache (default: `<results-dir>/merged`) |
| `--force-refresh` | Ignore the directory-scan cache, the metadata cache and the merged-statistics manifest for this run; all caches are rebuilt |

## Inputs and conversion

- Metadata Excel workbooks (`*.xlsx`, first `LID` and `BARCODE` columns) are
  converted to headerless per-batch TSVs with numeric barcodes zero-padded to
  the longest observed length. Conversion happens **only** when the batch TSV
  is missing or unreadable; existing readable TSVs are never rewritten.
  Readability is checked concurrently (32 workers) because the metadata
  directory lives on a high-latency NAS share.
- Metadata TSVs map sample_id → barcode from the first two tab-separated
  columns of every non-empty line. Unchanged TSVs are served from the cache;
  fresh files are read concurrently.
- Regenerate sample IDs are the keys of each batch's `statistic_fullbatch.json`;
  batches without that JSON fall back to listing their sample directories.
- FASTA sample IDs are file stems under the active and archive batch roots;
  archive results add to the same per-batch view.

## Caching and performance

`<output-dir>/scan_cache.json` (schema version 3) caches every category —
`regenerate`, `fasta`, `archive_fasta` and `metadata` — per source file/directory
keyed by path, type, `mtime_ns`, `ctime_ns` and size. Repeat runs only stat
unchanged sources; a changed, added or removed source is re-read, and an
incompatible or corrupt cache is ignored and rebuilt once. The merged-statistics
JSON additionally carries a manifest (below); it is reused verbatim while all
source batch JSONs are unchanged and no source is newer than the output.
Sample-level duplicate IDs keep the **last batch** entry in `merged_statistics.json`
(batches merge in batch-name order); the merged TSV may contain the same sample
ID under several batches as separate rows.

## Outputs

| File | Contents |
|---|---|
| `merged_metadata.tsv` | `batch`, `sample_id`, `barcode` for every scanned sample (regenerate ∪ FASTA) that exists in metadata |
| `merged_statistics.json` | Sample-ID-keyed merge of all per-batch `statistic_fullbatch.json` records, each annotated with `batch` |
| `merged_statistics.manifest.json` | `created_at`, sorted `excluded_batches`, per-source signature (`path`, `mtime_ns`, `ctime_ns`, `size`), `batch_count` |
| `batch_statistics.tsv` | Per-batch Successful/Metadata/Regenerate/Fasta counts, status and a Total row |
| `sample_discrepancies.tsv` | One row per unique sample ID across the merged TSV and JSON: `discrepancies` flags (`duplicate_rows`, `in_tsv_not_json`, `in_json_not_tsv`, `matched`), merged-TSV row count and batches, JSON batch |
| `scan_cache.json` | Shared scan and metadata cache described above |

`Successful` counts samples present in both the scan results and metadata.
`Status` is `No Meta` when a batch has no metadata rows, `No FASTA` without
FASTA samples, `No Pass` when nothing scanned is in metadata, `OK` when
metadata and successful counts match, else `Mismatch`. The run ends with a
discrepancy report comparing unique TSV and JSON sample IDs in both directions;
duplicate TSV rows are counted separately. `sample_discrepancies.tsv` tracks the
same information as data: one row per unique sample ID with every applicable
flag, the TSV row count and batches, and the JSON batch. Flags can combine — a
sample duplicated across batches and absent from JSON carries
`duplicate_rows;in_tsv_not_json` — and each flag's row count equals the matching
console-log count.

Failures are explicit: unreadable regenerate JSON falls back to a directory
listing per batch with a warning; a corrupt batch statistics JSON is logged,
its batch skipped, and the merge still written; unexpected `OSError` exits 1.

## Python API

Path-dependent functions take a frozen `MergePaths` contract; build it from the
same roots the CLI resolves:

```python
from pathlib import Path

from src.modules.statistics.merge import (
    DEFAULT_DATA_DIR, DEFAULT_RESULTS_DIR, DirectoryScanner, MergePaths,
    analyze_sample_discrepancies, batch_statistics, create_merged_tsv,
    get_successful_samples, merge_statistics_json, process_excel_files,
    read_metadata_files,
)

excluded = ["20250429_mtDNA_100"]
paths = MergePaths(
    data_dir=DEFAULT_DATA_DIR,
    results_dir=DEFAULT_RESULTS_DIR,
    output_dir=Path("results/merged-dry-run"),  # CLI default: <results_dir>/merged
)

process_excel_files(excluded, paths)
scanner = DirectoryScanner(paths, excluded)
scan_result = scanner.scan_all_directories()          # or force_refresh=True
successful_samples = get_successful_samples(scan_result)
metadata = read_metadata_files(excluded, paths)       # or force_refresh=True
samples_per_batch, total = create_merged_tsv(successful_samples, metadata, paths)
batch_statistics(samples_per_batch, successful_samples, metadata, scan_result, paths)
merged_statistics, n_batches = merge_statistics_json(excluded, paths)
analyze_sample_discrepancies(merged_statistics, total, paths)
```

Conversions are also available directly:
`convert_excel_to_tsv(excel_file, tsv_path)`. Importing the package requires
no environment variables; all locations come from `MergePaths` or the CLI.
