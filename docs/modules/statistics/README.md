# Statistics Module

> Query and filter mtDNA variant data with multi-filter support.

## Purpose

`src/modules/statistics/` provides utilities for querying, aggregating, and summarizing mtDNA variant data. `query_multi_filters.py` supports multi-filter queries with region-based windowing and variant formatting. The `merge/` package consolidates data and per-batch results into merged outputs (`python -m src.modules.statistics.merge`). It also contains the blinded TNLS test-data helpers: `blind_copy.py`, `add_variants_to_samples.py`, and `add_variants_to_metadata.py`.

Run the merger with `uv run python -m src.modules.statistics.merge`.

## Merge package

`src/modules/statistics/merge/` merges per-batch metadata and results into `<results>/merged/` by default: `merged_metadata.tsv`, `merged_statistics.json` (+ freshness manifest), `batch_statistics.tsv`, and `sample_discrepancies.tsv` (per-sample TSV-vs-JSON discrepancy tracking). Input and output locations are CLI arguments like the population module (`--data-dir`, `--results-dir`, `--output-dir`; defaults are the production NAS paths, and `.env` values can be passed explicitly). Directory scans and metadata TSV reads are cached in `scan_cache.json` keyed by source file metadata, so repeat runs only stat unchanged sources; `--force-refresh` bypasses both caches. NAS metadata TSVs are read concurrently (32 workers), and Excel workbooks are re-converted only when their batch TSV is missing or unreadable. See [merge.md](merge.md) for the full CLI, inputs, outputs, caching behavior and API.

Prepare the blinded TNLS test data with `bash scripts/modules/statistics/test.sh`; it recreates `/home/tan/workspaces/mtdna_raw/temp/test`. Maintain its local, ignored blind-ID mapping at `data/modules/statistics/sample_mapping.tsv`.

## Population frequency and haplotype diversity

Run `uv run python -m src.modules.statistics.population` against the current
upstream-approved `merged_statistics.json` export (or set `--input`). Without
`--output-dir`, reports are written to
`<results>/modules/statistics/population/current-cohort-<timestamp>`.
This module uses callable denominators and complete HV1+HV2+HV3 profiles; it does not
use the older coverage-independent matching or whole-cohort variant denominator.
See [population.md](population.md) for input requirements, output schemas, methods and API.

## Public API

| Function | Description |
|----------|-------------|
| `variants_in_windows(variants, regions)` | Group variants by genomic region windows |
| `format_analyzed_range(intervals)` | Format analyzed intervals for display |
| `format_variants(variants)` | Format variant list for display |

## Cross-References

- [core/models.md](../../core/models.md) — Sample and Variant data models
- [variants.md](../../variants.md) — Variant utilities
