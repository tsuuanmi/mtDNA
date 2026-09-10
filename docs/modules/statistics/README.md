# Statistics Module

> Query and filter mtDNA variant data with multi-filter support.

## Purpose

`src/modules/statistics/` provides utilities for querying, aggregating, and summarizing mtDNA variant data. `query_multi_filters.py` supports multi-filter queries with region-based windowing and variant formatting. `merge.py` consolidates data and per-batch results into merged outputs. It also contains the blinded TNLS test-data helpers: `blind_copy.py`, `add_variants_to_samples.py`, and `add_variants_to_metadata.py`.

Run the merger with `bash scripts/modules/statistics/run.sh`.

Prepare the blinded TNLS test data with `bash scripts/modules/statistics/test.sh`; it recreates `/home/tan/workspaces/mtdna_raw/temp/test`. Maintain its local, ignored blind-ID mapping at `data/modules/statistics/sample_mapping.tsv`.

## Public API

| Function | Description |
|----------|-------------|
| `variants_in_windows(variants, regions)` | Group variants by genomic region windows |
| `format_analyzed_range(intervals)` | Format analyzed intervals for display |
| `format_variants(variants)` | Format variant list for display |

## Cross-References

- [core/models.md](../../core/models.md) — Sample and Variant data models
- [variants.md](../../variants.md) — Variant utilities
