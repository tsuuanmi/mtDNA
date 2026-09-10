# Batch Unification

> Merge Mutation Surveyor and Sequencher per-region JSONs into unified profiles.

## Purpose

`src/core/unify.py` provides the `unify_batch()` function and CLI for discovering per-region JSON files from Mutation Surveyor and Sequencher directories, merging them by region via `merge_by_regions()`, and producing a single unified `{sample_id}.json` per LID plus a `statistic_fullbatch.json` batch summary.

## Pipeline Position

```
Mutation Surveyor per-region JSONs (HV1_{LID}.json, HV2-3_{LID}.json)
Sequencher per-region JSONs (HV1_{LID}.json, HV2-3_{LID}.json)
         │
         ▼
   unify_batch()
    ├── Discover per-region JSONs for each LID
    ├── Merge by region via merge_by_regions() (Sequencher-wins priority)
    ├── Produce unified {LID}.json per sample
    └── Write statistic_fullbatch.json
```

## Public API

### `unify_batch(ms_dir, seq_dir, output_dir, batch_id, ref_path=None) → int`

Merge MS and Sequencher per-region JSONs into unified profiles.

| Parameter | Type | Description |
|-----------|------|-------------|
| `ms_dir` | `Path` | Directory with Mutation Surveyor per-region JSONs |
| `seq_dir` | `Path` | Directory with Sequencher per-region JSONs |
| `output_dir` | `Path` | Base output directory for unified results |
| `batch_id` | `str` | Batch identifier |
| `ref_path` | `Path \| None` | Path to rCRS reference FASTA (defaults to settings) |

Returns the number of unified samples produced.

### CLI Entry Point

```bash
python src/core/unify.py \
    --ms-dir <mutation_surveyor_dir> \
    --seq-dir <sequencher_dir> \
    --output-dir <output_dir> \
    --batch-id <batch_id> \
    [--ref-path <reference_path>]
```

## Key Concepts

### Region File Discovery

Per-region JSON files follow the naming convention:
- `HV1_{LID}.json` — HV1 region data
- `HV2-3_{LID}.json` — HV2+HV3 region data

Files are discovered in `{tool_dir}/{LID}/` subdirectories.

### Merge Priority

When both tools provide data for the same region, Sequencher wins (manual review is more reliable than MS auto-pass). The winning tool's variants replace the other's; `region_sources` records the winning tool per region.

### Output Format

Each unified LID produces `{LID}.json` in the standard `sample_to_dict()` format with:
- `source_tool = "unified"`
- `region_sources` in `information` dict tracking per-region provenance
- `statistic_fullbatch.json` at the batch level aggregating all samples

## Configuration

- `src/config.py` — Reference FASTA path from settings
- `src/core/mtdna_merger.py` — `merge_by_regions()` for region-based merging
- `src/core/region.py` — Region file patterns and JSON I/O
- `src/core/batch.py` — Batch aggregation and serialization

## Cross-References

- [core/mtdna_merger.md](../core/mtdna_merger.md) — `merge_by_regions()` and merge logic
- [core/region.md](../core/region.md) — Per-region JSON I/O and naming conventions
- [core/batch.md](../core/batch.md) — Batch serialization
- [ARCHITECTURE.md](../ARCHITECTURE.md) — Merge pipeline overview
