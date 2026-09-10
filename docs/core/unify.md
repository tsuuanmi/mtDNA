# Batch Unification

> Module: `src/core/unify.py`

Merges Mutation Surveyor and Sequencher per-region JSONs into unified mtDNA profiles.

## Purpose

`unify_batch()` discovers per-region intermediate JSON files from Mutation Surveyor and Sequencher, merges them by region via `merge_by_regions()`, and writes unified per-sample JSON plus `statistic_fullbatch.json`.

This is region-aware unification, not simple file concatenation. Region provenance is preserved through `information.region_sources`.

## Expected Input Layout

Each tool directory is expected to contain a `regions/` subtree:

```text
<tool_dir>/
└── regions/
    ├── HV1/
    │   └── HV1_<LID>.json
    └── HV2-3/
        └── HV2-3_<LID>.json
```

The region filename patterns come from `src.core.region.REGION_FILE_PATTERNS`:

| Region group | Pattern |
|---|---|
| `HV1` | `HV1_{lid}.json` |
| `HV2-3` | `HV2-3_{lid}.json` |

## Public API

### `unify_batch(ms_dir, seq_dir, output_dir, batch_id, ref_path=None) -> int`

| Argument | Description |
|---|---|
| `ms_dir` | Mutation Surveyor directory containing per-region JSONs |
| `seq_dir` | Sequencher directory containing per-region JSONs |
| `output_dir` | Destination directory for unified outputs |
| `batch_id` | Batch identifier for logging and batch output |
| `ref_path` | Optional rCRS FASTA path; defaults to settings reference path |

Returns the number of unified samples produced.

## Processing Steps

1. Resolve the reference FASTA path from `ref_path` or settings.
2. Discover all matching region JSON files under both tool directories.
3. Group discovered files by sample ID, region group, and tool name.
4. Load each region JSON with `read_region_json()`.
5. Merge each sample with `merge_by_regions()`.
6. Write each unified per-sample JSON during merge.
7. Write `statistic_fullbatch.json` for all unified samples.

## Output Layout

```text
<output_dir>/
├── statistic_fullbatch.json
├── SAMPLE_1/SAMPLE_1.json
└── SAMPLE_2/SAMPLE_2.json
```

Each unified sample uses `source_tool = "unified"` and includes per-region source provenance in `information.region_sources`.

## CLI

```bash
python -m src.core.unify \
  --ms-dir results/tools/mutation_surveyor/<batch> \
  --seq-dir results/tools/sequencher/<batch> \
  --output-dir results/merged/<batch> \
  --batch-id <batch> \
  --ref-path ref/rCRS.fasta
```

Arguments:

| Flag | Required | Description |
|---|---:|---|
| `--ms-dir` | Yes | Mutation Surveyor per-region output directory |
| `--seq-dir` | Yes | Sequencher per-region output directory |
| `--output-dir` | Yes | Unified output directory |
| `--batch-id` | Yes | Batch identifier |
| `--ref-path` | No | Reference FASTA override |

## Cross-References

- [region.md](region.md) — per-region JSON I/O and filename patterns
- [mtdna_merger.md](mtdna_merger.md) — region-aware merge behavior
- [batch.md](batch.md) — `statistic_fullbatch.json` writing
