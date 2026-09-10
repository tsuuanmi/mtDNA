# Output Paths

> Module: `src/core/paths.py`

Defines the standard output directory layout used by migrated tool pipelines.

## Purpose

`OutputPaths` centralizes the three output categories produced by tool pipelines:

```text
base_dir/
├── json/        # Batch.write() output
├── regions/     # per-region intermediate JSONs
└── preprocess/  # tool-specific preprocessing artifacts
```

BLASTn and Tracy use this layout. Mutation Surveyor and Sequencher are not fully migrated to it yet.

## Public API

### `OutputPaths(base_dir: Path)`

Constructs paths under a batch/tool output base directory.

| Attribute | Path | Purpose |
|---|---|---|
| `base_dir` | `<base_dir>` | Root output directory passed by the caller |
| `json_dir` | `<base_dir>/json` | Standard batch JSON output from `Batch.write()` |
| `regions_dir` | `<base_dir>/regions` | Per-region intermediate JSONs from `write_region_jsons()` |
| `preprocess_dir` | `<base_dir>/preprocess` | Tool-specific preprocessing artifacts; Tracy stores each trace set and its `qc_report.json` under `<sample_id>/` |

### `create_dirs() -> None`

Creates `json_dir`, `regions_dir`, and `preprocess_dir` with `parents=True, exist_ok=True`.

## Usage

```python
from pathlib import Path
from src.core.paths import OutputPaths

paths = OutputPaths(Path("results/tools/blastn/20260615_batch"))
paths.create_dirs()

# paths.json_dir       -> results/tools/blastn/20260615_batch/json
# paths.regions_dir    -> results/tools/blastn/20260615_batch/regions
# paths.preprocess_dir -> results/tools/blastn/20260615_batch/preprocess
```

## Pipeline Role

Tool pipelines should pass these subdirectories to the lower-level writers:

- `Batch.write(paths.json_dir, batch_id, nest_batch_id=False)`
- `write_region_jsons(sample, paths.regions_dir, ...)`
- preprocessing output goes under `paths.preprocess_dir`; Tracy writes artifacts and `qc_report.json` in `paths.preprocess_dir / <sample_id>`

This avoids mixing final JSON, per-region intermediate files, and preprocessing artifacts in one directory.

## Cross-References

- [batch.md](batch.md) — `Batch.write()` JSON output
- [region.md](region.md) — per-region JSON output
- [regenerate.md](regenerate.md) — final filtered JSON output
