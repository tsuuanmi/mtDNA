# mtDNA Merger

> Single-pass merger for per-sample mtDNA JSON profiles with strict validation.

## Purpose

`src/core/mtdna_merger.py` provides `MtDnaMerger` for merging multiple per-sample JSON profiles of the same sample and `merge_by_regions()` for region-aware unification. It enforces strict consistency checks and uses `Sample` objects internally for compatibility with the rest of the pipeline.

Merged results are marked with `source_tool = "unified"` (via `Tool.UNIFIED`) and include a `region_sources` dict in the `information` dict tracking which tool produced each region's data (e.g., `{"HV1": "sequencher", "HV2": "sequencher", "HV3": "mutation_surveyor"}`).

## Pipeline Position

```
Multiple per-sample JSON files (from different tools)
         │
         ▼
   MtDnaMerger.merge() / merge_by_regions()
    ├── Load each profile into flat variant dict
    ├── Validate individual profiles (ref allele, intervals)
    ├── Auto-merge variants (hard error on conflicts)
    ├── Union intervals across profiles
    ├── Validate merged profile (overlap consistency)
    ├── Set source_tool = Tool.UNIFIED, region_sources = {"HV1": "...", "HV2": "...", "HV3": "..."}
    └── Output merged Sample → JSON
```

## Public API

### `MergeValidationError`

```python
class MergeValidationError(ValueError)
```

Raised when merged variants fail validation checks or contain overlapping conflicts.

### `MergeResult`

```python
class MergeResult(sample: Sample, warnings: list[str] | None = None)
```

Result of the auto-merge phase. Contains the merged `Sample` object and a list of warning strings from validation.

### `MtDnaMerger`

```python
MtDnaMerger(ref_path: str | Path | None = None)
```

Single-pass merger for per-sample mtDNA JSON profiles. If `ref_path` is not provided, uses `Settings.directories.ref / "rCRS.fasta"`.

#### Methods

| Method | Description |
|--------|-------------|
| `merge(json_files: list[str | Path], sample_id: str, output_dir: str | Path = ".") → MergeResult` | Merge multiple JSON profiles into a single unified profile, writing `{sample_id}/{sample_id}.json` |

#### Methods (Additional)

| Method | Description |
|--------|-------------|
| `build_merged_sample(merged_data, sample_id, ref_path) → Sample` | Build a Sample from merged data (variants, intervals, flags, info) |
| `merge_profiles(flat_variants_list) → dict[Position, dict]` | Auto-merge variant dicts from multiple profiles, resolving conflicts |

#### Private Methods

| Method | Description |
|--------|-------------|
| `_load_profile(path: Path, _sample_id: str) → dict` | Load one JSON file in unwrapped per-sample format |
| `_flat_variant_dict(profile)` | Extract flat variant dict from a profile |
| `_union_intervals(profiles)` | Compute union of intervals across profiles |
| `_load_rCRS_reference()` | Load rCRS reference genome for validation |
| `_validate_individual_profile(profile, filename)` | Validate single profile (ref allele, interval bounds) |
| `_check_overlap_consistency(merged, profiles)` | Check explicit overlap conflicts (implicit silence is not treated as a conflict) |
| `_validate_merged_profile(merged, profiles)` | Validate the merged profile |
| `_load_and_flatten_profiles(json_files)` | Load and flatten multiple JSON profiles into variant dicts |
| `_collect_flags_and_info(profiles, source_tools)` | Collect sample flags and info from multiple profiles |

### `merge_samples()`

```python
merge_samples(samples: list[Sample], ref_path: str | Path | None = None, output_dir: str | Path = ".") → MergeResult
```

Convenience function for merging `Sample` objects directly. Creates a `MtDnaMerger` instance, merges variant dicts, intervals, flags, and information, and writes `{output_dir}/{sample_id}/{sample_id}.json`.

### CLI Entry Point

```bash
python -m src.core.mtdna_merger \
    --files FILE1 FILE2 ... \
    --sample-id ID \
    --output-dir DIR \
    [--ref-path PATH]
```

Output: `{output_dir}/{sample_id}/{sample_id}.json` (unwrapped per-sample format).

## Key Concepts

### Input Format

`_load_profile` accepts both wrapped and unwrapped per-sample JSON formats:

- **Wrapped**: `{sample_id: {no_snps, HV1, variants, ...}}` — top-level key matches `--sample-id`
- **Unwrapped**: `{no_snps, HV1, variants, ...}` — profile data is directly at top level (no sample_id wrapper)

When the top-level key matches `--sample-id`, the file is treated as wrapped; otherwise, it is treated as unwrapped and the entire JSON is used as the profile dict.

### Strict Validation

The merger enforces strict consistency (no non-strict mode currently):

1. **Interval bounds**: All variants must lie inside their own profile's declared intervals
2. **Ref allele consistency**: Reference alleles must exactly match the rCRS reference
3. **Overlap consistency**: When multiple profiles explicitly report different alleles at the same position, alleles must agree (explicit conflicts detected). When only one profile reports a variant at a position, the explicit call is trusted — the other tool may simply not detect/call that variant (e.g., different deletion/insertion calling sensitivity).

### Auto-Merge Strategy

When multiple profiles contain variants at the same position:
- If alleles agree → keep the variant, union `file[]` and `quality[]`
- If alleles disagree → `MergeValidationError` (hard error)
- If only one profile reports a variant → trust the explicit call (the other tool may not detect it)

> **Note**: Quality-based conflict resolution is not yet implemented. The Regenerate Gate (see ARCHITECTURE.md) is planned to add strategies such as `best_quality` and `prefer_automated`.

### Source Tool Tracking

Merged samples use `Tool.UNIFIED` (`"unified"`) as their `source_tool`, regardless of which tools contributed. The `region_sources` key in the `information` dict tracks which tool produced each region's data:

```python
Sample(
    source_tool=Tool.UNIFIED,  # always "unified" for merged results
    information={
        "source_tool": "unified",
        "region_sources": {"HV1": "mutation_surveyor", "HV2": "sequencher", "HV3": "sequencher"},  # per-region provenance
        ...
    }
)
```

### HV Sequence Contract

The merger rebuilds HV strings from the merged flat variant dict via `generate_sequence()` (called inside `sample_to_dict()`), so no single input file dominates the output. This is consistent with the Batch class's approach.

## Configuration

- `src/config.py` — `Settings.directories.ref` for rCRS reference path
- `src/core/models.py` — `Sample`, `Variant`, `Tool` (including `Tool.UNIFIED`)
- `src/core/sample.py` — `sample_to_dict()` for output serialization
- `src/core/sample.py` — `generate_sequence()` for HV sequence generation (called inside `sample_to_dict()`)
- `src/core/variants.py` — `variant_to_dict()`, `split_variant_types()` for variant normalization

## Cross-References

- [regenerate.md](regenerate.md) — interval-scoped final JSON generation
- [models.md](models.md) — `Sample`, `Variant`, `Tool` models
- [sample.md](sample.md) — `sample_to_dict()` for output format
- [variants.md](variants.md) — Position normalization and interval utilities
- [unify.md](unify.md) — batch-level region unification entry point

---

## Region-Filtered Merge Pipeline

The region-filtered merge pipeline combines Mutation Surveyor and Sequencher per-region results into unified profiles. Sequencher writes per-region JSONs for all regions it has data for; Mutation Surveyor writes per-region JSONs for auto-pass regions only. When both tools cover the same region, **Sequencher takes priority** (manual review is more reliable than auto-pass, which can miss variants).

### Pipeline Flow

```
Step 2: MS core pipeline → per-sample JSONs (unwrapped)
         │
         ▼
Step 5 (inline): auto-pass region filter (inside etl.build_samples)
    ├── For each MS sample (after ETL):
    │   ├── Read sample_flags → determine auto-pass regions
    │   │   ("Autopass HV1" → HV1, "Autopass HV2 and HV3" → HV2+HV3)
    │   ├── Filter MS sample → keep only auto-pass regions
    │   └── No auto-pass or no variants → skip (not written to batch)
    │
    └── Output: filtered Samples written to batch JSON (no separate CLI step)

Step 5b: Sequencher per-region JSON output (inline in Sequencher pipeline)
    ├── For each Sequencher sample:
    │   ├── Detect region group from TXT filename
    │   ├── Write per-region JSON for ALL regions Sequencher has data for
    │   └── (No complementary filtering — Sequencher writes everything it has)
    │
    └── Output: per-region JSONs in Sequencher results dir
         │
         ▼
Step 6: `python -m src.core.unify` (batch-level)
    ├── Discover per-region JSONs under MS and Sequencher `regions/` dirs
    ├── Group by sample ID and region group (`HV1`, `HV2-3`)
    ├── For each sample, call merge_by_regions()
    └── Output: unified per-sample JSONs plus statistic_fullbatch.json
```

### `filter_sample_by_regions()`

```python
filter_sample_by_regions(sample: Sample, keep_regions: List[str]) → Optional[Sample]
```

Tool-agnostic post-ETL region filter. Takes a Sample and a list of region names to keep (e.g., `["HV1"]`, `["HV2", "HV3"]`). Returns a new Sample with only variants whose positions fall within the kept regions' intervals. Returns `None` when no matching intervals exist.

- Variant positions are checked against kept-region intervals using `pos_base()`
- The returned `Sample.intervals` dict contains **only** kept region keys (Shape A)
- HV sequences for dropped regions are set to `None`
- `source_tool` is preserved (filtering does not change tool attribution)
- `sample_flags`, `variant_flags`, `information`, and `batch_id` are carried forward

### Auto-Pass Region Detection

| Sample Flag | Regions Kept (MS) | Regions from Sequencher |
|-------------|-------------------|---------------------------|
| `"Autopass HV1"` | HV1 | HV1, HV2, HV3 (all regions) |
| `"Autopass HV2 and HV3"` | HV2, HV3 | HV1, HV2, HV3 (all regions) |
| Both flags | HV1, HV2, HV3 | HV1, HV2, HV3 (all regions) |
| Neither flag | None (no MS output) | HV1, HV2, HV3 (all regions) |

### Batch Unify CLI

```bash
python -m src.core.unify \
  --ms-dir results/tools/mutation_surveyor/<batch> \
  --seq-dir results/tools/sequencher/<batch> \
  --output-dir results/merged/<batch> \
  --batch-id <batch> \
  --ref-path ref/rCRS.fasta
```

The unify step can merge samples that have one or more discovered region groups. The merged output is marked `source_tool = "unified"` and includes `information.region_sources` for the regions that contributed data.

### Constraints

- When both tools cover the same region, Sequencher takes priority (manual review > auto-pass); MS variants for that region are discarded, `region_sources` shows only `"sequencher"`
- Sequencher writes per-region JSONs for all regions it has data for (no complementary filtering)
- MS writes per-region JSONs for auto-pass regions only (unchanged)
- Filter produces no file when nothing to keep (pipeline checks file existence)
- ETL modules stay pure — filtering is post-ETL
- Region filtering uses whole regions (HV1, HV2, HV3), not sub-ranges

## Region-Based Merge

The `merge_by_regions()` function provides per-region provenance tracking
with a **Sequencher-wins priority rule**. It accepts a dict mapping region
group names to lists of Sample objects, merges them into a unified Sample,
and records which tool produced each region's data in `region_sources`.

When both Sequencher and Mutation Surveyor provide data for the same
region group, Sequencher takes priority because manual review is more
reliable than MS auto-pass (which can miss variants). The winning tool's
variants are used; the other tool's variants for that region are discarded.
`region_sources` shows only the winning tool name, not a comma-separated list.

Tool priority order (lower number = higher priority):
1. `sequencher` — manual review, most reliable
2. `mutation_surveyor` — auto-pass, can miss variants
3. `tracy` — variant caller
4. `blastn` — alignment-based

```python
from src.core.mtdna_merger import merge_by_regions

# Example: Both tools cover HV2-3 — Sequencher wins
result = merge_by_regions(
    region_samples={
        "HV2-3": [ms_hv23_sample, seq_hv23_sample],  # Both tools → Sequencher wins
        "HV1": [seq_hv1_sample],                       # Sequencher only
    },
    sample_id="LID_001",
    ref_path="ref/rCRS.fasta",
    output_dir=output_dir,
)
```

The resulting `information` dict contains — note `region_sources` shows only
the winning tool per region:

```json
{
  "source_tool": "unified",
  "region_sources": {
    "HV2": "sequencher",
    "HV3": "sequencher",
    "HV1": "sequencher"
  }
}
```

### Per-Region Intermediate JSONs

Per-region JSONs are **temporary intermediate files** written by each tool
pipeline. They use the naming convention `HV1_{LID}.json` and
`HV2-3_{LID}.json` under region-group subdirectories:

```text
{tool_dir}/regions/HV1/HV1_{LID}.json
{tool_dir}/regions/HV2-3/HV2-3_{LID}.json
```

The unify step discovers these per-region JSONs and merges them via
`merge_by_regions()`. The final output is a single unified `{LID}.json`
per LID in the standard format — no change to the downstream format.

See [region.md](region.md) for the per-region JSON I/O and
range/interval QC validation module.
