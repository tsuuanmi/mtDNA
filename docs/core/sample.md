# Sample Serialization

> Canonical serialization and deserialization between `Sample` objects and on-disk JSON format.

## Purpose

`src/core/sample.py` provides the bridge between in-memory `Sample` model objects and the on-disk JSON format used throughout the pipeline. It handles both grouped variant format (from ETL) and flat variant format (from legacy scripts).

## Pipeline Position

```
ETL Module (tool-specific) → Sample (Pydantic)
         │
         ▼
   sample_to_dict() → canonical JSON
         │
         ▼   load_sample_batch() ← on-disk JSON
   Batch / Comparison / FASTA Gen
```

## Public API

### `sample_to_dict(sample: Sample, ref_path: str = "ref/rCRS.fasta") → dict[str, Any]`

Convert a `Sample` object to the canonical on-disk JSON dict format.

**Steps:**
1. Call `generate_sequence()` to compute HV sequences from `sample.intervals`
2. Call `statistic_variants()` to compute variant counts (SNPs, insertions, deletions, identicals, unread)
3. Convert variant list to grouped format: `{snps: [...], insertions: [...], deletions: [...]}`
4. Call `flag_variants()` for per-variant flags
5. Recompute core sample-level flags against the serialized variant set with `recompute_sample_flags()`
6. Build `information` with `source_tool`, `batch_id`, additional metadata, and optional `computation_error`
7. Return dict ready for JSON serialization

**Output format:**
```json
{
  "no_snps": 5,
  "no_ins": 1,
  "no_dels": 0,
  "no_identicals": 3,
  "no_unread": 0,
  "HV1": "ACGT...",
  "HV2": "ACGT...",
  "HV3": "ACGT...",
  "variants": {"snps": [...], "insertions": [...], "deletions": [...]},
  "intervals": {"HV1": [[16024, 16365]], ...},
  "sample_flags": ["No 315.1 variant"],
  "variant_flags": {"73|A|G": ["Low quality"], ...},
  "information": {"batch_id": "...", "source_tool": "sequencher", ...}
}
```

### `load_sample_batch(path: str | Path) → list[Sample]`

Load a batch JSON file back into `Sample` objects. Handles:
- **Type-grouped format**: `{snps: [...], insertions: [...], deletions: [...]}` — flattened to `List[Variant]` (with `variant_type` tags)
- **Region-grouped format**: `{HV1: [...], HV2: [...], HV3: [...]}` — flattened to `List[Variant]` (legacy/manual-regenerate JSON)
- **Flat variant format**: `[{pos, ref, seq, ...}]` — normalized directly

Also handles legacy key mapping:
- `HV1` → `hv1`, `HV2` → `hv2`, `HV3` → `hv3`
- `alt` → `seq` (bidirectional compatibility)

### `load_sample(path: str | Path) → Sample`

Load a single per-sample JSON file (`<json_dir>/<sample_id>/<sample_id>.json`) into one `Sample` object; the file stem is used as the `sample_id`. Shares the same normalization as `load_sample_batch()` (type-grouped and region-grouped variants, flat variants, legacy key mapping).

### `_normalize_keyed_sample(raw: dict[str, Any]) → dict[str, Any]`

Normalize a keyed sample entry from grouped variant format to a flat list. Accepts type-grouped variants `{snps, insertions, deletions}` (each tagged with `variant_type`), region-grouped variants `{HV1, HV2, HV3}`, or an already-flat list.

### `_parse_variant_dict(v: dict[str, Any], sample_id: str) → dict[str, Any]`

Normalize a single variant dict for `Sample` model validation:
- Ensures `pos` field exists (falls back to `position`)
- Maps `alt` ↔ `seq` for backward compatibility
- Normalizes `quality`: scalar values become single-item lists; lists are preserved; missing/`None` becomes `[]`
- Sets defaults for missing fields

## Key Concepts

### Format Bridge

The module bridges two representations:
1. **Internal**: `Sample` Pydantic model with typed `List[Variant]`
2. **On-disk**: JSON with grouped variants `{snps, insertions, deletions}`

`sample_to_dict()` converts internal → on-disk. `load_sample_batch()` converts on-disk → internal.

### Error Handling

If sequence generation or statistics computation fails in `sample_to_dict()`:
- `information.computation_error` is set to `true`
- HV sequences default to empty strings
- Variant counts default to zeros
- Processing continues for remaining samples

This ensures downstream consumers can detect incomplete data rather than silently operating on defaults.

## Cross-References

- [core/models.md](models.md) — `Sample` and `Variant` model definitions
- [core/batch.md](batch.md) — `Batch.write()` calls `sample_to_dict()`
- [core/flagging.md](flagging.md) — `flag_variants()` used in serialization
- [variants.md](variants.md) — `statistic_variants()` and position utilities

---

### `filter_sample_by_regions(sample: Sample, keep_regions: list[str]) → Sample | None`

Tool-agnostic post-ETL region filter. Returns a new `Sample` containing only variants and intervals for the specified regions, or `None` when no matching intervals exist.

**Usage:**
```python
from src.core.sample import filter_sample_by_regions

# Keep only HV1 (e.g., when HV1 auto-passes in Mutation Surveyor)
filtered = filter_sample_by_regions(sample, ["HV1"])

# Keep HV2+HV3 (e.g., complementary regions for Sequencher)
filtered = filter_sample_by_regions(sample, ["HV2", "HV3"])

# Empty list returns None
result = filter_sample_by_regions(sample, [])  # → None
```

**Behaviour:**
- Variants are kept when their position (using `pos_base()` for insertions) falls within any interval of a kept region
- The returned `Sample.intervals` dict contains **only** kept region keys (omitted if empty — Shape A)
- HV sequences for dropped regions are set to `None`; kept regions retain their value
- `source_tool` is **preserved** — filtering does not change tool attribution (only merging sets `Tool.UNIFIED`)
- `sample_flags`, `variant_flags`, `information`, and `batch_id` are carried forward unchanged
- Raises `ValueError` if any region name is not in `{"HV1", "HV2", "HV3"}`

**Pipeline role:** Used by per-region JSON writing and regenerate/unify paths to scope a `Sample` to selected HV regions after ETL.

### Region-Level Provenance: `region_sources`

When samples are merged by region via `merge_by_regions()`, the `information` dict in the serialized output includes `region_sources` — a mapping of each HV region to the source tool that contributed data for that region:

```json
{
  "information": {
    "source_tool": "unified",
    "batch_id": "MS_210426_003",
    "region_sources": {
      "HV1": "mutation_surveyor",
      "HV2": "sequencher",
      "HV3": "sequencher"
    }
  }
}
```

`region_sources` replaces the earlier `merged_from` field (sample-level tool list) with per-region provenance tracking. Each region key maps to the winning tool name (not a list). For single-source samples, all regions map to that tool.

**Pipeline role:** Set by `merge_by_regions()` in `src/core/mtdna_merger.py`. Consumed by `unify_batch()` in `src/core/unify.py` for logging and debugging.

**Serialization:** `region_sources` is included in the `information` dict written by `sample_to_dict()`. It is read back by `_normalize_keyed_sample()` during deserialization but does not affect the `Sample` model directly — it remains in `information`.
