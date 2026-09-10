# Variants & Positions

> Canonical variant and position domain utilities for the mtDNA pipeline.

## Purpose

`src/core/variants.py` is the **single source of truth** for all variant-related and position-related operations. Other modules must import from here (or from `src.core.models` for the `Position` type) rather than redefining these operations.

## Pipeline Position

Variants operates at the center of the pipeline — every module that works with variant data uses it.

```
src/core/models.py (Position type)
        │
        ▼
src/core/variants.py ◄── src/tools/tracy/
       │    ▲           src/tools/blastn/
       │    │           src/core/batch.py
       ▼    │           src/core/comparison.py
src/core/sample.py     src/core/flagging.py
                       src/core/mtdna_merger.py
```

## Public API

### Position Type & Helpers

| Function | Signature | Description |
|----------|-----------|-------------|
| `pos_base` | `(pos: float | str) → int` | Extract integer base position. `pos_base(309)` → `309`, `pos_base("309.1")` → `309` |
| `normalize_position` | `(pos: Union[int,float,str]) → Union[int,str]` | Normalize to canonical type: `int` for bases, `str` for insertions. Handles commas, whitespace, and floats |
| `pos_sort_key` | `(pos: float | str) → tuple[int, int]` | Sort key: `(base_position, insertion_index)`. Base positions sort before their insertions |

### Sequence Generation

> Sequence generation has moved to `src/core/sample.py`. See [core/sample.md](core/sample.md) for `generate_sequence()` and `sample_to_dict()`. The older tuple-returning `gen_seq_from_variants_regions()` has been removed in favor of the structured `generate_sequence()`.

### Variant Statistics

| Function | Signature | Description |
|----------|-----------|-------------|
| `statistic_variants` | `(variants, hv_seq_refs, hv_seqs) → tuple[int, int, int, int, int]` | Count variant types by category |
| `split_variant_types` | `(variants: list) → tuple[list, list, list]` | Split variant list into SNPs, insertions, and deletions |

### Variant Formatting

| Function | Signature | Description |
|----------|-----------|-------------|
| `format_variants_simplified` | `(variants: List) → List[str]` | Format variants as simplified position strings |
| `format_variants_list` | `(variants_list: List[str], separator: str = " ", empty_value: str = "None") → str` | Format variant strings for display (space or comma separated) |
| `parse_variant_position_allele` | `(variant: str) → Tuple[str, str]` | Parse variant string like `"152C"` or `"315.1C"` into `(position, allele)` |
| `variant_to_dict` | `(variant) → Dict` | Convert Variant model object to canonical dict format `{pos, ref, seq, file, peaks?, quality}` (legacy key order; peak heights and quality coerced to `int`) |

### Region & Interval Lookup

| Function | Signature | Description |
|----------|-----------|-------------|
| `get_hv_region_for_position` | `(pos, regions=None) → Optional[str]` | Map a position to its HV region name (`"HV1"`, `"HV2"`, `"HV3"`, or `None`) |
| `is_position_in_intervals` | `(pos, intervals) → bool` | Check whether a position falls within covered reference-base intervals; insertions after an interval endpoint are outside coverage |
| `get_overlap` | `(a: List[List[int]], b: List[List[int]]) → List[List[int]]` | Compute overlapping intervals between two interval sets |
| `is_special_position` | `(pos, ref=None, seq=None) → bool` | Check if position should be ignored in concordance (polyC and special positions) |

## Key Concepts

### Position Contract

The mtDNA pipeline uses a strict position type system:

- **Base positions**: `int` — e.g., `73`, `309`, `16569`
- **Insertion positions**: `str` — e.g., `"309.1"`, `"217.10"`

This is defined in `src/core/models.py` as the `Position` type alias. The contract is:
1. All external inputs (JSON files, Excel, CLI args) must be normalized via `validate_position()` or `normalize_position()`
2. `float` positions are **rejected** at validation boundaries to prevent precision loss
3. `pos_base()` extracts the integer base for comparisons
4. `pos_sort_key()` provides deterministic ordering: `(309, -1)` for base, `(309, 1)` for first insertion

### Sequence Generation

`generate_sequence()` (in `src/core/sample.py`) is the canonical HV consensus generator, called inside `sample_to_dict()` by Batch. ETL modules must not call it directly — see [core/sample.md](core/sample.md). The older tuple-returning `gen_seq_from_variants_regions()` has been removed.

### Variant Dict Format

The canonical variant dict uses simple keys:
```json
{
  "pos": 309,
  "ref": "C",
  "seq": "Y",
  "file": ["sample1_HV1F.ab1"],
  "quality": [35.0],
  "peaks": [[120, 80, 60, 90]]
}
```

## Configuration

- `src/config.py` — Region definitions (`REGIONS`, `POLYC_REGIONS`)
- `src/core/models.py` — `Position` type, `Variant` and `Sample` models

## Cross-References

- [ARCHITECTURE.md](ARCHITECTURE.md) §1.1 — Pipeline steps and sequence generation ownership
- [core/models.md](core/models.md) — Position type and Variant model definitions
- [core/batch.md](core/batch.md) — Batch calls `sample_to_dict()` (which runs `generate_sequence()` and `statistic_variants()`)
- [core/flagging.md](core/flagging.md) — Uses `pos_base`, `normalize_position`, `is_special_position`
