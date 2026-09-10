# Variant Utilities

> Module: `src/core/variants.py`

Centralized variant, position, IUPAC, formatting, and interval helpers.

## Position Contract

The canonical position type is defined in `src.core.models` and re-exported here for convenience:

```python
Position = int | str
```

Rules:

- Base positions are `int`, e.g. `73`, `309`, `16569`.
- Insertion positions are `str`, e.g. `"309.1"`, `"217.10"`.
- Use `validate_position()` at external boundaries when float precision loss must be rejected.
- Use `normalize_position()` when legacy/backward-compatible input may include floats or formatted strings.

## Position Helpers

| Function | Description |
|---|---|
| `pos_base(pos)` | Extract integer base position from int, float, or insertion string |
| `normalize_position(pos)` | Normalize base positions to `int`; preserve insertion decimals as `str` |
| `pos_sort_key(pos)` | Sort by base position, then insertion index; base sorts before insertions |

Examples:

```python
normalize_position("309")      # 309
normalize_position("309.1")    # "309.1"
pos_base("309.1")              # 309
pos_sort_key("217.10")         # (217, 10)
```

## Variant Classification and Statistics

| Function | Description |
|---|---|
| `split_variant_types(variants)` | Normalize positions and split variants into SNPs, insertions, deletions |
| `statistic_variants(variants, hv_seq_refs, hv_seqs)` | Count SNPs, insertions, deletions, identical bases, and unread bases |

`split_variant_types()` expects dicts with `pos`, `ref`, and `seq` keys and returns sorted lists: `(snps, insertions, deletions)`.

Variant type rules:

| Rule | Type |
|---|---|
| `ref == "-"` and `seq != "-"` | insertion |
| `ref != "-"` and `seq == "-"` | deletion |
| otherwise | SNP |

## IUPAC Helpers

| Function | Description |
|---|---|
| `nucleotide_to_iupac(bases)` | Convert one or more bases/IUPAC codes to a single IUPAC code |
| `iupac_to_bases(iupac)` | Expand an IUPAC code to concrete bases |

Unknown or unsupported combinations return `"N"` from `nucleotide_to_iupac()`; unknown IUPAC inputs return an empty list from `iupac_to_bases()`.

## Formatting Helpers

| Function | Description |
|---|---|
| `format_variants_simplified(variants_list)` | Format variant dicts as display strings like `73G 315.1C`; empty list returns `None` |
| `format_variants_list(variants_list, separator=" ", empty_value="None")` | Join preformatted variant strings using a separator |
| `parse_variant_position_allele(variant)` | Split strings like `152C`, `315.1C`, `524DEL` into `(position, allele)` |

`format_variants_list()` treats `separator="space"` and `separator="comma"` as named modes. Other separator values are used literally.

## Region and Interval Helpers

| Function | Description |
|---|---|
| `get_hv_region_for_position(pos, regions=None)` | Return the HV region containing a position, or `None` |
| `get_overlap(a, b)` | Return overlapping intervals between two interval lists |
| `is_position_in_intervals(pos, intervals)` | Return whether a position falls inside any interval |
| `is_special_position(pos, ref=None, seq=None)` | Return whether a concordance-excluded special position applies |

Special positions are `455`, `463`, `573`, and `309`, including insertion positions like `309.1`. When `ref` and `seq` are provided, only INDELs at those positions are treated as special.

## Conversion Helper

### `variant_to_dict(variant) -> dict | None`

Converts a `Variant` model or compatible dict to the canonical simple variant dict:

```python
{
  "pos": 309,
  "ref": "C",
  "seq": "T",
  "file": ["trace.ab1"],
  "quality": [42],
  "peaks": [[...]],  # optional
}
```

Dict input may use `position` instead of `pos`, and `alt` instead of `seq`. Variants with no position return `None` and log a warning.

## Cross-References

- [models.md](models.md) — canonical `Position` and `Variant` models
- [sample.md](sample.md) — serialization and sequence generation consumers
- [comparison.md](comparison.md) — concordance and special-position handling consumers
