# Pairwise Sample Comparison

> Tool-agnostic pairwise comparison of Sample objects with binary concordance and flagging.

## Purpose

`src/core/comparison.py` provides `SampleComparator` for comparing two `Sample` objects (typically from different analysis tools) to determine concordance and identify unique variants. It produces `ComparisonResult` objects and exports comparison data to TSV/Excel formats.

## Pipeline Position

```
Two Sample objects (e.g., Sequencher vs Mutation Surveyor)
         │
         ▼
   SampleComparator.compare()
    ├── Find unique variants in each
    ├── Determine concordance (Y/N/N/A)
    ├── Compute flag levels
    └── Produce ComparisonResult
         │
         ▼
   TSV/Excel export + JSON output
```

## Public API

### `SampleComparator`

```python
SampleComparator  # All methods are static — stateless comparison
```

#### Core Comparison Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `compare` | `(result_a: Sample, result_b: Sample) → ComparisonResult` | Full pairwise comparison of two samples, returning `ComparisonResult` with per-tool data and unique variants |
| `is_concordant` | `(result_a: Sample, result_b: Sample) → str` | Binary concordance: `"Y"`, `"N"`, or `"N/A"` |
| `find_unique_variants` | `(variants_a, variants_b, *, ignore_special=False) → list[Variant]` | Find variants in `variants_a` not present in `variants_b`. If `ignore_special=True`, skip special positions |
| `format_intervals` | `(intervals: dict[str, list[list[int]]] | None) → str` | Format intervals dict to a readable string. Returns `"FULL REGION"` if all predefined regions are covered |
| `compare_batches` | `(batch_a_dir, batch_b_dir, ...) → None` | Compare all samples in two batch directories |
| `summary_row` | `(result: ComparisonResult, ...) → Dict` | Generate a summary row dict for TSV/Excel output |

### Module-Level Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `create_variant_level_data` | `(result: ComparisonResult, batch_id: str = "") → list[dict[str, Any]]` | Create per-variant level data rows for TSV/Excel |
| `write_comparison_json` | `(batch_result: dict[str, Any], output_path: str | Path) → Path` | Write comparison batch result to JSON file |
| `write_comparison_tsv` | `(batch_result: dict[str, Any], output_path: str | Path) → Path` | Write comparison summary to TSV file |
| `write_variant_level_tsv` | `(batch_result: dict[str, Any], output_path: str | Path) → Path` | Write variant-level comparison data to TSV |
| `write_comparison_excel` | `(batch_result: dict[str, Any], output_path: str | Path) → Path` | Write comparison data to Excel file with summary and variant sheets |
| `compare_batch_files` | `(batch_a_dir, batch_b_dir, ...) → None` | CLI entry point for batch comparison |

## Key Concepts

### Concordance Logic

`is_concordant()` determines binary concordance between two samples:

1. If either sample has empty variants → `"N/A"`
2. Filter out special positions (polyC, always-concordant positions)
3. Check all shared positions for disagreements → `"N"` if any differ
4. Check unique variants in each direction:
   - A-only variants in HV regions → `"N"`
   - B-only variants in A's analyzed ranges → `"N"`
5. Otherwise → `"Y"`

### Always-Concordant Positions

The following positions are always treated as concordant regardless of content:
- `573.1`, `309.1`, `309.2`

### Range Filtering

B-only variants that fall outside A's analyzed ranges are ignored in concordance calculation, matching the original behavior where Sequencher-only variants outside the pipeline's ranges are excluded.

### Comparison Output

The `ComparisonResult` model contains:
- `tool_a` and `tool_b`: `ToolResult` objects with variants, flags, analyzed ranges
- `unique_variants_a` and `unique_variants_b`: Variants unique to each tool
- `concordance`: `"Y"`, `"N"`, or `"N/A"`

The TSV/Excel summary and variant-level column headers are parameterized by each batch's source tool (e.g. `Variants (sequencher)`, `Variant Flags (mutation_surveyor)`) instead of the generic `(A)`/`(B)` convention. A batch without tool attribution (empty batch or `Tool.UNKNOWN`) displays `Unknown`. `compare_batches` exposes the resolved names as `tool_name_a`/`tool_name_b` and the column lists as `summary_fields`/`variant_fields` in the returned result dict.

### Legacy Output Mode

`compare_batches` / `compare_batch_files` accept an opt-in `legacy_format` flag (CLI: `--legacy-format`) that reproduces the legacy 13-column comparison summary layout used by the older reviewing system (see `Batch_MS_280326_001.tsv`). When enabled, the summary TSV, Excel, and JSON outputs use the legacy layout; the variant-level TSV is unchanged. The default 16-column layout is byte-identical when the flag is absent.

Legacy columns (hardcoded labels, in order): `Sample ID`, `Batch`, `Analyzed Range (Pipeline)`, `Variants (Pipeline)`, `Variants Unique (Pipeline)`, `Analyzed Range (Sequencher)`, `Variants (Sequencher)`, `Variants Unique (Sequencher)`, `Variants Count (Sequencher)`, `Concordant`, `Flag`, `Flagged (Pipeline)`, `Flagged (Sequencher)`.

The legacy `Flagged (<tool>)` column merges each tool's sample and variant flags into one deduplicated string joined by `", "` (`No` when empty). Variant-style flags (per-variant reasons, then the consolidated `16180-16193 region (...)` entry) come first; sample-style flags (`No 315.1 variant`, `Has-... variant`, `Consecutive indels`, `Lowercase`, range flags) come last. Insertion flags keep their current `Insertion {seq} at {pos}` form. Labels are hardcoded to `Pipeline` (A) / `Sequencher` (B) regardless of the actual source tools, and only the `Variants Count (Sequencher)` count column is emitted.

`compare_batches` exposes the selected layout via `legacy_format` and the `summary_fields`/`summary_rows` keys in the returned result dict.

## Configuration

- `src/config.py` — `REGIONS` for interval formatting and range checks
- `src/core/models.py` — `Sample`, `Variant`, `ComparisonResult`, `ToolResult`
- `src/core/variants.py` — Position normalization, `is_special_position`, `variant_to_dict`

## Cross-References

- [core/models.md](models.md) — `ComparisonResult`, `ToolResult`, `Sample` model definitions
- [core/flagging.md](flagging.md) — Flag levels and severity classification
- [variants.md](variants.md) — Position utilities and `is_special_position`
