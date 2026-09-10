# Per-Region Intermediate JSON I/O and Range QC

> Module: `src/core/region.py`

Handles per-region intermediate JSON I/O and range/interval QC validation
for the mtDNA pipeline.

Per-region JSONs are **temporary intermediate files** — the final output is
always a single unified JSON per LID in the standard format (see ARCHITECTURE.md).

## Naming Convention

| Region Group | Filename Pattern | Contains |
|---|---|---|
| HV1 | `HV1_{LID}.json` | HV1 variants and intervals |
| HV2-3 | `HV2-3_{LID}.json` | HV2 + HV3 variants and intervals |

Files are placed in region-group subdirectories under the supplied output directory:

```text
{output_dir}/HV1/HV1_{LID}.json
{output_dir}/HV2-3/HV2-3_{LID}.json
```

## Region Group Mappings

| Autopass Flag | Region Group | Region Keys |
|---|---|---|
| `Autopass HV1` | `HV1` | `["HV1"]` |
| `Autopass HV2 and HV3` | `HV2-3` | `["HV2", "HV3"]` |

## JSON I/O Functions

### `write_region_jsons(sample, output_dir, ref_path, region_groups)`

Write per-region intermediate JSONs for auto-pass regions. For each region
group that has intervals or variants, writes a separate JSON file containing
only that region's data. Adds `region_sources` to the `information` dict for
provenance tracking by the unify step.

Returns: `dict[str, Path]` mapping region group name → output file path.

### `read_region_json(path, sample_id)`

Read a per-region intermediate JSON file. Returns a `Sample` object.
Provenance is carried via `region_sources` in `information`.
Legacy `region_group` string fields are stripped on read.

### `region_groups_from_flags(sample_flags)`

Map auto-pass sample flags to region file groups. Returns
`Dict[str, list[str]]` mapping group names to region key lists.

## Range/Interval QC Functions

### `validate_intervals_in_region(intervals, region_group, region_keys)`

Validate that all intervals fall within the canonical region boundaries.
Returns a list of warning flag strings for out-of-range intervals.

### `validate_region_group_intervals(intervals, region_group)`

Convenience wrapper that looks up region keys from the group name.
Returns a list of warning flag strings.

## Per-Region JSON Format

Each per-region JSON is an unwrapped per-sample dict (same as current `{LID}.json`),
but containing only variants for the specified region(s). The `intervals` dict
contains only the relevant region keys. A `region_sources` dict is added to
`information`, mapping each region key to its source tool:

```json
{
  "no_snps": 3,
  "HV1": "ATCG...",
  "HV2": "",
  "HV3": "",
  "variants": { "snps": [...], "insertions": [...], "deletions": [...] },
  "intervals": { "HV1": [[16024, 16365]] },
  "sample_flags": ["Autopass HV1"],
  "variant_flags": {},
  "information": {
    "source_tool": "mutation_surveyor",
    "batch_id": "MS_210426_003",
    "region_sources": {
      "HV1": "mutation_surveyor"
    }
  }
}
```

## Warning Flags

Intervals that exceed region boundaries produce warning flags like:

```
Range extends beyond HV1 boundaries: [16020, 16370] vs [16024, 16365]
```

These are warnings, not hard errors — the pipeline continues processing.
