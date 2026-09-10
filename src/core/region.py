"""Per-region intermediate JSON I/O and range/interval QC for mtDNA pipeline.

This module combines two concerns for per-region intermediate JSONs:

1. **Region JSON I/O**: Read, write, and discover per-region intermediate
   JSON files (HV1_{LID}.json, HV2-3_{LID}.json) that serve as temporary
   interchange between tool pipelines and the unify step. These files are
   NOT the final output — the unify step merges them into a single per-LID
   unified JSON in the standard format.

2. **Range/interval QC**: Validate that intervals declared in a per-region
   JSON actually fall within the claimed HV region boundaries. Uses warning
   flags (not hard errors) to maintain consistency with the existing
   flag-based QC approach.

Region file naming convention:
    HV1_{LID}.json    — HV1 region data
    HV2-3_{LID}.json  — HV2+HV3 region data

Files are placed in region group subdirectories:
    {batch_dir}/HV1/HV1_{LID}.json    — HV1 region data
    {batch_dir}/HV2-3/HV2-3_{LID}.json — HV2+HV3 region data
Each region JSON includes `region_sources` provenance tracking in `information`.
"""

import json
from pathlib import Path

from loguru import logger

from src.config import get_settings
from src.core.models import Sample, Tool
from src.core.sample import _normalize_keyed_sample, filter_sample_by_regions, sample_to_dict

# ---------------------------------------------------------------------------
# Region group definitions
# ---------------------------------------------------------------------------

REGION_FILE_PATTERNS: dict[str, str] = {
    "HV1": "HV1_{lid}.json",
    "HV2-3": "HV2-3_{lid}.json",
}

REGION_TO_KEYS: dict[str, list[str]] = {
    "HV1": ["HV1"],
    "HV2-3": ["HV2", "HV3"],
}

AUTOPASS_FLAG_TO_REGION_GROUPS: dict[str, list[str]] = {
    "Autopass HV1": ["HV1"],
    "Autopass HV2 and HV3": ["HV2-3"],
}


# ---------------------------------------------------------------------------
# Region group helpers
# ---------------------------------------------------------------------------


def region_groups_from_flags(sample_flags: list[str]) -> dict[str, list[str]]:
    """Map auto-pass sample flags to region file groups.

    Returns:
        Dict mapping region group names to region key lists.
        E.g. {"HV1": ["HV1"], "HV2-3": ["HV2", "HV3"]}

        Only includes region groups that are auto-pass.

    """
    groups: dict[str, list[str]] = {}
    for flag in sample_flags:
        group_names = AUTOPASS_FLAG_TO_REGION_GROUPS.get(flag)
        if group_names:
            for group_name in group_names:
                if group_name not in groups:
                    groups[group_name] = REGION_TO_KEYS[group_name]
    return groups


# ---------------------------------------------------------------------------
# Per-region JSON I/O
# ---------------------------------------------------------------------------


def write_region_jsons(
    sample: Sample,
    output_dir: Path,
    ref_path: str = "ref/rCRS.fasta",
    region_groups: dict[str, list[str]] | None = None,
) -> dict[str, Path]:
    """Write per-region intermediate JSONs for auto-pass regions.

    For each region group that has intervals and variants, write a separate
    JSON file containing only that region's data.

    Per-region JSONs are temporary intermediate files — the unify step merges
    them into a single per-LID unified JSON in the standard format.

    Region JSONs are organized by region group in subdirectories:
        {output_dir}/HV1/HV1_{LID}.json
        {output_dir}/HV2-3/HV2-3_{LID}.json

    Args:
        sample: Source Sample (already filtered for auto-pass regions).
        output_dir: Directory for region JSON output (typically the batch directory).
        ref_path: Path to rCRS reference FASTA.
        region_groups: Mapping of group name → region keys. If None, derived
            from sample_flags via region_groups_from_flags().

    Returns:
        Dict mapping region group name → output file path.
        Only includes regions that have intervals and variants.

    """
    if region_groups is None:
        region_groups = region_groups_from_flags(sample.sample_flags)

    written: dict[str, Path] = {}
    skipped: dict[str, str] = {}

    for group_name, region_keys in region_groups.items():
        filtered = filter_sample_by_regions(sample, region_keys)
        if filtered is None:
            skipped[group_name] = "filter returned None"
            continue
        if not filtered.variants and not filtered.intervals:
            skipped[group_name] = "no variants or intervals"
            continue

        json_dir = output_dir / group_name
        json_dir.mkdir(parents=True, exist_ok=True)

        filename = REGION_FILE_PATTERNS[group_name].format(lid=sample.sample_id)
        filepath = json_dir / filename

        sample_dict = sample_to_dict(filtered, ref_path)
        # Add region_sources to information for provenance
        if "information" not in sample_dict or sample_dict["information"] is None:
            sample_dict["information"] = {}
        sample_dict["information"]["region_sources"] = dict.fromkeys(region_keys, sample.source_tool.value)

        with filepath.open("w", encoding="utf-8") as fh:
            json.dump(sample_dict, fh, indent=2)

        written[group_name] = filepath

    logger.info(
        "Region JSON write {}: {}/{} groups written ({}), {} skipped ({})",
        sample.sample_id,
        len(written),
        len(region_groups),
        ", ".join(written.keys()) if written else "none",
        len(skipped),
        ", ".join(f"{k}: {v}" for k, v in skipped.items()) if skipped else "none",
    )

    return written


def read_region_json(path: Path, sample_id: str) -> Sample:
    """Read a per-region intermediate JSON file.

    Provenance is carried via the ``region_sources`` dict in ``information``,
    which maps each region key to its source tool (e.g. {"HV1": "mutation_surveyor"}).

    Args:
        path: Path to the per-region JSON file.
        sample_id: Expected sample ID (for validation).

    Returns:
        Sample object with ``information.region_sources`` provenance.

    """
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)

    entry = _normalize_keyed_sample(data)
    entry["sample_id"] = sample_id

    information = data.get("information", {})
    # Remove stale region_group if present in legacy files
    if isinstance(information, dict):
        information.pop("region_group", None)
    source_tool_val = entry.pop("source_tool", None) or information.get("source_tool")
    if isinstance(source_tool_val, str) and source_tool_val:
        try:
            entry["source_tool"] = Tool(source_tool_val)
        except ValueError:
            logger.warning("Unknown source_tool '{}' for {}, marking as unknown", source_tool_val, sample_id)
            entry["source_tool"] = Tool.UNKNOWN
    else:
        # No tool attribution in the region JSON: surface as "Unknown".
        entry["source_tool"] = Tool.UNKNOWN

    batch_id_val = entry.pop("batch_id", None) or information.get("batch_id")
    if batch_id_val:
        entry["batch_id"] = batch_id_val

    return Sample.model_validate(entry)


# ---------------------------------------------------------------------------
# Range/interval QC validation
# ---------------------------------------------------------------------------


def _region_boundaries() -> dict[str, tuple[int, int]]:
    """Return canonical region boundaries from settings."""
    regions = get_settings().regions.REGIONS
    return {name: (bounds[0], bounds[1]) for name, bounds in regions.items()}


def validate_intervals_in_region(
    intervals: dict[str, list[list[int]]],
    region_group: str,
    region_keys: list[str],
) -> list[str]:
    """Validate that all intervals fall within the claimed region boundaries.

    For each interval in the intervals dict, checks that start and end
    positions fall within the canonical boundaries of the declared region.

    Intervals that exceed region boundaries get a warning flag added to
    the returned list, but do NOT raise an error — consistent with the
    existing flag-based QC approach.

    Args:
        intervals: Dict mapping region names to list of [start, end] pairs.
        region_group: The region group name (e.g. "HV1", "HV2-3").
        region_keys: The region keys this group covers (e.g. ["HV1"] or ["HV2", "HV3"]).

    Returns:
        List of warning flag strings for out-of-range intervals.
        Empty list means all intervals are within boundaries.

    """
    boundaries = _region_boundaries()
    warnings: list[str] = []

    for region_name in region_keys:
        if region_name not in intervals:
            continue

        region_bounds = boundaries.get(region_name)
        if region_bounds is None:
            continue

        start_bound, end_bound = region_bounds

        for iv in intervals[region_name]:
            iv_start, iv_end = iv[0], iv[1]

            # Allow slight extension beyond canonical boundaries (common in sequencing)
            # but flag significant excursions
            if iv_start < start_bound or iv_end > end_bound:
                flag = (
                    f"Range extends beyond {region_name} boundaries: "
                    f"[{iv_start}, {iv_end}] vs [{start_bound}, {end_bound}]"
                )
                warnings.append(flag)
                logger.warning(
                    "Region QC: {} — interval [{}, {}] extends beyond {} [{}, {}]",
                    region_group,
                    iv_start,
                    iv_end,
                    region_name,
                    start_bound,
                    end_bound,
                )

    return warnings


def validate_region_group_intervals(
    intervals: dict[str, list[list[int]]],
    region_group: str,
) -> list[str]:
    """Validate intervals for a region group using the canonical mapping.

    Convenience wrapper that looks up region keys from the group name.

    Args:
        intervals: Dict mapping region names to list of [start, end] pairs.
        region_group: The region group name (e.g. "HV1", "HV2-3").

    Returns:
        List of warning flag strings.

    """
    region_keys = REGION_TO_KEYS.get(region_group, [])
    if not region_keys:
        logger.warning("Unknown region group: {}", region_group)
        return [f"Unknown region group: {region_group}"]

    return validate_intervals_in_region(intervals, region_group, region_keys)
