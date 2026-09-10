"""Region and interval utilities for Sequencher.

Filename parsing and analysis-range validation (plus the pre-ETL QC gate) live
in ``quality_control.py``. This module keeps the region/interval helpers and
region-group detection used by the ETL and pipeline layers.
"""

from loguru import logger

from src.config import get_settings
from src.core.variants import (
    get_overlap,
)
from src.tools.sequencher.quality_control import parse_filename_ranges

NUM_REGIONS = 3

# ---------------------------------------------------------------------------
# Range / interval utilities
# ---------------------------------------------------------------------------


def parse_intervals_string(intervals_string: str) -> list[list[int]]:
    """
    Parse intervals string from either format into list of [start, end] pairs.

    Handles:
    - Old format: "73-340, 438-573, 16024-16365" (comma-separated)
    - New format: "73-302 316-340 438-576 16024-16365" (space-separated)

    Args:
        intervals_string: String containing interval ranges

    Returns:
        List of [start, end] pairs

    Examples:
        parse_intervals_string("73-340, 438-573, 16024-16365")
        -> [[73, 340], [438, 573], [16024, 16365]]

        parse_intervals_string("73-302 316-340 438-576 16024-16365")
        -> [[73, 302], [316, 340], [438, 576], [16024, 16365]]
    """
    if not intervals_string or intervals_string == "FULL REGION":
        return []

    if "," in intervals_string:
        # Old format: comma-separated
        logger.debug(f"Parsing old format intervals: {intervals_string}")
        interval_parts = [part.strip() for part in intervals_string.split(",")]
    else:
        # New format: space-separated
        logger.debug(f"Parsing new format intervals: {intervals_string}")
        interval_parts = intervals_string.split()

    intervals = []
    for part in interval_parts:
        if "-" in part:
            try:
                start, end = part.split("-")
                intervals.append([int(start), int(end)])
            except ValueError as e:
                logger.warning(f"Could not parse interval '{part}': {e}")
                continue

    return intervals


def get_intervals_from_manual(seq_regions: list[list[int]]) -> dict:
    """Get genomic region intervals from manual sequencing regions.

    Maps each HV region to its overlap with the provided sequencing regions.

    Args:
        seq_regions: List of [start, end] intervals from manual sequencing.

    Returns:
        Dictionary mapping region names to their overlapping intervals.
    """
    regions = get_settings().regions.REGIONS
    intervals = {}
    for region in regions:
        intervals[region] = get_overlap([regions[region]], seq_regions)
    return intervals


def is_full_region(intervals: list[list[int]]) -> str:
    """Check if intervals match predefined HV regions.

    Args:
        intervals: List of [start, end] pairs.

    Returns:
        ``"FULL REGION"`` if all 3 HV regions are present,
        space-separated range string otherwise, or ``"None"`` for empty input.
    """
    regions = get_settings().regions.REGIONS

    if not intervals:
        return "None"

    # Check if we have all three predefined regions
    if len(intervals) == NUM_REGIONS:
        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        predefined_ranges = sorted(regions.values(), key=lambda x: x[0])

        match_all = True
        for i, (start, end) in enumerate(sorted_intervals):
            pred_start, pred_end = predefined_ranges[i]
            if start != pred_start or end != pred_end:
                match_all = False
                break

        if match_all:
            return "FULL REGION"

    # Format intervals as space-separated string (no commas)
    interval_strs = [f"{start}-{end}" for start, end in intervals]
    return " ".join(interval_strs)


# ---------------------------------------------------------------------------
# Region group detection from filename
# ---------------------------------------------------------------------------


def _detect_region_group_from_ranges(ranges: list[list[int]]) -> str | None:
    """Analyze range boundaries to determine which HV region group they cover.

    Args:
        ranges: List of [start, end] pairs from the filename.

    Returns:
        "HV1" if only HV1, "HV2-3" if only HV2/HV3, None if full region or indeterminate.
    """
    regions = get_settings().regions.REGIONS
    hv1_start, hv1_end = regions["HV1"]
    hv2_start, hv2_end = regions["HV2"]
    hv3_start, hv3_end = regions["HV3"]

    has_hv1 = any(hv1_start <= s <= hv1_end or hv1_start <= e <= hv1_end for s, e in ranges)
    has_hv2 = any(hv2_start <= s <= hv2_end or hv2_start <= e <= hv2_end for s, e in ranges)
    has_hv3 = any(hv3_start <= s <= hv3_end or hv3_start <= e <= hv3_end for s, e in ranges)

    if has_hv1 and not has_hv2 and not has_hv3:
        return "HV1"
    if (has_hv2 or has_hv3) and not has_hv1:
        return "HV2-3"
    return None


def detect_region_group_from_filename(filename: str) -> str | None:
    """Detect the region group from a Sequencher TXT filename.

    Sequencher TXT filenames follow the pattern:
        LID.TXT                    → full region (all 3 HV regions)
        LID-a1-b1 a2-b2 ... .TXT  → ranges in filename
        HV1_LID.TXT                → HV1 region
        HV2-3_LID.TXT              → HV2-3 region
        HV1_LID-ranges.TXT         → HV1 region (with explicit ranges)
        HV2-3_LID-ranges.TXT       → HV2-3 region (with explicit ranges)

    The region group is determined by:
    1. The region prefix (HV1_ or HV2-3_) if present in the filename.
    2. Otherwise, by analyzing the ranges against canonical HV region boundaries.

    Args:
        filename: The filename to analyze (e.g. "2513827-16024-16365.TXT").

    Returns:
        "HV1", "HV2-3", or None for full-region files.
    """
    parse_result = parse_filename_ranges(filename)

    if parse_result.region_prefix is not None:
        return parse_result.region_prefix

    if parse_result.is_full_region or parse_result.has_format_error:
        return None

    ranges = parse_result.ranges
    if not ranges or not isinstance(ranges, list):
        return None

    return _detect_region_group_from_ranges(ranges)
