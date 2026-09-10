"""Sequencher ETL module — TXT → Sample.

This module provides the ETL (Extract, Transform, Load) pipeline for Sequencher
TXT output files. It parses raw Sequencher data, validates analysis ranges,
and produces Sample objects with full source provenance.

ETL is pure: it populates only variants, intervals, and flags. HV sequences
and statistics are computed by Batch at aggregation time.

The entry point is `process()` which takes a single TXT file and returns a
Sample. For batch processing, use ``src.tools.sequencher.pipeline.process_batch``.
"""

from pathlib import Path
from typing import Any

from loguru import logger

from src.config import get_settings
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Sample, Tool, Variant
from src.core.variants import (
    normalize_position,
)
from src.tools.sequencher.quality_control import (
    _expected_regions_for_prefix,
    parse_filename_ranges,
    validate_analysis_ranges,
)
from src.tools.sequencher.utils import get_intervals_from_manual

_MIN_VARIANT_PARTS = 3


def _default_seq_regions() -> list[list[int]]:
    """Return the full HV1/HV2/HV3 region list from settings."""
    regions = get_settings().regions.REGIONS
    return [list(regions["HV1"]), list(regions["HV2"]), list(regions["HV3"])]


def _parse_sequencher_txt(file_path: str) -> dict[str, dict[str, Any]]:
    """Parse Sequencher TXT file to extract variants.

    Expected format::

        Pos         Seq             Con           Required Edit
        73          A               G            Change base
        309.1       :               C            Insert base
        16,093      T               -            Delete base

    Args:
        file_path: Path to the TXT file.

    Returns:
        Dictionary mapping normalized position to variant info dict.
        Each variant dict has keys: pos, ref, seq, file, quality.
    """
    variants: dict[str, dict[str, Any]] = {}

    with Path(file_path).open() as f:
        lines = f.readlines()

    # Find the data section (after the header)
    data_start = 0
    for i, line in enumerate(lines):
        if "Pos" in line and "Seq" in line and "Con" in line:
            data_start = i + 2  # Skip header and blank line
            break

    # Parse each variant line
    for raw_line in lines[data_start:]:
        stripped = raw_line.strip()
        if not stripped:
            continue

        parts = stripped.split()
        if len(parts) < _MIN_VARIANT_PARTS:
            continue

        pos_str = parts[0].replace(",", "")
        pos = normalize_position(pos_str)
        if pos is None:
            continue

        # Parse reference and alternate bases
        ref_base = parts[1] if parts[1] != ":" else "-"
        alt_base = parts[2] if parts[2] != ":" else "-"

        variants[str(pos)] = {
            "pos": normalize_position(pos),
            "ref": ref_base,
            "seq": alt_base,
            "file": [],
            "quality": [],
        }

    return variants


def _resolve_seq_regions(
    ranges: list[list[int]] | str,
    filename: str,
    region_prefix: str | None = None,
) -> tuple[list[list[int]], list[str]]:
    """Resolve sequencing regions from parsed filename ranges.

    Args:
        ranges: Parsed ranges from filename (list of int pairs, or "FULL REGION").
        filename: Original filename (used in warning messages).

    Returns:
        Tuple of (seq_regions, range_flags).
    """
    range_flags: list[str] = []

    if ranges == "FULL REGION" or (isinstance(ranges, list) and len(ranges) == 0):
        return _default_seq_regions(), range_flags

    if not isinstance(ranges, list):
        return _default_seq_regions(), range_flags

    should_display, flag = validate_analysis_ranges(ranges, region_prefix=region_prefix, filename=filename)
    if not should_display:
        logger.warning("File {}: invalid ranges {}", filename, ranges)
        if flag:
            range_flags.append(flag)
        return _default_seq_regions(), range_flags

    if flag:
        range_flags.append(flag)

    return ranges, range_flags


def process(
    sample_id: str,
    input_path: Path,
    batch_id: str | None = None,
) -> Sample:
    """Main ETL: Sequencher TXT → Sample.

    ETL is pure: returns Sample with variants, intervals, and flags only.
    HV sequences (hv1/hv2/hv3) and statistics (no_snps etc.) are computed
    by Batch at aggregation time.

    Args:
        sample_id: Sample identifier (LID from filename).
        input_path: Path to Sequencher TXT file.
        batch_id: Optional batch identifier.

    Returns:
        Sample with parsed variants, intervals, source provenance, and flags.
    """
    input_path = Path(input_path)

    # Parse filename for ranges
    parse_result = parse_filename_ranges(input_path.name)
    seq_regions, range_flags = _resolve_seq_regions(
        parse_result.ranges, input_path.name, region_prefix=parse_result.region_prefix
    )
    if parse_result.has_format_error and parse_result.format_error and parse_result.format_error not in range_flags:
        range_flags.append(parse_result.format_error)

    # Parse TXT file
    variants_dict = _parse_sequencher_txt(str(input_path))

    # Build intervals dict
    intervals = get_intervals_from_manual(seq_regions)
    intervals_dict = {}
    for region, spans in intervals.items():
        intervals_dict[region] = [list(span) for span in spans]

    # Build variant list
    variant_list: list[Variant] = []
    for var_data in variants_dict.values():
        # Handle quality field
        quality_raw = var_data.get("quality", [])
        quality_list: list[int] = []
        if quality_raw:
            if isinstance(quality_raw, list):
                quality_list = [int(q) for q in quality_raw]
            elif isinstance(quality_raw, (int, float)):
                quality_list = [int(quality_raw)]

        variant = Variant(
            pos=var_data["pos"],
            ref=var_data["ref"],
            seq=var_data["seq"],
            quality=quality_list,
            files=var_data.get("file", []),
            peaks=var_data.get("peaks"),
        )
        variant_list.append(variant)

    # Keep ALL parsed variants — including positions outside the analyzed
    # intervals — so the comparison layer sees the full picture and decides
    # which variants matter for concordance. The Sample still carries
    # ``intervals`` so downstream range checks remain available.

    # Compute flags: sample-level (range/QC + analyze) and per-variant
    sample_flags_list = list(range_flags or [])
    variant_flags_dict = flag_variants(variant_list)

    # Also run sample-level flag analysis for flags not captured per-variant
    # (e.g. "No 315.1 variant", "Consecutive indels (5+)")
    sequenced_regions = _expected_regions_for_prefix(parse_result.region_prefix or "")
    flagging_intervals = {region: spans for region, spans in intervals_dict.items() if region in sequenced_regions}
    flagger = SampleFlagger(variant_list, intervals=flagging_intervals or None)
    _, analyze_reasons, _ = flagger.analyze()
    for reason in analyze_reasons:
        if reason not in sample_flags_list:
            sample_flags_list.append(reason)

    # Per-variant flags belong only in variant_flags; drop any sample-level
    # flag that duplicates a per-variant flag string.
    sample_flags_list = deduplicate_sample_flags(sample_flags_list, variant_flags_dict)

    # Build Sample — ETL is pure: only variants, intervals, and flags.
    # HV sequences and statistics are computed by Batch.
    return Sample(
        sample_id=sample_id,
        variants=variant_list,
        source_tool=Tool.SEQUENCHER,
        intervals=intervals_dict or None,
        batch_id=batch_id,
        sample_flags=sample_flags_list or [],
        variant_flags=variant_flags_dict,
        information=None,
        hv1=None,
        hv2=None,
        hv3=None,
        no_snps=None,
        no_ins=None,
        no_dels=None,
        no_identicals=None,
        no_unread=None,
    )
