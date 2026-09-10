"""Tool-agnostic pairwise comparison of Sample objects.

Produces ComparisonResult with binary concordance and domain-specific flagging.
Reads pre-computed ``sample_flags`` and ``variant_flags`` from each Sample directly.
"""

import csv
import json
import sys
from argparse import ArgumentParser
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from src.config import get_settings
from src.core.models import (
    ComparisonResult,
    Sample,
    Tool,
    ToolResult,
    Variant,
)
from src.core.sample import load_sample_batch
from src.core.variants import (
    format_variants_simplified,
    get_hv_region_for_position,
    is_position_in_intervals,
    is_special_position,
    normalize_position,
    pos_base,
    pos_sort_key,
    variant_to_dict,
)

# Positions always treated as concordant regardless of content.
ALWAYS_CONCORDANT_POSITIONS = ["573.1", "309.1", "309.2"]

RANGE_FLAGS = {"Range - Possibly wrong", "Range - Missing"}
MANUAL_TOOLS = {Tool.SEQUENCHER, Tool.MUTATION_SURVEYOR}

# Display name used when a batch carries no tool attribution (e.g. empty batch
# or samples without a source_tool). The comparison output column headers are
# built dynamically from the source tool of each batch instead of the static
# "Tool A" / "Tool B" convention.
UNKNOWN_TOOL_NAME = "Unknown"


def _tool_display_name(tool: Tool | None) -> str:
    """Return the column display name for a tool, or "Unknown" if absent.

    The ``Tool.UNKNOWN`` sentinel (a batch/sample with no tool attribution)
    also displays as "Unknown".
    """
    if tool is None or tool is Tool.UNKNOWN:
        return UNKNOWN_TOOL_NAME
    return tool.value


def _batch_tool_name(batch: Sequence[Sample]) -> str:
    """Derive a single tool display name for a whole batch.

    Each batch corresponds to one analysis tool, so the first sample's
    ``source_tool`` is used. An empty batch (no samples to read the tool
    from) yields "Unknown".
    """
    if not batch:
        return UNKNOWN_TOOL_NAME
    return _tool_display_name(batch[0].source_tool)


def _summary_fields(tool_name_a: str, tool_name_b: str) -> list[str]:
    """Column order for the comparison summary TSV/Excel, parameterized by tool."""
    return [
        "Sample ID",
        "Batch",
        f"Analyzed Range ({tool_name_a})",
        f"Variants ({tool_name_a})",
        f"Variants Unique ({tool_name_a})",
        f"Variants Count ({tool_name_a})",
        f"Analyzed Range ({tool_name_b})",
        f"Variants ({tool_name_b})",
        f"Variants Unique ({tool_name_b})",
        f"Variants Count ({tool_name_b})",
        "Concordant",
        "Flag",
        f"Sample Flags ({tool_name_a})",
        f"Sample Flags ({tool_name_b})",
        f"Variant Flags ({tool_name_a})",
        f"Variant Flags ({tool_name_b})",
    ]


def _variant_fields(tool_name_a: str, tool_name_b: str) -> list[str]:
    """Column order for the variant-level TSV, parameterized by tool."""
    return [
        "Sample ID",
        "Batch",
        "Region",
        "Position",
        f"{tool_name_a} Variant",
        f"Variant Flags ({tool_name_a})",
        f"{tool_name_b} Variant",
        f"Variant Flags ({tool_name_b})",
        f"Sample Flags ({tool_name_a})",
        f"Sample Flags ({tool_name_b})",
    ]


def _tool_range_error(tool_result: ToolResult) -> bool:
    """Range error display only applies to manual tools."""
    if tool_result.tool not in MANUAL_TOOLS:
        return False
    return any(r in RANGE_FLAGS for r in tool_result.flag_reasons)


def _format_tool_columns(
    tool_result: ToolResult,
    unique_variants: list[dict[str, Any]],
    *,
    range_error: bool,
) -> dict[str, Any]:
    if range_error:
        return {
            "range": "Error",
            "variants": "Error",
            "unique": "Error",
            "count": 0,
        }
    return {
        "range": tool_result.analyzed_range,
        "variants": format_variants_simplified(tool_result.variants),
        "unique": format_variants_simplified(unique_variants) if unique_variants else "",
        "count": len(tool_result.variants),
    }


# Legacy comparison summary layout (Batch_MS_*.tsv): hardcoded tool labels,
# a single Variants Count column (Sequencher/B only), and merged per-tool
# "Flagged (<tool>)" columns combining sample + variant flags.


def _legacy_summary_fields() -> list[str]:
    """Column order for the legacy 13-column comparison summary TSV/Excel."""
    return [
        "Sample ID",
        "Batch",
        "Analyzed Range (Pipeline)",
        "Variants (Pipeline)",
        "Variants Unique (Pipeline)",
        "Analyzed Range (Sequencher)",
        "Variants (Sequencher)",
        "Variants Unique (Sequencher)",
        "Variants Count (Sequencher)",
        "Concordant",
        "Flag",
        "Flagged (Pipeline)",
        "Flagged (Sequencher)",
    ]


def _format_legacy_flagged(
    flag_reasons: list[str],
    variant_flags: dict[str, list[str]],
) -> str:
    """Merge sample + variant flags into one legacy "Flagged" column.

    Deduplicated union joined by ", "; "No" when empty. Variant-style flags
    come first (per-variant reasons, then the consolidated "16180-16193
    region (...)" entry), matching the legacy file ordering; sample-style flags
    (No 315.1 variant, Has-... variant, Consecutive indels, Lowercase, Range)
    come last. After ETL deduplication the two sets do not overlap. Insertion
    flags keep their current "Insertion {seq} at {pos}" form.
    """
    parts: list[str] = []
    seen: set[str] = set()
    variant_str = _format_variant_flags(variant_flags)
    if variant_str:
        for reason in variant_str.split("; "):
            if reason and reason not in seen:
                parts.append(reason)
                seen.add(reason)
    for reason in flag_reasons:
        if reason and reason not in seen:
            parts.append(reason)
            seen.add(reason)
    return ", ".join(parts) if parts else "No"


def _legacy_summary_row(result: ComparisonResult, batch_id: str = "") -> dict[str, Any]:
    """Build one legacy 13-column summary row from a ComparisonResult.

    Hardcodes Pipeline (A) / Sequencher (B) labels, emits a single Variants
    Count (Sequencher) column, and merges each tool's flags into a
    "Flagged (<tool>)" column. Range-error "Error" display for manual tools
    is preserved via the shared _format_tool_columns helper.
    """
    a_cols = _format_tool_columns(
        result.tool_a,
        result.unique_variants_a,
        range_error=_tool_range_error(result.tool_a),
    )
    b_cols = _format_tool_columns(
        result.tool_b,
        result.unique_variants_b,
        range_error=_tool_range_error(result.tool_b),
    )
    overall_flag = "Yes" if (result.tool_a.flagged or result.tool_b.flagged) else "No"
    return {
        "Sample ID": result.sample_id,
        "Batch": batch_id,
        "Analyzed Range (Pipeline)": a_cols["range"],
        "Variants (Pipeline)": a_cols["variants"],
        "Variants Unique (Pipeline)": a_cols["unique"],
        "Analyzed Range (Sequencher)": b_cols["range"],
        "Variants (Sequencher)": b_cols["variants"],
        "Variants Unique (Sequencher)": b_cols["unique"],
        "Variants Count (Sequencher)": b_cols["count"],
        "Concordant": result.concordance,
        "Flag": overall_flag,
        "Flagged (Pipeline)": _format_legacy_flagged(result.tool_a.flag_reasons, result.tool_a.variant_flags),
        "Flagged (Sequencher)": _format_legacy_flagged(result.tool_b.flag_reasons, result.tool_b.variant_flags),
    }


class SampleComparator:
    """Pairwise comparison of two Sample objects.

    Stateless — all methods are pure functions over their inputs.
    """

    @staticmethod
    def format_intervals(intervals: dict[str, list[list[int]]] | None) -> str:
        if not intervals:
            return "None"
        all_ranges = [rng for ranges in intervals.values() for rng in ranges]
        if not all_ranges:
            return "None"
        # Check if intervals cover all predefined regions (FULL REGION)
        regions = get_settings().regions.REGIONS
        if len(intervals) == len(regions):
            predefined_ranges = sorted(regions.values(), key=lambda x: x[0])
            interval_list = sorted(all_ranges, key=lambda r: r[0])
            if len(interval_list) == len(predefined_ranges) and all(
                interval_list[i][0] == predefined_ranges[i][0] and interval_list[i][1] == predefined_ranges[i][1]
                for i in range(len(predefined_ranges))
            ):
                return "FULL REGION"
        sorted_ranges = sorted(all_ranges, key=lambda r: r[0])
        return " ".join(f"{s}-{e}" for s, e in sorted_ranges)

    @staticmethod
    def find_unique_variants(
        variants_a: Sequence[Variant],
        variants_b: Sequence[Variant],
        *,
        ignore_special: bool = False,
    ) -> list[Variant]:
        """Variants in *variants_a* not present in *variants_b*."""
        if not variants_a:
            return []
        if not variants_b:
            if ignore_special:
                return [v for v in variants_a if not is_special_position(v.pos, v.ref, v.seq)]
            return list(variants_a)

        b_keys = {(str(normalize_position(v.pos)), v.ref, v.seq) for v in variants_b}
        result = []
        for v in variants_a:
            if ignore_special and is_special_position(v.pos, v.ref, v.seq):
                continue
            if (str(normalize_position(v.pos)), v.ref, v.seq) not in b_keys:
                result.append(v)
        return result

    @staticmethod
    def _non_special_variants(sample: Sample) -> list[Variant]:
        """Filter out special positions from a sample's variants."""
        return [v for v in sample.variants if not is_special_position(v.pos, v.ref, v.seq)]

    @staticmethod
    def _collect_intervals(sample: Sample) -> list[tuple[int, int]]:
        """Collect analyzed intervals from a sample, warning if absent."""
        if not sample.intervals:
            logger.warning("Sample {} has no intervals; concordance check may be unreliable", sample.sample_id)
            return []
        return [(span[0], span[1]) for spans in sample.intervals.values() for span in spans]

    @staticmethod
    def _shared_positions_differ(
        variants_a: list[Variant],
        variants_b: list[Variant],
    ) -> bool:
        """Check if any shared position has differing ref/seq between A and B."""
        a_dict = {str(normalize_position(v.pos)): v for v in variants_a}
        b_dict = {str(normalize_position(v.pos)): v for v in variants_b}
        for pos in set(a_dict) & set(b_dict):
            if pos in ALWAYS_CONCORDANT_POSITIONS:
                continue
            if (a_dict[pos].ref, a_dict[pos].seq) != (b_dict[pos].ref, b_dict[pos].seq):
                return True
        return False

    @staticmethod
    def is_concordant(result_a: Sample, result_b: Sample) -> str:
        """Binary concordance: 'Y', 'N', or 'N/A' (empty variants).

        Checks concordance between two Sample results. B-only variants
        that fall outside A's analyzed ranges are ignored —
        Sequencher-only variants outside the pipeline's ranges are excluded
        from concordance.
        """
        if not result_a.variants or not result_b.variants:
            return "N/A"

        non_special_a = SampleComparator._non_special_variants(result_a)
        non_special_b = SampleComparator._non_special_variants(result_b)
        a_intervals = SampleComparator._collect_intervals(result_a)

        if SampleComparator._shared_positions_differ(non_special_a, non_special_b):
            return "N"

        unique_a = SampleComparator.find_unique_variants(non_special_a, non_special_b, ignore_special=True)
        unique_b = SampleComparator.find_unique_variants(non_special_b, non_special_a, ignore_special=True)

        # A-only variants in HV regions → discordant
        if any(get_hv_region_for_position(v.pos) is not None for v in unique_a):
            return "N"

        # B-only variants: only discordant if within A's analyzed ranges
        for v in unique_b:
            if get_hv_region_for_position(v.pos) is None:
                continue
            if not a_intervals or is_position_in_intervals(v.pos, a_intervals):
                return "N"

        return "Y"

    @staticmethod
    def compare(result_a: Sample, result_b: Sample) -> ComparisonResult:
        """Full pairwise comparison returning a ComparisonResult."""
        non_special_a = [v for v in result_a.variants if not is_special_position(v.pos, v.ref, v.seq)]
        non_special_b = [v for v in result_b.variants if not is_special_position(v.pos, v.ref, v.seq)]

        unique_a = SampleComparator.find_unique_variants(non_special_a, non_special_b, ignore_special=True)
        unique_b = SampleComparator.find_unique_variants(non_special_b, non_special_a, ignore_special=True)

        concordance = SampleComparator.is_concordant(result_a, result_b)

        return ComparisonResult(
            sample_id=result_a.sample_id,
            tool_a=ToolResult(
                tool=result_a.source_tool,
                analyzed_range=SampleComparator.format_intervals(result_a.intervals),
                variants=[d for d in (variant_to_dict(v) for v in result_a.variants) if d is not None],
                flagged=bool(result_a.sample_flags) or bool(result_a.variant_flags),
                flag_reasons=result_a.sample_flags,
                variant_flags=result_a.variant_flags,
            ),
            tool_b=ToolResult(
                tool=result_b.source_tool,
                analyzed_range=SampleComparator.format_intervals(result_b.intervals),
                variants=[d for d in (variant_to_dict(v) for v in result_b.variants) if d is not None],
                flagged=bool(result_b.sample_flags) or bool(result_b.variant_flags),
                flag_reasons=result_b.sample_flags,
                variant_flags=result_b.variant_flags,
            ),
            unique_variants_a=[d for d in (variant_to_dict(v) for v in unique_a) if d is not None],
            unique_variants_b=[d for d in (variant_to_dict(v) for v in unique_b) if d is not None],
            concordance=concordance,
        )

    @staticmethod
    def compare_batches(
        batch_a: Sequence[Sample],
        batch_b: Sequence[Sample],
        batch_id: str = "",
        *,
        legacy_format: bool = False,
    ) -> dict[str, Any]:
        a_by_id = {s.sample_id: s for s in batch_a}
        b_by_id = {s.sample_id: s for s in batch_b}

        common_ids = sorted(set(a_by_id) & set(b_by_id))
        a_only = sorted(set(a_by_id) - set(b_by_id))
        b_only = sorted(set(b_by_id) - set(a_by_id))

        comparison_models = [SampleComparator.compare(a_by_id[sid], b_by_id[sid]) for sid in common_ids]

        # Column headers are parameterized by the source tool of each batch so
        # the output names the actual tools instead of the generic "A"/"B".
        # A batch without samples (or without tool attribution) falls back to
        # "Unknown".
        tool_name_a = _batch_tool_name(batch_a)
        tool_name_b = _batch_tool_name(batch_b)

        variant_rows: list[dict[str, Any]] = []
        for result in comparison_models:
            variant_rows.extend(
                create_variant_level_data(
                    result,
                    batch_id=batch_id,
                    tool_name_a=tool_name_a,
                    tool_name_b=tool_name_b,
                )
            )

        if legacy_format:
            summary_rows = [_legacy_summary_row(r, batch_id=batch_id) for r in comparison_models]
            summary_fields = _legacy_summary_fields()
        else:
            summary_rows = [
                SampleComparator.summary_row(
                    r,
                    batch_id=batch_id,
                    tool_name_a=tool_name_a,
                    tool_name_b=tool_name_b,
                )
                for r in comparison_models
            ]
            summary_fields = _summary_fields(tool_name_a, tool_name_b)

        concordant = sum(1 for r in comparison_models if r.concordance == "Y")
        discordant = sum(1 for r in comparison_models if r.concordance == "N")
        flagged = sum(1 for r in comparison_models if r.tool_a.flagged or r.tool_b.flagged)

        return {
            "summary": {
                "matched": len(common_ids),
                "concordant": concordant,
                "discordant": discordant,
                "flagged": flagged,
                "a_only": a_only,
                "b_only": b_only,
            },
            "summary_rows": summary_rows,
            "variant_rows": variant_rows,
            "summary_fields": summary_fields,
            "variant_fields": _variant_fields(tool_name_a, tool_name_b),
            "tool_name_a": tool_name_a,
            "tool_name_b": tool_name_b,
            "legacy_format": legacy_format,
        }

    @staticmethod
    def summary_row(
        result: ComparisonResult,
        batch_id: str = "",
        tool_name_a: str = UNKNOWN_TOOL_NAME,
        tool_name_b: str = UNKNOWN_TOOL_NAME,
    ) -> dict[str, Any]:
        """Build one TSV summary row from a ComparisonResult.

        When a manual tool (Sequencher or Mutation Surveyor) has a range
        range flag ("Range - Possibly wrong" or "Range - Missing"),
        that tool's Analyzed Range, Variants, Variants Unique, and Variants
        Count columns show "Error" instead of their values. Automated tools
        (pipeline) never get the "Error" display because their ranges are
        always considered authoritative.
        """
        overall_flag = "Yes" if (result.tool_a.flagged or result.tool_b.flagged) else "No"

        a_cols = _format_tool_columns(
            result.tool_a,
            result.unique_variants_a,
            range_error=_tool_range_error(result.tool_a),
        )
        b_cols = _format_tool_columns(
            result.tool_b,
            result.unique_variants_b,
            range_error=_tool_range_error(result.tool_b),
        )

        return {
            "Sample ID": result.sample_id,
            "Batch": batch_id,
            f"Analyzed Range ({tool_name_a})": a_cols["range"],
            f"Variants ({tool_name_a})": a_cols["variants"],
            f"Variants Unique ({tool_name_a})": a_cols["unique"],
            f"Variants Count ({tool_name_a})": a_cols["count"],
            f"Analyzed Range ({tool_name_b})": b_cols["range"],
            f"Variants ({tool_name_b})": b_cols["variants"],
            f"Variants Unique ({tool_name_b})": b_cols["unique"],
            f"Variants Count ({tool_name_b})": b_cols["count"],
            "Concordant": result.concordance,
            "Flag": overall_flag,
            f"Sample Flags ({tool_name_a})": _format_flag_reasons(
                flagged=bool(result.tool_a.flag_reasons), reasons=result.tool_a.flag_reasons
            ),
            f"Sample Flags ({tool_name_b})": _format_flag_reasons(
                flagged=bool(result.tool_b.flag_reasons), reasons=result.tool_b.flag_reasons
            ),
            f"Variant Flags ({tool_name_a})": _format_variant_flags(result.tool_a.variant_flags),
            f"Variant Flags ({tool_name_b})": _format_variant_flags(result.tool_b.variant_flags),
        }


def _format_variant_display(variant: dict[str, Any] | None) -> str:
    """Format a variant dict for display in variant-level TSV.

    Matches legacy convention: insertion shows seq, deletion shows DEL, SNP shows alt.
    """
    if not variant:
        return "None"
    ref = variant.get("ref", "-")
    alt = variant.get("seq", variant.get("alt", "-"))
    if ref == "-":
        return alt
    if alt == "-":
        return "DEL"
    return alt


def _get_position_flags(
    position: str,
    _variant: dict[str, Any] | None,
    variant_flags: dict[str, list[str]],
) -> str:
    """Look up variant-level flags for a position from the variant_flags dict.

    Searches for any key in variant_flags that matches the position,
    since keys are "pos|ref|alt" format and the position may appear
    with different ref/alt combinations.
    """
    if not variant_flags:
        return ""
    matching = []
    for key, reasons in variant_flags.items():
        key_pos = key.split("|")[0]
        if key_pos == str(position):
            matching.extend(reasons)
    return "; ".join(matching) if matching else ""


def create_variant_level_data(
    result: ComparisonResult,
    batch_id: str = "",
    tool_name_a: str = UNKNOWN_TOOL_NAME,
    tool_name_b: str = UNKNOWN_TOOL_NAME,
) -> list[dict[str, Any]]:
    """Create per-position variant rows from a ComparisonResult.

    Each row shows one position with both A and B variant data,
    plus any associated variant-level flags. Uses variant_flags dict
    keyed by "pos|ref|alt" for per-position flag lookup. Column keys are
    parameterized by the source tool names so the output names the actual
    tools instead of the generic "A"/"B".
    """
    a_variants = result.tool_a.variants
    b_variants = result.tool_b.variants

    a_positions = {str(normalize_position(v.get("pos", v.get("position", "")))) for v in a_variants}
    b_positions = {str(normalize_position(v.get("pos", v.get("position", "")))) for v in b_variants}
    all_positions = a_positions | b_positions

    a_dict = {str(normalize_position(v.get("pos", v.get("position", "")))): v for v in a_variants}
    b_dict = {str(normalize_position(v.get("pos", v.get("position", "")))): v for v in b_variants}

    a_variant_flags = result.tool_a.variant_flags
    b_variant_flags = result.tool_b.variant_flags

    a_sample_flags = "; ".join(result.tool_a.flag_reasons) if result.tool_a.flag_reasons else ""
    b_sample_flags = "; ".join(result.tool_b.flag_reasons) if result.tool_b.flag_reasons else ""

    rows: list[dict[str, Any]] = []
    for position in sorted(all_positions, key=pos_sort_key):
        pos_int = pos_base(position)
        region = get_hv_region_for_position(pos_int) or "Other"

        a_variant = a_dict.get(position)
        b_variant = b_dict.get(position)

        # Look up variant-level flags for this position
        a_flags = _get_position_flags(position, a_variant, a_variant_flags)
        b_flags = _get_position_flags(position, b_variant, b_variant_flags)

        rows.append(
            {
                "Sample ID": result.sample_id,
                "Batch": batch_id,
                "Region": region,
                "Position": position,
                f"{tool_name_a} Variant": _format_variant_display(a_variant),
                f"Variant Flags ({tool_name_a})": a_flags,
                f"{tool_name_b} Variant": _format_variant_display(b_variant),
                f"Variant Flags ({tool_name_b})": b_flags,
                f"Sample Flags ({tool_name_a})": a_sample_flags,
                f"Sample Flags ({tool_name_b})": b_sample_flags,
            },
        )
    return rows


def _format_variant_flags(variant_flags: dict[str, list[str]]) -> str:
    """Format variant_flags dict into a compact string for display.

    Consolidates "16180-16193 region" entries across positions into a single
    "16180-16193 region (pos1, pos2, ...)" entry, matching legacy format.
    Other per-variant flags (Deletion at X, Insertion {seq} at X, etc.) stay per-position.

    Each entry "pos|ref|alt": [reasons...] becomes "reason1; reason2".
    The position prefix is omitted because each reason already embeds its
    position (e.g. "Has 521A-G variant", "Heteroplasmy at 16129",
    "Insertion G at 309.2"). Multiple entries are separated by "; ".
    """
    if not variant_flags:
        return ""

    region_positions: list[str] = []
    non_region_parts: list[str] = []

    for key, reasons in variant_flags.items():
        pos = key.split("|")[0]
        for reason in reasons:
            if reason == "16180-16193 region":
                region_positions.append(pos)
            else:
                non_region_parts.append(reason)

    parts = list(non_region_parts)
    if region_positions:
        parts.append(f"16180-16193 region ({', '.join(sorted(set(region_positions)))})")

    return "; ".join(parts)


def _format_flag_reasons(*, flagged: bool, reasons: list[str]) -> str:
    if not flagged:
        return "No"
    return "; ".join(reasons) if reasons else "Yes"


def write_comparison_json(batch_result: dict[str, Any], output_path: str | Path) -> Path:
    json_path = Path(output_path)
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(batch_result, handle, indent=4, ensure_ascii=False)
    return json_path


def write_comparison_tsv(batch_result: dict[str, Any], output_path: str | Path) -> Path:
    tsv_path = Path(output_path)
    rows = batch_result.get("summary_rows", [])
    fieldnames = batch_result.get("summary_fields") or _summary_fields(UNKNOWN_TOOL_NAME, UNKNOWN_TOOL_NAME)
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return tsv_path


def write_variant_level_tsv(batch_result: dict[str, Any], output_path: str | Path) -> Path:
    tsv_path = Path(output_path)
    rows = batch_result.get("variant_rows", [])
    fieldnames = batch_result.get("variant_fields") or _variant_fields(UNKNOWN_TOOL_NAME, UNKNOWN_TOOL_NAME)
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return tsv_path


def write_comparison_excel(batch_result: dict[str, Any], output_path: str | Path) -> Path:
    """Write comparison summary to an Excel file with auto-adjusted column widths.

    Excel output format matches the original Sequencher comparison format.
    """
    excel_path = Path(output_path)
    rows = batch_result.get("summary_rows", [])

    if not rows:
        logger.warning("No summary rows to write to Excel")
        return excel_path

    results_df = pd.DataFrame(rows)
    # Enforce the tool-parameterized column order (matches the TSV header).
    fieldnames = batch_result.get("summary_fields") or _summary_fields(UNKNOWN_TOOL_NAME, UNKNOWN_TOOL_NAME)
    results_df = results_df.reindex(columns=fieldnames)
    results_df = results_df.sort_values(by=["Sample ID"])

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        results_df.to_excel(writer, sheet_name="Comparison Results", index=False)

    logger.success(f"Results saved to Excel file: {excel_path}")
    return excel_path


def _validate_input_file(path: Path, label: str) -> bool:
    """Validate that an input file exists and contains valid JSON.

    Args:
        path: Path to validate.
        label: Human-readable label for error messages (e.g. "Batch A").

    Returns:
        True if valid, False otherwise (errors logged via Loguru).
    """
    if not path.exists():
        logger.error("{} input file not found: {}", label, path)
        return False
    try:
        text = path.read_text(encoding="utf-8")
        json.loads(text)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        logger.error("Invalid JSON in {} input file {}: {}", label, path, exc)
        return False
    except OSError as exc:
        logger.error("Cannot read {} input file {}: {}", label, path, exc)
        return False
    return True


def compare_batch_files(
    path_a: str | Path,
    path_b: str | Path,
    output_dir: str | Path,
    batch_id: str,
    *,
    legacy_format: bool = False,
) -> dict[str, Any]:
    batch_result = SampleComparator.compare_batches(
        load_sample_batch(path_a),
        load_sample_batch(path_b),
        batch_id=batch_id,
        legacy_format=legacy_format,
    )
    batch_result["batch_id"] = batch_id
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_comparison_json(batch_result, out / f"{batch_id}.json")
    write_comparison_tsv(batch_result, out / f"{batch_id}.tsv")
    write_variant_level_tsv(batch_result, out / f"{batch_id}_variant_level.tsv")
    write_comparison_excel(batch_result, out / f"{batch_id}.xlsx")
    return batch_result


def main() -> None:
    parser = ArgumentParser(description="Compare two sample batch JSON files and produce JSON/TSV outputs")
    parser.add_argument("-a", "--batch-a", type=Path, required=True, help="Path to first batch JSON")
    parser.add_argument("-b", "--batch-b", type=Path, required=True, help="Path to second batch JSON")
    parser.add_argument("-o", "--output-dir", type=Path, required=True, help="Output directory for comparison files")
    parser.add_argument("--batch-id", type=str, required=True, help="Batch identifier used as filename prefix")
    parser.add_argument(
        "--legacy-format",
        action="store_true",
        help="Emit the legacy 13-column summary layout "
        "(Pipeline/Sequencher labels, single Variants Count, merged Flagged columns)",
    )
    args = parser.parse_args()

    if not _validate_input_file(args.batch_a, "Batch A") or not _validate_input_file(args.batch_b, "Batch B"):
        sys.exit(1)

    result = compare_batch_files(
        path_a=args.batch_a,
        path_b=args.batch_b,
        output_dir=args.output_dir,
        batch_id=args.batch_id,
        legacy_format=args.legacy_format,
    )
    s = result["summary"]
    logger.info(
        f"Matched: {s['matched']}, Concordant: {s['concordant']}, "
        f"Discordant: {s['discordant']}, Flagged: {s['flagged']}",
    )


if __name__ == "__main__":
    main()
