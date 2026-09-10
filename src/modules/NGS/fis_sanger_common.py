"""Shared helpers for staged FIS-to-Sanger comparison preparation."""

import json
from pathlib import Path
from typing import Any

from loguru import logger

from src.core.flagging import flag_variants, recompute_sample_flags
from src.core.models import Sample, Tool, Variant, validate_position
from src.core.variants import is_position_in_intervals, variant_to_dict


def load_json_object(path: Path) -> dict[str, Any]:
    """Load a JSON object and reject non-object roots."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        msg = f"Expected a JSON object in {path}"
        raise TypeError(msg)
    return data


def normalize_intervals(raw_intervals: object) -> dict[str, list[list[int]]] | None:
    """Normalize a per-region interval mapping into the Sample contract."""
    if not isinstance(raw_intervals, dict):
        return None

    normalized: dict[str, list[list[int]]] = {}
    for region, raw_ranges in raw_intervals.items():
        if not isinstance(raw_ranges, list):
            continue
        ranges: list[list[int]] = []
        for raw_range in raw_ranges:
            if not isinstance(raw_range, (list, tuple)) or len(raw_range) != 2:  # noqa: PLR2004
                continue
            try:
                start, end = int(raw_range[0]), int(raw_range[1])
            except (TypeError, ValueError):
                continue
            ranges.append([start, end])
        if ranges:
            normalized[str(region)] = ranges
    return normalized or None


def flatten_intervals(intervals: dict[str, list[list[int]]] | None) -> list[tuple[int, int]]:
    """Flatten Sample intervals for position filtering."""
    if not intervals:
        return []
    return [(start, end) for ranges in intervals.values() for start, end in ranges]


def normalize_variant(raw_variant: dict[str, Any], sample_id: str) -> Variant | None:
    """Normalize one external variant, warning and skipping invalid positions."""
    raw_pos = raw_variant.get("pos", raw_variant.get("position"))
    if raw_pos is None:
        logger.warning("Skipping variant without a position for {}", sample_id)
        return None
    try:
        # Legacy Sanger JSON stores insertion positions as JSON floats. Convert
        # explicitly at this boundary; the source has already lost trailing-zero precision.
        pos = validate_position(str(raw_pos) if isinstance(raw_pos, float) else raw_pos)
    except (TypeError, ValueError):
        logger.warning("Skipping invalid variant position for {}: {}", sample_id, raw_pos)
        return None

    ref = str(raw_variant.get("ref", "")).strip().upper()
    seq = str(raw_variant.get("seq", raw_variant.get("alt", ""))).strip().upper()
    if not ref or not seq or seq == "N":
        return None
    return Variant(pos=pos, ref=ref, seq=seq)


def filter_variants_to_intervals(
    variants: list[Variant],
    intervals: dict[str, list[list[int]]] | None,
) -> list[Variant]:
    """Keep variants within the declared intervals; no intervals means no comparable calls."""
    flat_intervals = flatten_intervals(intervals)
    if not flat_intervals:
        return []
    return [variant for variant in variants if is_position_in_intervals(variant.pos, flat_intervals)]


def build_comparison_sample(
    *,
    sample_id: str,
    variants: list[Variant],
    source_tool: Tool,
    intervals: dict[str, list[list[int]]] | None,
    batch_id: str | None = None,
    information: dict[str, Any] | None = None,
) -> Sample:
    """Build a canonical comparison-scoped Sample with core flags."""
    variant_flags = flag_variants(variants)
    sample_flags = recompute_sample_flags([], variants, intervals, variant_flags)
    return Sample(
        sample_id=sample_id,
        variants=variants,
        source_tool=source_tool,
        intervals=intervals,
        sample_flags=sample_flags,
        variant_flags=variant_flags,
        batch_id=batch_id,
        information=information,
    )


def write_sample_batch(samples: list[Sample], output_path: Path) -> Path:
    """Write Samples as a keyed batch JSON accepted by load_sample_batch()."""
    keyed: dict[str, dict[str, Any]] = {}
    for sample in samples:
        if sample.sample_id in keyed:
            msg = f"Duplicate comparison sample key: {sample.sample_id}"
            raise ValueError(msg)
        sample_data = sample.model_dump(mode="json", exclude={"sample_id"}, exclude_none=True)
        sample_data["variants"] = [
            variant_data
            for variant_data in (variant_to_dict(variant) for variant in sample.variants)
            if variant_data is not None
        ]
        keyed[sample.sample_id] = sample_data

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(keyed, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.success("Saved prepared sample batch: {} ({} samples)", output_path, len(samples))
    return output_path
