"""Apply Tracy noise ranges to completed per-trace samples.

Masking is separate from trace QC and variant calling: QC discovers ranges,
ETL creates canonical variants, and this module removes overlapping variants and
coverage before traces are merged into final output.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from src.core.flagging import flag_variants, recompute_sample_flags
from src.tools.tracy.quality_control import ExcludedVariant, NoiseRange
from src.tools.tracy.utils import variant_overlaps_range

if TYPE_CHECKING:
    from src.core.models import Sample, Variant


@dataclass(frozen=True)
class NoiseMaskResult:
    """A per-trace sample and variants excluded by its QC ranges."""

    sample: Sample
    excluded_variants: tuple[ExcludedVariant, ...]


def _matching_range(variant: Variant, ranges: tuple[NoiseRange, ...]) -> NoiseRange | None:
    """Return the first merged noise range overlapping a canonical variant."""
    return next(
        (
            noise_range
            for noise_range in ranges
            if variant_overlaps_range(
                variant.pos,
                variant.ref,
                noise_range.start,
                noise_range.end,
            )
        ),
        None,
    )


def _subtract_noise_ranges(
    intervals: dict[str, list[list[int]]] | None,
    ranges: tuple[NoiseRange, ...],
) -> dict[str, list[list[int]]] | None:
    """Return coverage intervals with inclusive noise ranges removed."""
    if intervals is None:
        return None

    sorted_ranges = sorted(
        (min(noise_range.start, noise_range.end), max(noise_range.start, noise_range.end))
        for noise_range in ranges
    )
    retained_by_region: dict[str, list[list[int]]] = {}
    for region, spans in intervals.items():
        retained = [list(span) for span in spans]
        for range_start, range_end in sorted_ranges:
            next_retained: list[list[int]] = []
            for span_start, span_end in retained:
                if range_end < span_start or range_start > span_end:
                    next_retained.append([span_start, span_end])
                    continue
                if span_start < range_start:
                    next_retained.append([span_start, range_start - 1])
                if range_end < span_end:
                    next_retained.append([range_end + 1, span_end])
            retained = next_retained
        if retained:
            retained_by_region[region] = retained
    return retained_by_region


def apply_noise_mask(sample: Sample, ranges: tuple[NoiseRange, ...]) -> NoiseMaskResult:
    """Exclude noisy-range variants and coverage, then recompute flags."""
    if not ranges:
        return NoiseMaskResult(sample, ())

    masked_intervals = _subtract_noise_ranges(sample.intervals, ranges)
    retained: list[Variant] = []
    excluded: list[ExcludedVariant] = []
    for variant in sample.variants:
        noise_range = _matching_range(variant, ranges)
        if noise_range is None:
            retained.append(variant)
            continue
        excluded.append(
            ExcludedVariant(
                pos=variant.pos,
                ref=variant.ref,
                seq=variant.seq,
                noise_range=noise_range,
            ),
        )

    if not excluded and masked_intervals == sample.intervals:
        return NoiseMaskResult(sample, ())

    variant_flags = flag_variants(retained)
    sample_flags = recompute_sample_flags(
        sample.sample_flags,
        retained,
        masked_intervals,
        variant_flags,
    )
    masked_sample = sample.model_copy(
        update={
            "variants": retained,
            "intervals": masked_intervals,
            "variant_flags": variant_flags,
            "sample_flags": sample_flags,
        },
    )
    return NoiseMaskResult(masked_sample, tuple(excluded))
