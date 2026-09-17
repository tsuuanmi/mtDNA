"""Callability from declared intervals and normalized variant alleles.

Profiles are computed once per distinct plain profile — identical sorted
normalized variant keys plus identical per-region declared intervals — and
aggregated with the group's member count. Callable bases and insertion anchors
are stored as inclusive run pairs plus the few excluded positions, never as
per-position sets; the run arithmetic is exactly equivalent to per-position
interval membership (validated against the previous per-position
implementation on the full current cohort).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from itertools import pairwise
from typing import TYPE_CHECKING

from src.core.models import Sample, validate_position
from src.core.variants import pos_base, pos_sort_key
from src.modules.statistics.population.models import REGIONS

if TYPE_CHECKING:
    from collections.abc import Mapping

RESOLVED = frozenset("ACGTRYSWKMBDHV")
UNRESOLVED = frozenset("Nacgtryswkmbdhvn")
VariantKey = tuple[int | str, str, str]
RegionIntervals = tuple[tuple[int, ...], ...]
PlainProfileSignature = tuple[tuple[VariantKey, ...], tuple[RegionIntervals, ...]]


def variant_key(key: VariantKey) -> VariantKey:
    """Normalize position formatting only; never reinterpret allele case or biology."""
    pos, ref, alt = key
    if isinstance(pos, bool) or not re.fullmatch(r"\s*\d+(?:\.\d+)?\s*", str(pos)):
        msg = "Variant positions must be positive integer bases or insertion indices"
        raise ValueError(msg)
    position = validate_position(pos)
    if isinstance(position, str):
        base, index = position.split(".")
        if int(index) < 1:
            msg = "Insertion indices must be positive"
            raise ValueError(msg)
        position = f"{int(base)}.{int(index)}"
    if pos_base(position) < 1:
        msg = "Variant base positions must be positive"
        raise ValueError(msg)
    return position, ref, alt


def key_order(key: VariantKey) -> tuple[tuple[int, int], str, str]:
    return pos_sort_key(key[0]), key[1], key[2]


def plain_profile_signature(sample: Sample) -> PlainProfileSignature:
    """Group key of the compact plain profile: sorted raw keys plus declared intervals."""
    keys = tuple(sorted(((variant.pos, variant.ref, variant.seq) for variant in sample.variants), key=key_order))
    intervals = tuple(
        tuple(sorted(tuple(pair) for pair in (sample.intervals or {}).get(region, ()))) for region in REGIONS
    )
    return keys, intervals


@dataclass
class _RegionCalls:
    """Verified callable variant keys plus range-based callability for one region."""

    keys: tuple[VariantKey, ...] = ()
    base_runs: tuple[tuple[int, int], ...] = ()
    anchor_runs: tuple[tuple[int, int], ...] = ()
    excluded_bases: frozenset[int] = frozenset()
    excluded_anchors: frozenset[int] = frozenset()
    complete: bool = False


@dataclass
class GroupProfile:
    """Compact callability shared by every sample of one identical plain profile.

    Callability, verified keys, completeness and reasons derive only from the
    variant keys and declared intervals, so identical plain profiles are
    computed once and aggregated with the group's member count. The finalized
    HV sequence strings are deliberately not consumed: historical sequence
    generation is known to disagree with the normalized variant lists for part
    of the cohort, so population statistics trust the pipeline's declared
    callable intervals and normalized variant alleles instead.
    """

    regions: dict[str, _RegionCalls]
    reasons: list[str]
    issues: list[tuple[str, str, str]] = field(default_factory=list)


@dataclass
class _VerifiedCalls:
    """Verified keys plus the positions their representation leaves unverifiable or unresolved."""

    keys: tuple[VariantKey, ...] = ()
    unverifiable_bases: set[int] = field(default_factory=set)
    unverifiable_anchors: set[int] = field(default_factory=set)
    unverifiable_keys: list[VariantKey] = field(default_factory=list)
    unresolved_bases: set[int] = field(default_factory=set)
    unresolved_anchors: set[int] = field(default_factory=set)
    boundary_insertion: bool = False


@dataclass(frozen=True)
class _ClippedRegion:
    """Interval runs clipped to one region plus full-span coverage facts."""

    base_runs: tuple[tuple[int, int], ...]
    anchor_runs: tuple[tuple[int, int], ...]
    full_base_span: bool
    full_anchor_span: bool


def _merged_runs(runs: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Merge inclusive runs into disjoint, non-adjacent runs with the same covered union.

    A sample's declared intervals may overlap or repeat, so per-interval runs
    must be merged before any per-position counting: overlapping runs would
    otherwise count their shared positions once per run instead of once per
    sample.
    """
    merged: list[tuple[int, int]] = []
    for run_start, run_end in runs:
        if merged and run_start <= merged[-1][1] + 1:
            if run_end > merged[-1][1]:
                merged[-1] = (merged[-1][0], run_end)
        else:
            merged.append((run_start, run_end))
    return tuple(merged)


def _clipped_region(
    region: str, intervals: RegionIntervals, start: int, end: int, cache: dict[object, _ClippedRegion]
) -> _ClippedRegion:
    cache_key = (region, intervals)
    clipped = cache.get(cache_key)
    if clipped is None:
        base_runs = []
        for interval_start, interval_end in intervals:
            low, high = max(start, interval_start), min(end, interval_end)
            if low <= high:
                base_runs.append((low, high))
        # Anchor runs are derived per interval before merging: merging first would
        # erase the boundary where one interval ends and another begins, where no
        # insertion anchor is provable.
        anchor_runs = [(low, high - 1) for low, high in base_runs if low < high]
        merged_base = _merged_runs(base_runs)
        merged_anchors = _merged_runs(anchor_runs)
        clipped = _ClippedRegion(
            base_runs=merged_base,
            anchor_runs=merged_anchors,
            full_base_span=sum(high - low + 1 for low, high in merged_base) == end - start + 1,
            full_anchor_span=sum(high - low + 1 for low, high in merged_anchors) == end - start,
        )
        cache[cache_key] = clipped
    return clipped


def _resolved_allele(alt: str) -> bool:
    """A resolved call is an accepted uppercase code or a called deletion."""
    return alt == "-" or all(base in RESOLVED for base in alt)


def _covered(runs: tuple[tuple[int, int], ...], position: int) -> bool:
    return any(start <= position <= end for start, end in runs)


def _canonical(key: VariantKey, cache: dict[VariantKey, VariantKey]) -> VariantKey:
    canonical = cache.get(key)
    if canonical is None:
        canonical = variant_key(key)
        cache[key] = canonical
    return canonical


def _verified_calls(
    keys: list[VariantKey], end: int, reference: str, key_cache: dict[VariantKey, VariantKey]
) -> _VerifiedCalls:
    """Verify normalized representation against the configured reference.

    A mis-referenced, conflicting, reference-equal or unsupported entry cannot be
    trusted at its position: it never counts as a carrier and leaves that position
    (or insertion anchor) non-callable for this group instead of silently passing.
    """
    grouped: dict[int | str, set[VariantKey]] = {}
    for key in keys:
        canonical = _canonical(key, key_cache)
        grouped.setdefault(canonical[0], set()).add(canonical)
    calls = _VerifiedCalls()
    verified: list[VariantKey] = []
    for position, gkeys in grouped.items():
        anchor = isinstance(position, str)
        for key in sorted(gkeys):
            _, ref, alt = key
            valid = len(gkeys) == 1
            valid = valid and ((ref == "-") if anchor else (len(ref) == 1 and ref == reference[position - 1]))
            valid = valid and bool(alt) and ((alt != "-") if anchor else len(alt) == 1)
            valid = valid and alt != ref
            valid = valid and (alt == "-" or set(alt) <= RESOLVED | UNRESOLVED)
            if valid:
                verified.append(key)
                if anchor and pos_base(position) == end:
                    calls.boundary_insertion = True
                if not _resolved_allele(alt):
                    base = pos_base(position)
                    if anchor:
                        # An unresolved insertion blocks both its anchor and its base position.
                        calls.unresolved_anchors.add(base)
                        calls.unresolved_bases.add(base)
                    else:
                        calls.unresolved_bases.add(position)
            elif anchor:
                calls.unverifiable_anchors.add(pos_base(position))
                calls.unverifiable_keys.append(key)
            else:
                calls.unverifiable_bases.add(position)
                calls.unverifiable_keys.append(key)
    calls.unverifiable_keys.sort(key=key_order)
    calls.keys = tuple(sorted(verified, key=key_order))
    return calls


def _callable_keys(
    calls: _VerifiedCalls,
    clipped: _ClippedRegion,
    excluded_bases: set[int],
    excluded_anchors: set[int],
) -> tuple[VariantKey, ...]:
    """Keep verified keys whose position is callable with a resolved allele."""
    carriers = []
    for key in calls.keys:
        if not _resolved_allele(key[2]):
            continue
        if isinstance(key[0], str):
            anchor = pos_base(key[0])
            if (
                _covered(clipped.anchor_runs, anchor)
                and anchor not in excluded_bases
                and anchor + 1 not in excluded_bases
                and anchor not in excluded_anchors
            ):
                carriers.append(key)
        elif _covered(clipped.base_runs, key[0]) and key[0] not in excluded_bases:
            carriers.append(key)
    return tuple(carriers)


def _interval_overlaps(intervals: RegionIntervals) -> list[str]:
    """Describe interval pairs that overlap or duplicate each other."""
    ordered = sorted(intervals)
    return [
        f"[{current[0]}, {current[1]}] overlaps or duplicates [{previous[0]}, {previous[1]}]"
        for previous, current in pairwise(ordered)
        if current[0] <= previous[1]
    ]


def _unverifiable_details(keys: list[VariantKey], reference: str) -> str:
    """Describe each unverifiable entry so the upstream record can be found and corrected."""
    parts = []
    for position, ref, alt in keys:
        if isinstance(position, str):
            parts.append(f"{position} {ref}>{alt} (insertion representation)")
        else:
            parts.append(f"{position} {ref}>{alt} (reference base {reference[position - 1]})")
    return "; ".join(parts)


def _region_issues(
    region: str, region_intervals: RegionIntervals, calls: _VerifiedCalls, reference: str
) -> list[tuple[str, str, str]]:
    """Collect fixable representation issues for one region of one group profile."""
    issues = []
    overlaps = _interval_overlaps(region_intervals)
    if overlaps:
        issues.append((region, "duplicated_or_overlapping_intervals", "; ".join(overlaps)))
    if calls.unverifiable_bases or calls.unverifiable_anchors:
        issues.append(
            (
                region,
                "unverifiable_variant_representation",
                _unverifiable_details(calls.unverifiable_keys, reference),
            )
        )
    return issues


def group_profile(
    signature: PlainProfileSignature,
    regions: Mapping[str, list[int]],
    reference: str,
    *,
    clip_cache: dict[object, _ClippedRegion],
    key_cache: dict[VariantKey, VariantKey],
) -> GroupProfile:
    """Derive compact callability for one distinct plain profile.

    Callable bases are the clipped interval runs minus excluded positions. A
    base position is an insertion anchor when one interval spans it and the
    next base, which the per-interval anchor runs encode exactly. Region
    completeness requires the full configured base span, every internal anchor
    and no represented insertion at the configured right boundary. Fixable
    upstream representation issues — missing, duplicated or overlapping
    intervals and unverifiable variant entries — are collected per region for
    the data-issues report.
    """
    keys, interval_sets = signature
    buckets: dict[str, list[VariantKey]] = {region: [] for region in REGIONS}
    for key in keys:
        base = pos_base(key[0])
        for region in REGIONS:
            start, end = regions[region]
            if start <= base <= end:
                buckets[region].append(key)
                break
    profile = GroupProfile({}, [])
    for index, region in enumerate(REGIONS):
        start, end = regions[region]
        calls = _verified_calls(buckets[region], end, reference, key_cache)
        region_intervals = interval_sets[index]
        for pair in region_intervals:
            if len(pair) != 2 or pair[0] < 1 or pair[0] > pair[1]:  # noqa: PLR2004 - endpoint pair
                msg = "Invalid declared callable interval"
                raise ValueError(msg)
        if not region_intervals:
            profile.regions[region] = _RegionCalls()
            profile.reasons.append(f"{region}:missing_intervals")
            profile.issues.append((region, "missing_intervals", "declared callable intervals are absent or empty"))
            continue
        clipped = _clipped_region(region, region_intervals, start, end, clip_cache)
        profile.issues.extend(_region_issues(region, region_intervals, calls, reference))
        excluded_bases = calls.unverifiable_bases | calls.unresolved_bases
        excluded_anchors = calls.unverifiable_anchors | calls.unresolved_anchors
        anchor_span_excluded = any(start <= anchor <= end - 1 for anchor in excluded_anchors)
        complete = (
            clipped.full_base_span
            and not excluded_bases
            and clipped.full_anchor_span
            and not anchor_span_excluded
            and not calls.boundary_insertion
        )
        profile.regions[region] = _RegionCalls(
            keys=_callable_keys(calls, clipped, excluded_bases, excluded_anchors),
            base_runs=clipped.base_runs,
            anchor_runs=clipped.anchor_runs,
            excluded_bases=frozenset(excluded_bases),
            excluded_anchors=frozenset(excluded_anchors),
            complete=complete,
        )
        unverifiable_reason = f"{region}:unverifiable_variant_representation"
        if calls.unverifiable_bases or calls.unverifiable_anchors:
            profile.reasons.append(unverifiable_reason)
        if not complete and unverifiable_reason not in profile.reasons:
            profile.reasons.append(f"{region}:non_callable_positions_or_insertion")
    return profile
