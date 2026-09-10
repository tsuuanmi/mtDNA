"""Indexed mtDNA-profile family grouping for TNLS verification."""

from collections.abc import Iterable
from itertools import combinations, product
from time import monotonic
from typing import Any

from loguru import logger

from src.core.variants import get_overlap, iupac_bases_compatible, normalize_position, pos_base

type Interval = tuple[int, int]
type Coverage = tuple[Interval, ...]
type VariantCall = tuple[str, str, str, int]
type Profile = tuple[VariantCall, ...]
type PositionKey = tuple[str, ...]
type ProfileIndex = dict[PositionKey, dict[Profile, list[str]]]

_PROGRESS_INTERVAL_SECONDS = 10.0


class _DisjointSet:
    """Disjoint-set structure used to build transitive profile families."""

    def __init__(self, sample_ids: Iterable[str]) -> None:
        self._parent = {sample_id: sample_id for sample_id in sample_ids}

    def _root(self, sample_id: str) -> str:
        root = sample_id
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[sample_id] != root:
            self._parent[sample_id], sample_id = root, self._parent[sample_id]
        return root

    def union(self, left: str, right: str) -> None:
        left_root = self._root(left)
        right_root = self._root(right)
        if left_root != right_root:
            self._parent[left_root] = right_root

    def union_groups(self, left: list[str], right: list[str]) -> None:
        """Union two non-empty groups connected by compatible profiles."""
        anchor = left[0]
        for sample_id in left[1:]:
            self.union(anchor, sample_id)
        for sample_id in right:
            self.union(anchor, sample_id)

    def components(self) -> list[list[str]]:
        grouped: dict[str, list[str]] = {}
        for sample_id in self._parent:
            grouped.setdefault(self._root(sample_id), []).append(sample_id)
        families = [sorted(members) for members in grouped.values()]
        families.sort(key=lambda members: members[0])
        return families


def _sample_coverage(sample: dict[str, Any], regions: list[str]) -> Coverage:
    intervals = ((start, end) for region in regions for start, end in sample.get("intervals", {}).get(region, []))
    return tuple(sorted(intervals))


def _coverage_bp(coverage: Coverage) -> int:
    return sum(end - start + 1 for start, end in coverage)


def sample_coverage_bp(sample: dict[str, Any], regions: list[str]) -> int:
    """Return the total sequenced base pairs in selected regions."""
    return _coverage_bp(_sample_coverage(sample, regions))


def _coverage_overlap(left: Coverage, right: Coverage) -> Coverage:
    overlap = get_overlap(
        [list(interval) for interval in left],
        [list(interval) for interval in right],
    )
    return tuple(sorted((start, end) for start, end in overlap))


def _profile(sample: dict[str, Any], coverage: Coverage) -> Profile:
    """Build the canonical typed variant profile inside a coverage range."""
    calls: dict[str, VariantCall] = {}
    for variant_type, records in sample.get("variants", {}).items():
        for record in records:
            position = record.get("pos")
            if position is None:
                continue
            normalized = normalize_position(position)
            base_position = pos_base(normalized)
            if not any(start <= base_position <= end for start, end in coverage):
                continue
            key = f"{variant_type}_{normalized}"
            calls[key] = (key, variant_type, str(record.get("seq", "")), base_position)
    return tuple(sorted(calls.values()))


def _restrict_profile(profile: Profile, coverage: Coverage) -> Profile:
    """Project a canonical profile onto a common coverage range."""
    return tuple(call for call in profile if any(start <= call[3] <= end for start, end in coverage))


def _positions(profile: Profile) -> PositionKey:
    return tuple(call[0] for call in profile)


def _profiles_compatible(left: Profile, right: Profile) -> bool:
    """Check exact typed positions with IUPAC-tolerant SNP alleles."""
    if _positions(left) != _positions(right):
        return False
    for left_call, right_call in zip(left, right, strict=True):
        _, left_type, left_sequence, _ = left_call
        _, right_type, right_sequence, _ = right_call
        if left_type != right_type:
            return False
        if left_type == "snps":
            if not iupac_bases_compatible(left_sequence, right_sequence):
                return False
        elif left_sequence != right_sequence:
            return False
    return True


def _index_profiles(
    sample_ids: list[str],
    profiles: dict[str, Profile],
    coverage: Coverage,
) -> ProfileIndex:
    index: ProfileIndex = {}
    for sample_id in sample_ids:
        profile = _restrict_profile(profiles[sample_id], coverage)
        index.setdefault(_positions(profile), {}).setdefault(profile, []).append(sample_id)
    return index


def _link_profile_indexes(
    families: _DisjointSet,
    left: ProfileIndex,
    right: ProfileIndex,
    *,
    same_index: bool,
) -> int:
    comparisons = 0
    for position_key in left.keys() & right.keys():
        left_profiles = left[position_key]
        right_profiles = right[position_key]
        profile_pairs = (
            combinations(left_profiles.items(), 2)
            if same_index
            else product(left_profiles.items(), right_profiles.items())
        )
        for (left_profile, left_ids), (right_profile, right_ids) in profile_pairs:
            comparisons += 1
            if _profiles_compatible(left_profile, right_profile):
                families.union_groups(left_ids, right_ids)
    return comparisons


def _collapse_exact_profiles(
    samples: dict[str, Any],
    regions: list[str],
    minimum: int,
    families: _DisjointSet,
) -> tuple[dict[Coverage, list[str]], dict[str, Profile]]:
    exact_groups: dict[tuple[Coverage, Profile], list[str]] = {}
    profiles: dict[str, Profile] = {}
    for sample_id, sample in samples.items():
        coverage = _sample_coverage(sample, regions)
        profile = _profile(sample, coverage)
        profiles[sample_id] = profile
        exact_groups.setdefault((coverage, profile), []).append(sample_id)

    coverage_groups: dict[Coverage, list[str]] = {}
    for (coverage, profile), members in exact_groups.items():
        if _coverage_bp(coverage) < minimum:
            continue
        # Invalid/unknown SNP codes are not compatible even with themselves.
        # Keep those samples separate to preserve the matching contract.
        if not _profiles_compatible(profile, profile):
            coverage_groups.setdefault(coverage, []).extend(members)
            continue
        representative = members[0]
        coverage_groups.setdefault(coverage, []).append(representative)
        for member in members[1:]:
            families.union(representative, member)
    return coverage_groups, profiles


def _coverage_jobs(
    coverage_groups: dict[Coverage, list[str]],
    minimum: int,
) -> list[tuple[Coverage, Coverage, Coverage]]:
    coverages = list(coverage_groups)
    jobs = [(coverage, coverage, coverage) for coverage in coverages]
    for left, right in combinations(coverages, 2):
        overlap = _coverage_overlap(left, right)
        if _coverage_bp(overlap) >= minimum:
            jobs.append((left, right, overlap))
    return jobs


def group_families(
    samples: dict[str, Any],
    regions: list[str],
    minimum: int,
) -> list[list[str]]:
    """Group samples by common-overlap mtDNA profile.

    Exact coverage/profile duplicates are collapsed first. Remaining work is
    indexed by coverage, typed variant positions, and canonical profiles, so
    compatibility checks operate on unique profiles rather than every sample
    pair. Compatible links remain transitive through the disjoint set.
    """
    families = _DisjointSet(samples)
    coverage_groups, profiles = _collapse_exact_profiles(samples, regions, minimum, families)
    representative_count = sum(len(group) for group in coverage_groups.values())
    jobs = _coverage_jobs(coverage_groups, minimum)
    logger.info(
        f"Family grouping: indexed {len(samples)} samples as "
        f"{representative_count} profile representatives across {len(coverage_groups)} coverage groups"
    )

    comparisons = 0
    last_progress = monotonic()
    for job_index, (left_coverage, right_coverage, overlap) in enumerate(jobs, start=1):
        left_index = _index_profiles(coverage_groups[left_coverage], profiles, overlap)
        same_index = left_coverage == right_coverage
        right_index = left_index if same_index else _index_profiles(coverage_groups[right_coverage], profiles, overlap)
        comparisons += _link_profile_indexes(
            families,
            left_index,
            right_index,
            same_index=same_index,
        )
        now = monotonic()
        if now - last_progress >= _PROGRESS_INTERVAL_SECONDS or job_index == len(jobs):
            logger.info(
                f"  Family progress: {job_index}/{len(jobs)} coverage comparisons ({job_index / len(jobs):.0%})"
            )
            last_progress = now

    logger.info(f"Family grouping: checked {comparisons} unique profile pairs")
    return families.components()
