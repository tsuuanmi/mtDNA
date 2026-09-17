"""Deterministic exact profiles and count-based grouping, never pairwise matching."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

from src.core.variants import format_variants_simplified
from src.modules.statistics.population.coverage import GroupProfile, VariantKey, key_order
from src.modules.statistics.population.models import (
    COMBINED,
    REGIONS,
    HaplotypePopulationStat,
    SampleHaplotype,
)
from src.modules.statistics.population.statistics import frequency_interval


def display(keys: tuple[VariantKey, ...]) -> str:
    return format_variants_simplified([{"pos": pos, "ref": ref, "seq": alt} for pos, ref, alt in keys])


def haplotype_id(keys: tuple[VariantKey, ...], scope: str, fingerprint: str) -> str:
    payload = json.dumps([scope, fingerprint, keys], separators=(",", ":"))
    return "HAP-" + hashlib.sha256(payload.encode()).hexdigest()


def group_haplotypes(
    groups: Sequence[tuple[GroupProfile, int]],
    sample_groups: Sequence[tuple[str, int]],
    fingerprint: str,
    confidence_level: float | None,
) -> tuple[list[SampleHaplotype], dict[str, list[HaplotypePopulationStat]]]:
    """Count haplotypes per distinct profile group with member multiplicity.

    Per-group artifacts (combined keys, display strings, stable IDs,
    completeness) are computed once and reused for every member sample mapping.
    """
    counts: dict[str, Counter[tuple[VariantKey, ...]]] = {scope: Counter() for scope in (*REGIONS, COMBINED)}
    region_displays: list[dict[str, str | None]] = []
    region_completes: list[dict[str, bool]] = []
    combined_displays: list[str | None] = []
    combined_ids: list[str | None] = []
    complete_flags: list[bool] = []
    for profile, multiplicity in groups:
        region_keys = {region: profile.regions[region].keys for region in REGIONS}
        keys = tuple(sorted((key for region in REGIONS for key in region_keys[region]), key=key_order))
        complete = all(profile.regions[region].complete for region in REGIONS)
        for region in REGIONS:
            if profile.regions[region].complete:
                counts[region][region_keys[region]] += multiplicity
        if complete:
            counts[COMBINED][keys] += multiplicity
        region_displays.append(
            {region: display(region_keys[region]) if profile.regions[region].complete else None for region in REGIONS}
        )
        region_completes.append({region: profile.regions[region].complete for region in REGIONS})
        combined_displays.append(display(keys) if complete else None)
        combined_ids.append(haplotype_id(keys, COMBINED, fingerprint) if complete else None)
        complete_flags.append(complete)
    mappings = [
        SampleHaplotype(
            sample_id=sample_id,
            hv1_haplotype=region_displays[group_index]["HV1"],
            hv2_haplotype=region_displays[group_index]["HV2"],
            hv3_haplotype=region_displays[group_index]["HV3"],
            combined_haplotype=combined_displays[group_index],
            haplotype_id=combined_ids[group_index],
            complete_hv_profile=complete_flags[group_index],
            region_complete=region_completes[group_index],
            incomplete_reasons=groups[group_index][0].reasons,
        )
        for sample_id, group_index in sample_groups
    ]
    tables = {}
    for scope, scope_groups in counts.items():
        total = sum(scope_groups.values())
        tables[scope] = [
            HaplotypePopulationStat(
                scope=scope,
                haplotype_id=haplotype_id(keys, scope, fingerprint),
                canonical_haplotype=display(keys),
                count=count,
                frequency=count / total,
                is_singleton=count == 1,
                **frequency_interval(count, total, confidence_level).model_dump(),
            )
            for keys, count in sorted(scope_groups.items(), key=lambda item: (-item[1], tuple(map(key_order, item[0]))))
        ]
    return mappings, tables
