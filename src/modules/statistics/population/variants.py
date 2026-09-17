"""Population carrier frequencies using site-specific callable denominators."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from src.core.models import Variant
from src.core.variants import format_variants_simplified
from src.modules.statistics.population.coverage import GroupProfile, VariantKey, _covered, key_order
from src.modules.statistics.population.models import REGIONS, AnalysisOptions, VariantDiversity, VariantPopulationStat
from src.modules.statistics.population.statistics import frequency_interval

if TYPE_CHECKING:
    from src.modules.statistics.population.reference import ReferenceProvider

DOUBLETON_COUNT = 2
AF_THRESHOLDS = (0.001, 0.005, 0.01)


@dataclass
class VariantCounts:
    """Fixed-scope denominator counters, aggregated once per distinct profile group.

    Per-position callable denominators are kept as counters over inclusive run
    signatures plus per-position exclusion counts, expanded into per-position
    totals only at aggregation time.
    """

    carriers: dict[str, Counter[VariantKey]] = field(default_factory=lambda: {region: Counter() for region in REGIONS})
    base_runs: dict[str, Counter[tuple[tuple[int, int], ...]]] = field(
        default_factory=lambda: {region: Counter() for region in REGIONS}
    )
    anchor_runs: dict[str, Counter[tuple[tuple[int, int], ...]]] = field(
        default_factory=lambda: {region: Counter() for region in REGIONS}
    )
    base_excluded: dict[str, Counter[int]] = field(default_factory=lambda: {region: Counter() for region in REGIONS})
    anchor_excluded: dict[str, Counter[int]] = field(default_factory=lambda: {region: Counter() for region in REGIONS})

    def add(self, profile: GroupProfile, multiplicity: int = 1) -> None:
        for region in REGIONS:
            calls = profile.regions[region]
            carriers = self.carriers[region]
            for key in calls.keys:
                carriers[key] += multiplicity
            self.base_runs[region][calls.base_runs] += multiplicity
            self.anchor_runs[region][calls.anchor_runs] += multiplicity
            self._excluded(self.base_excluded[region], calls.base_runs, calls.excluded_bases, multiplicity)
            # A non-callable base removes its own insertion anchor and the previous one.
            killed_anchors = calls.excluded_anchors | {
                anchor for base in calls.excluded_bases for anchor in (base, base - 1)
            }
            self._excluded(self.anchor_excluded[region], calls.anchor_runs, killed_anchors, multiplicity)

    @staticmethod
    def _excluded(
        counter: Counter[int],
        runs: tuple[tuple[int, int], ...],
        positions: frozenset[int],
        multiplicity: int,
    ) -> None:
        for position in positions:
            if _covered(runs, position):
                counter[position] += multiplicity


def _position_counts(run_counts: Counter[tuple[tuple[int, int], ...]], excluded: Counter[int]) -> Counter[int]:
    """Expand run-signature counts into per-position callable counts."""
    positions: Counter[int] = Counter()
    for runs, multiplicity in run_counts.items():
        for start, end in runs:
            for position in range(start, end + 1):
                positions[position] += multiplicity
    for position, multiplicity in excluded.items():
        positions[position] -= multiplicity
    return positions


def aggregate_variants(
    counts: VariantCounts,
    options: AnalysisOptions,
    reference_sha256: str,
    regions: dict[str, list[int]],
    provider: ReferenceProvider | None = None,
) -> list[VariantPopulationStat]:
    rows = []
    for region in REGIONS:
        base_counts = _position_counts(counts.base_runs[region], counts.base_excluded[region])
        insertion_counts = _position_counts(counts.anchor_runs[region], counts.anchor_excluded[region])
        for key, count in sorted(counts.carriers[region].items(), key=lambda item: key_order(item[0])):
            position, ref, alt = key
            total = (
                insertion_counts[int(position.split(".")[0])] if isinstance(position, str) else base_counts[position]
            )
            if not 0 < count <= total:
                msg = "Carrier/callable count invariant violated"
                raise ValueError(msg)
            frequency = count / total
            annotations = (
                {}
                if provider is None
                else provider.annotate(
                    Variant(pos=position, ref=ref, seq=alt),
                    reference_sha256=reference_sha256,
                    region_definitions=regions,
                )
            )
            rows.append(
                VariantPopulationStat(
                    variant=format_variants_simplified([{"pos": position, "ref": ref, "seq": alt}]),
                    region=region,
                    position=position,
                    reference=ref,
                    alternate=alt,
                    carrier_count=count,
                    callable_sample_count=total,
                    population_frequency=frequency,
                    is_singleton=count == 1,
                    is_doubleton=count == DOUBLETON_COUNT,
                    carrier_count_le_2=count <= DOUBLETON_COUNT,
                    frequency_lt_0_001=frequency < AF_THRESHOLDS[0],
                    frequency_lt_0_005=frequency < AF_THRESHOLDS[1],
                    frequency_lt_0_01=frequency < AF_THRESHOLDS[2],
                    additional_frequency_flags={
                        str(threshold): frequency < threshold
                        for threshold in sorted(set(options.additional_frequency_thresholds))
                    },
                    reference_annotations=annotations,
                    **frequency_interval(count, total, options.confidence_level).model_dump(),
                )
            )
    return sorted(rows, key=lambda row: key_order((row.position, row.reference, row.alternate)))


def summarize_variants(rows: list[VariantPopulationStat], options: AnalysisOptions) -> VariantDiversity:
    return VariantDiversity(
        n_distinct_variants=len(rows),
        n_singleton_variants=sum(row.is_singleton for row in rows),
        n_doubleton_variants=sum(row.is_doubleton for row in rows),
        n_variants_carrier_count_le_2=sum(row.carrier_count_le_2 for row in rows),
        n_variants_af_lt_0_001=sum(row.frequency_lt_0_001 for row in rows),
        n_variants_af_lt_0_005=sum(row.frequency_lt_0_005 for row in rows),
        n_variants_af_lt_0_01=sum(row.frequency_lt_0_01 for row in rows),
        additional_frequency_counts={
            str(threshold): sum(row.population_frequency < threshold for row in rows)
            for threshold in sorted(set(options.additional_frequency_thresholds))
        },
    )
