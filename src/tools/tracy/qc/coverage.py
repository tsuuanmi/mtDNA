"""Artifact-independent operations on trusted trace coverage and variants."""

from dataclasses import dataclass

from src.core.flagging import flag_variants, recompute_sample_flags
from src.core.models import Sample, Variant
from src.tools.tracy.qc.models import QCFinding
from src.tools.tracy.utils import variant_overlaps_range


def subtract_ranges(
    intervals: dict[str, list[list[int]]] | None,
    ranges: tuple[tuple[int, int], ...],
) -> dict[str, list[list[int]]] | None:
    """Subtract inclusive exclusions, retaining explicit empty coverage as empty."""
    if intervals is None:
        return None
    retained_by_region: dict[str, list[list[int]]] = {}
    for region, spans in intervals.items():
        retained = [list(span) for span in spans]
        for start, end in sorted(ranges):
            next_retained: list[list[int]] = []
            for left, right in retained:
                if end < left or start > right:
                    next_retained.append([left, right])
                    continue
                if left < start:
                    next_retained.append([left, start - 1])
                if end < right:
                    next_retained.append([end + 1, right])
            retained = next_retained
        if retained:
            retained_by_region[region] = retained
    return retained_by_region


@dataclass(frozen=True)
class ExcludedObservation:
    """Accepted ETL representation removed by an evidence finding."""

    variant: Variant
    finding: QCFinding


@dataclass(frozen=True)
class CoverageResult:
    """Trusted trace and its excluded canonical observations."""

    sample: Sample
    excluded: tuple[ExcludedObservation, ...]


def apply_findings(sample: Sample, findings: tuple[QCFinding, ...]) -> CoverageResult:
    """Apply only EXCLUDE decisions; REVIEW never changes trusted evidence."""
    exclusions = tuple(finding for finding in findings if finding.decision == "exclude")
    if not exclusions:
        return CoverageResult(sample, ())
    intervals = subtract_ranges(sample.intervals, tuple((f.start, f.end) for f in exclusions))
    retained: list[Variant] = []
    excluded: list[ExcludedObservation] = []
    for variant in sample.variants:
        finding = next(
            (f for f in exclusions if variant_overlaps_range(variant.pos, variant.ref, f.start, f.end)),
            None,
        )
        if finding is None:
            retained.append(variant)
        else:
            excluded.append(ExcludedObservation(variant, finding))
    if not excluded and intervals == sample.intervals:
        return CoverageResult(sample, ())
    variant_flags = flag_variants(retained)
    sample_flags = recompute_sample_flags(sample.sample_flags, retained, intervals, variant_flags)
    trusted = sample.model_copy(
        update={
            "variants": retained,
            "intervals": intervals,
            "variant_flags": variant_flags,
            "sample_flags": sample_flags,
        }
    )
    return CoverageResult(trusted, tuple(excluded))
