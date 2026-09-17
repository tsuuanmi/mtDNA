"""Population aggregation orchestration over an approved cohort snapshot."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.modules.statistics.population.coverage import (
    GroupProfile,
    PlainProfileSignature,
    VariantKey,
    _ClippedRegion,
    group_profile,
    plain_profile_signature,
)
from src.modules.statistics.population.haplotypes import group_haplotypes
from src.modules.statistics.population.models import (
    COMBINED,
    REGIONS,
    AnalysisOptions,
    DataIssue,
    PopulationCohort,
    PopulationResult,
    PopulationSummary,
)
from src.modules.statistics.population.statistics import diversity
from src.modules.statistics.population.variants import VariantCounts, aggregate_variants, summarize_variants

if TYPE_CHECKING:
    from src.modules.statistics.population.reference import ReferenceProvider

_PROGRESS_INTERVAL = 10_000


def _grouped_profiles(
    cohort: PopulationCohort, regions: dict[str, list[int]], reference: str
) -> tuple[list[tuple[GroupProfile, int]], list[tuple[str, int]]]:
    """Compute each distinct plain profile once and track its member count.

    Identical plain profiles — same sorted normalized variant keys plus same
    per-region declared intervals — share one computed group profile, so later
    analysis works with the group and its multiplicity instead of specific
    samples. Errors name the first sample that surfaces them.
    """
    key_cache: dict[VariantKey, VariantKey] = {}
    clip_cache: dict[object, _ClippedRegion] = {}
    profiles: list[GroupProfile] = []
    member_counts: list[int] = []
    index_of: dict[PlainProfileSignature, int] = {}
    sample_groups: list[tuple[str, int]] = []
    total = len(cohort.samples)
    for index, sample in enumerate(sorted(cohort.samples, key=lambda sample: sample.sample_id), start=1):
        try:
            signature = plain_profile_signature(sample)
        except (ValueError, TypeError) as error:
            msg = f"Sample {sample.sample_id}: {error}"
            raise type(error)(msg) from error
        group_index = index_of.get(signature)
        if group_index is None:
            try:
                profile = group_profile(signature, regions, reference, clip_cache=clip_cache, key_cache=key_cache)
            except (ValueError, TypeError) as error:
                msg = f"Sample {sample.sample_id}: {error}"
                raise type(error)(msg) from error
            group_index = len(profiles)
            index_of[signature] = group_index
            profiles.append(profile)
            member_counts.append(0)
        member_counts[group_index] += 1
        sample_groups.append((sample.sample_id, group_index))
        if index % _PROGRESS_INTERVAL == 0 or index == total:
            logger.info("Population cohort: built callability profiles for {}/{} samples", index, total)
    return list(zip(profiles, member_counts, strict=True)), sample_groups


def _sample_data_issues(
    groups: list[tuple[GroupProfile, int]], sample_groups: list[tuple[str, int]]
) -> tuple[list[DataIssue], dict[str, int], int]:
    """Expand per-group representation issues into one row per sample and issue.

    Rows keep sample-id order so upstream records can be found and corrected;
    counts report the number of affected samples per issue type.
    """
    rows = [
        DataIssue(sample_id=sample_id, issue=issue, region=region, details=details)
        for sample_id, group_index in sample_groups
        for region, issue, details in groups[group_index][0].issues
    ]
    counts = Counter(issue for _, issue in {(row.sample_id, row.issue) for row in rows})
    return rows, dict(sorted(counts.items())), len({row.sample_id for row in rows})


def analyze_population(  # noqa: C901, PLR0912, PLR0915 - orchestration includes explicit output invariants
    cohort: PopulationCohort,
    *,
    options: AnalysisOptions | None = None,
    reference_path: str | Path | None = None,
    reference_provider: ReferenceProvider | None = None,
    analysis_timestamp: datetime | None = None,
) -> PopulationResult:
    options = options or AnalysisOptions()
    if any(not math.isfinite(value) or not 0 < value <= 1 for value in options.additional_frequency_thresholds):
        msg = "Additional frequency thresholds must be finite values in (0, 1]"
        raise ValueError(msg)
    settings = get_settings()
    regions = {region: list(settings.regions.REGIONS[region]) for region in REGIONS}
    path = Path(reference_path) if reference_path is not None else settings.directories.ref / "rCRS.fasta"
    reference_bytes = path.read_bytes()
    reference_hash = hashlib.sha256(reference_bytes).hexdigest()
    reference = str(SeqIO.read(StringIO(reference_bytes.decode()), "fasta").seq).upper()
    for start, end in regions.values():
        if not 1 <= start <= end <= len(reference):
            msg = "Configured HV scope is outside the reference"
            raise ValueError(msg)
    occupied: set[int] = set()
    for start, end in regions.values():
        positions = set(range(start, end + 1))
        if occupied & positions:
            msg = "Configured HV scopes overlap"
            raise ValueError(msg)
        occupied.update(positions)
    region_hash = hashlib.sha256(json.dumps(regions, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    ids = [sample.sample_id for sample in cohort.samples]
    if len(ids) != len(set(ids)):
        msg = "Duplicate sample IDs in approved cohort"
        raise ValueError(msg)
    groups, sample_groups = _grouped_profiles(cohort, regions, reference)
    counts = VariantCounts()
    for profile, multiplicity in groups:
        counts.add(profile, multiplicity)
    mappings, tables = group_haplotypes(groups, sample_groups, reference_hash + region_hash, options.confidence_level)
    variants = aggregate_variants(counts, options, reference_hash, regions, reference_provider)
    haplotypes = diversity([row.count for row in tables[COMBINED]])
    for rows in tables.values():
        if rows and not math.isclose(math.fsum(row.frequency for row in rows), 1.0, abs_tol=1e-12):
            msg = "Haplotype frequencies do not sum to one"
            raise ValueError(msg)
    complete = sum(mapping.complete_hv_profile for mapping in mappings)
    if sum(row.count for row in tables[COMBINED]) != complete:
        msg = "Haplotype count invariant violated"
        raise ValueError(msg)
    # One row per sample and issue, in sample-id order, so upstream records can be found and corrected.
    data_issues, issue_counts, affected_samples = _sample_data_issues(groups, sample_groups)
    timestamp = analysis_timestamp or datetime.now(UTC)
    if timestamp.tzinfo is None:
        msg = "Analysis timestamp must be timezone-aware"
        raise ValueError(msg)
    try:
        software_version = version("mtDNA")
    except PackageNotFoundError:
        software_version = "uninstalled-source"
    databases = sorted({name for row in variants for name in row.reference_annotations})
    summary = PopulationSummary(
        analysis_timestamp=timestamp.astimezone(UTC).isoformat(),
        software_version=software_version,
        input_sha256=cohort.input_sha256,
        reference_sha256=reference_hash,
        region_configuration_sha256=region_hash,
        region_definitions=regions,
        regions_analyzed=list(REGIONS),
        eligibility_basis=cohort.eligibility_basis,
        n_samples_analyzed=len(cohort.samples),
        n_samples_complete_haplotype=complete,
        n_samples_incomplete_haplotype=len(cohort.samples) - complete,
        data_issue_counts=dict(sorted(issue_counts.items())),
        variants={
            "all": summarize_variants(variants, options),
            **{
                region: summarize_variants([row for row in variants if row.region == region], options)
                for region in REGIONS
            },
        },
        haplotypes=haplotypes,
        region_haplotypes={region: diversity([row.count for row in tables[region]]) for region in REGIONS},
        reference_not_observed_counts={
            name: sum(
                name in row.reference_annotations and row.reference_annotations[name].observed is False
                for row in variants
            )
            for name in databases
        },
        options=options,
        methods={
            "eligibility": "Supplied export is upstream-approved; no independent QC or biological deduplication",
            "callability": (
                "Declared intervals and resolved normalized variant alleles; "
                "finalized HV sequence strings are not consumed"
            ),
            "ambiguity": "Approved uppercase IUPAC exact codes; N/lowercase unresolved; no read-fraction averaging",
            "matching": "Exact normalized variants over complete configured HV1+HV2+HV3, no ignored special positions",
            "pairwise": "sum(n*(n-1))/(N*(N-1)); undefined N<2",
            "confidence_intervals": (
                "Marginal Wilson binomial intervals; no relatedness/ascertainment adjustment "
                "or match-probability bounds"
            ),
        },
    )
    logger.info(
        "Population analysis: {} complete / {} incomplete combined profiles, {} distinct variants, "
        "{} distinct complete haplotypes",
        complete,
        len(cohort.samples) - complete,
        len(variants),
        haplotypes.n_distinct_haplotypes,
    )
    if data_issues:
        logger.info("Population analysis: {} samples with upstream data issues", affected_samples)
    return PopulationResult(
        summary=summary,
        variants=variants,
        haplotypes=tables[COMBINED],
        region_haplotypes=[row for region in REGIONS for row in tables[region]],
        sample_haplotypes=mappings,
        data_issues=data_issues,
    )
