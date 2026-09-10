"""Explicit orchestration of optional per-trace evidence assessors."""

from dataclasses import replace
from typing import Any

from loguru import logger

from src.core.models import Sample
from src.core.variants import pos_base
from src.tools.tracy.evidence import TraceEvidence
from src.tools.tracy.qc.coverage import apply_findings
from src.tools.tracy.qc.models import ExcludedVariant, QCConfig, QCFinding, TraceQCResult
from src.tools.tracy.qc.policy import decide
from src.tools.tracy.qc.polyc import assess_directional_polyc, assess_homopolymers
from src.tools.tracy.qc.read_edge import assess_read_edge


def evaluate_trace(
    sample: Sample,
    evidence: TraceEvidence,
    result: TraceQCResult,
    config: QCConfig,
    data: dict[str, Any],
) -> tuple[Sample, TraceQCResult]:
    """Run enabled assessors, apply policy, and retain source-level diagnostics."""
    raw_coverage = sample.intervals
    if not config.enabled:
        return sample, replace(result, raw_coverage=raw_coverage, trusted_coverage=raw_coverage)
    findings: list[QCFinding] = []
    filename = result.filename or evidence.filename
    if config.signal_noise_enabled:
        findings.extend(
            QCFinding(
                code="SIGNAL_NOISE",
                region=None,
                start=r.start,
                end=r.end,
                traces=(filename,),
                decision="review",
                metrics=r.to_dict(),
            )
            for r in result.ranges
        )
    if config.polyc_enabled:
        findings.extend(
            assess_directional_polyc(
                raw_coverage or {},
                filename,
                polyc_created=evidence.has_16189_t_c,
            )
        )
        findings.extend(assess_homopolymers(result))
    if config.read_edge_enabled:
        findings.extend(assess_read_edge(data, filename, config.noise))
    # ETL records quality-rejected observations separately from normalization.
    # Rejecting their evidence must not turn their covered base into reference.
    for observation in evidence.observations:
        if observation.get("disposition") not in {"untrusted", "invalid"}:
            continue
        from_positions = observation.get("source_positions", [])
        for position in from_positions:
            anchor = pos_base(position)
            findings.append(
                QCFinding(
                    code="UNTRUSTED_OBSERVATION",
                    region=None,
                    start=anchor,
                    end=anchor + max(1, len(observation.get("ref", ""))) - 1,
                    traces=(filename,),
                    decision="review",
                    metrics={"reason": observation.get("reason"), "candidate_id": observation.get("candidate_id")},
                )
            )
    decisions = tuple(decide(finding, noise_mask_enabled=config.noise.mask_enabled) for finding in findings)
    masked = apply_findings(sample, decisions)
    observations = [dict(observation) for observation in evidence.observations]
    for excluded in masked.excluded:
        variant, finding = excluded.variant, excluded.finding
        observations.append(
            {
                "pos": variant.pos,
                "ref": variant.ref,
                "seq": variant.seq,
                "disposition": "excluded",
                "reason": finding.code,
                "finding": finding.to_dict(),
            }
        )
        logger.warning(
            "Sample {} trace {}: Exclude variant {} {}>{} ({})",
            sample.sample_id,
            filename,
            variant.pos,
            variant.ref,
            variant.seq,
            finding.code,
        )
    noise_excluded = tuple(
        ExcludedVariant(v.variant.pos, v.variant.ref, v.variant.seq, noise_range)
        for v in masked.excluded
        for noise_range in result.ranges
        if v.finding.code == "SIGNAL_NOISE"
        and noise_range.start == v.finding.start
        and noise_range.end == v.finding.end
    )
    report = replace(
        result,
        findings=decisions,
        raw_coverage=raw_coverage,
        trusted_coverage=masked.sample.intervals,
        excluded_variants=noise_excluded,
        excluded_observations=tuple(observations),
    )
    return masked.sample, report
