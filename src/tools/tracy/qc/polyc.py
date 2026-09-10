"""Validated directional exclusions and exploratory homopolymer evidence."""

from src.core.polyc import directional_hv1_polyc_suppression_reason, directional_hv2_polyc_exclusion_reason
from src.tools.tracy.qc.models import QCFinding, TraceQCResult
from src.tools.tracy.utils import detect_primer_type

_POLYC_REVIEW_REGION = (303, 315)
_POLYC_REVIEW_PURITY = 0.75
_PHASE_REVIEW_REGION = (574, 576)
_PHASE_SECOND_PEAK_FRACTION = 0.10


def assess_directional_polyc(
    intervals: dict[str, list[list[int]]],
    filename: str,
    *,
    polyc_created: bool,
) -> tuple[QCFinding, ...]:
    """Classify aligned positions using the authoritative shared core policy."""
    primers = detect_primer_type(filename)
    findings: list[QCFinding] = []
    for region, spans in intervals.items():
        for start, end in spans:
            run_start: int | None = None
            for position in range(start, end + 2):
                reasons = (
                    ()
                    if position > end
                    else (
                        directional_hv2_polyc_exclusion_reason(position, "forward" if primers.is_hv2f_hv3f else None),
                        directional_hv2_polyc_exclusion_reason(position, "reverse" if primers.is_hv2r_hv3r else None),
                        directional_hv1_polyc_suppression_reason(
                            position,
                            "forward" if primers.is_hv1f else None,
                            polyc_created=polyc_created,
                        ),
                        directional_hv1_polyc_suppression_reason(
                            position,
                            "reverse" if primers.is_hv1r else None,
                            polyc_created=polyc_created,
                        ),
                    )
                )
                excluded = any(reason is not None for reason in reasons)
                if excluded and run_start is None:
                    run_start = position
                elif not excluded and run_start is not None:
                    findings.append(
                        QCFinding(
                            code="POLYC_DIRECTIONAL_UNCERTAINTY",
                            region=region,
                            start=run_start,
                            end=position - 1,
                            traces=(filename,),
                            decision="review",
                            metrics={"polyc_created": polyc_created},
                        )
                    )
                    run_start = None
    return tuple(findings)


def assess_homopolymers(result: TraceQCResult) -> tuple[QCFinding, ...]:
    """Report same-batch exploratory signals without declaring missing coverage."""
    if result.primer not in {"HV2R", "HV3R"}:
        return ()
    findings: list[QCFinding] = []
    for code, region, interval, metric, threshold in (
        ("HOMOPOLYMER_ALIGNMENT_UNCERTAINTY", "HV2", _POLYC_REVIEW_REGION, "median_purity", _POLYC_REVIEW_PURITY),
        (
            "POST_HOMOPOLYMER_PHASE_UNCERTAINTY",
            "HV3",
            _PHASE_REVIEW_REGION,
            "high_second_peak_fraction",
            _PHASE_SECOND_PEAK_FRACTION,
        ),
    ):
        windows = [w for w in result.windows if w.start <= interval[1] and w.end >= interval[0]]
        matching = [
            w
            for w in windows
            if (w.median_purity < threshold if metric == "median_purity" else w.high_second_peak_fraction >= threshold)
        ]
        if matching:
            findings.append(
                QCFinding(
                    code=code,
                    region=region,
                    start=interval[0],
                    end=interval[1],
                    traces=(result.filename or result.primer,),
                    decision="review",
                    metrics={"threshold": threshold, "metric": metric, "windows": [w.to_dict() for w in matching]},
                )
            )
    return tuple(findings)
