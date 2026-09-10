"""Unvalidated, REVIEW-only assessment of short usable alignment edges.

This exploratory heuristic is not an exclusion policy or a validated edge trim.
It reuses noise thresholds without changing the sustained signal-noise assessor.
"""

from __future__ import annotations

from statistics import median
from typing import Any

from src.config import get_settings
from src.tools.tracy.qc.models import BaseMetric, NoiseConfig, QCFinding
from src.tools.tracy.qc.signal_noise import extract_base_metrics
from src.tools.tracy.utils import alignment_bounds, detect_primer_type, normalize_ref_positions

_EDGE_COLUMNS = 10
_LOCAL_COLUMNS = 5
_MIN_BAD_BASES = 2


def _problematic(metric: BaseMetric, config: NoiseConfig) -> bool:
    """Require very low SNR or at least two existing per-base warning checks."""
    checks = (
        metric.snr < config.snr_threshold,
        metric.purity < config.purity_threshold,
        metric.background_ratio > config.background_threshold,
        metric.second_peak_ratio > config.second_peak_threshold,
        metric.quality is not None and metric.quality < config.quality_threshold,
        metric.signal < config.signal_threshold,
    )
    return metric.snr < config.low_snr_threshold or sum(checks) >= _MIN_BAD_BASES


def _edge_runs(
    metrics: dict[int, BaseMetric],
    edge: int,
    stop: int,
    step: int,
    config: NoiseConfig,
) -> list[list[BaseMetric]]:
    """Find consecutive measured bad runs supported by a short local window."""
    indices = list(range(edge, stop, step))[:_EDGE_COLUMNS]
    supported: set[int] = set()
    for offset in range(max(1, len(indices) - _LOCAL_COLUMNS + 1)):
        window = [metrics[index] for index in indices[offset : offset + _LOCAL_COLUMNS] if index in metrics]
        bad = [metric for metric in window if _problematic(metric, config)]
        if len(bad) >= _MIN_BAD_BASES and len(bad) / len(window) >= config.bad_base_fraction:
            supported.update(metric.alignment_index for metric in bad)

    runs: list[list[BaseMetric]] = []
    run: list[BaseMetric] = []
    for index in indices:
        if index in supported:
            run.append(metrics[index])
        else:
            if len(run) >= _MIN_BAD_BASES:
                runs.append(run)
            run = []
    if len(run) >= _MIN_BAD_BASES:
        runs.append(run)
    return runs


def assess_read_edge(data: dict[str, Any], filename: str, config: NoiseConfig) -> tuple[QCFinding, ...]:
    """Review measured bad runs in the first/last ten usable alignment columns.

    Five-column windows need at least two problematic bases and the configured
    bad-base fraction. Only consecutive problematic reference bases are reported;
    clean, missing, insertion, and deletion columns break runs. Windows never
    expand a finding to include their clean flanks. Distances are zero-based
    alignment-column distances from the usable edge, not canonical region edges.

    Metrics record the heuristic/version and unvalidated status, edge side,
    window/search sizes, measured alignment bounds, edge distances, base count,
    and medians over the reported (canonical-region-clipped) bad run. Unknown
    primer direction is not guessed. Invalid/missing evidence yields no finding,
    not an assertion that the trace is clean. Neither input nor coverage mutates.
    """
    primers = detect_primer_type(filename)
    if not (primers.is_forward or primers.is_reverse or primers.is_hv2f):
        return ()
    trace_data = dict(data)  # normalize_ref_positions replaces only this mapping's ref1pos.
    try:
        ref_start = normalize_ref_positions(trace_data)
        bounds = alignment_bounds(trace_data)
        extracted = extract_base_metrics(trace_data, ref_start=ref_start, is_reverse=primers.is_reverse)
    except (IndexError, KeyError, TypeError, ValueError):
        return ()
    if bounds.start_index is None or bounds.end_index is None:
        return ()
    metrics = {
        metric.alignment_index: metric for metric in extracted if trace_data["ref1align"][metric.alignment_index] != "-"
    }
    findings: list[QCFinding] = []
    seen: set[tuple[str, int, int]] = set()
    for side, edge, stop, step in (
        ("left", bounds.start_index, bounds.end_index + 1, 1),
        ("right", bounds.end_index, bounds.start_index - 1, -1),
    ):
        for run in _edge_runs(metrics, edge, stop, step, config):
            for region, (region_start, region_end) in get_settings().regions.REGIONS.items():
                clipped = [metric for metric in run if region_start <= metric.reference_position <= region_end]
                if not clipped:
                    continue
                start = min(metric.reference_position for metric in clipped)
                end = max(metric.reference_position for metric in clipped)
                if (region, start, end) in seen:
                    continue
                seen.add((region, start, end))
                qualities = [metric.quality for metric in clipped if metric.quality is not None]
                findings.append(
                    QCFinding(
                        code="READ_EDGE_UNCERTAINTY",
                        region=region,
                        start=start,
                        end=end,
                        traces=(filename,),
                        decision="review",
                        metrics={
                            "heuristic": "short_edge_bad_run_v1",
                            "validated": False,
                            "edge": side,
                            "edge_columns": _EDGE_COLUMNS,
                            "window_columns": _LOCAL_COLUMNS,
                            "bad_base_fraction_threshold": config.bad_base_fraction,
                            "alignment_start": min(metric.alignment_index for metric in clipped),
                            "alignment_end": max(metric.alignment_index for metric in clipped),
                            "edge_distance_min": min(abs(metric.alignment_index - edge) for metric in clipped),
                            "edge_distance_max": max(abs(metric.alignment_index - edge) for metric in clipped),
                            "valid_bases": len(clipped),
                            "median_snr": median(metric.snr for metric in clipped),
                            "median_purity": median(metric.purity for metric in clipped),
                            "median_background_ratio": median(metric.background_ratio for metric in clipped),
                            "median_second_peak_ratio": median(metric.second_peak_ratio for metric in clipped),
                            "median_signal": median(metric.signal for metric in clipped),
                            "median_quality": median(qualities) if qualities else None,
                        },
                    )
                )
    return tuple(sorted(findings, key=lambda finding: (finding.start, finding.end)))
