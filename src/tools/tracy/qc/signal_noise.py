"""Window-based trace noise quality control for Tracy decompose JSON.

This module calculates peak-separation metrics and reports likely-noisy ranges.
It does not call or filter variants and has no dependency on Tracy ETL models.
"""

from __future__ import annotations

import json
from math import isfinite
from statistics import median
from typing import TYPE_CHECKING, Any

from src.tools.tracy.qc.models import (
    BaseMetric,
    NoiseConfig,
    NoiseRange,
    TraceQCResult,
    TraceStatus,
    WindowResult,
    WindowStatus,
    error_result,
    unavailable_result,
)
from src.tools.tracy.utils import (
    alignment_bounds,
    alignment_reference_position,
    detect_primer_type,
    extract_peak_data,
    normalize_ref_positions,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

_MIN_MODERATE_NOISE_CHECKS = 2


def _percentile(values: list[float], percentile: float) -> float:
    """Return a linearly interpolated percentile without a numeric dependency."""
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (len(ordered) - 1) * percentile
    lower = int(rank)
    upper = min(lower + 1, len(ordered) - 1)
    weight = rank - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def _fraction(metrics: list[BaseMetric], predicate: Callable[[BaseMetric], bool]) -> float:
    """Return the fraction of metrics matching a predicate."""
    return sum(bool(predicate(metric)) for metric in metrics) / len(metrics)


def _quality_fraction(metrics: list[BaseMetric], threshold: float) -> float:
    """Return the low-quality fraction among bases with quality values."""
    quality_metrics = [metric for metric in metrics if metric.quality is not None]
    if not quality_metrics:
        return 0.0
    return sum(metric.quality < threshold for metric in quality_metrics if metric.quality is not None) / len(
        quality_metrics,
    )


def _classify_window(
    *,
    median_snr: float,
    low_end_snr: float,
    median_purity: float,
    low_snr_fraction: float,
    low_purity_fraction: float,
    high_background_fraction: float,
    high_second_peak_fraction: float,
    low_quality_fraction: float,
    weak_signal_fraction: float,
    config: NoiseConfig,
) -> WindowStatus:
    """Classify a window using sustained and independent noise evidence."""
    strong_low_snr = low_snr_fraction >= config.bad_base_fraction and low_end_snr < config.low_snr_threshold
    moderate_checks = (
        median_snr < config.snr_threshold,
        median_purity < config.purity_threshold,
        high_background_fraction >= config.bad_base_fraction,
        high_second_peak_fraction >= config.bad_base_fraction,
        low_quality_fraction >= config.bad_base_fraction,
        weak_signal_fraction >= config.bad_base_fraction,
    )
    if strong_low_snr or sum(moderate_checks) >= _MIN_MODERATE_NOISE_CHECKS:
        return "likely_noisy"

    suspicious_checks = (
        low_snr_fraction >= config.suspicious_base_fraction,
        low_purity_fraction >= config.suspicious_base_fraction,
        high_background_fraction >= config.suspicious_base_fraction,
        high_second_peak_fraction >= config.suspicious_base_fraction,
        low_quality_fraction >= config.suspicious_base_fraction,
        weak_signal_fraction >= config.suspicious_base_fraction,
    )
    return "suspicious" if any(suspicious_checks) else "clean"


def _window_result(metrics: list[BaseMetric], config: NoiseConfig) -> WindowResult:
    """Aggregate and classify one non-empty metrics window."""
    snr_values = [metric.snr for metric in metrics]
    purity_values = [metric.purity for metric in metrics]
    quality_values = [metric.quality for metric in metrics if metric.quality is not None]
    median_snr = median(snr_values)
    low_end_snr = _percentile(snr_values, 0.10)
    median_purity = median(purity_values)
    low_snr_fraction = _fraction(metrics, lambda metric: metric.snr < config.snr_threshold)
    low_purity_fraction = _fraction(metrics, lambda metric: metric.purity < config.purity_threshold)
    high_background_fraction = _fraction(
        metrics,
        lambda metric: metric.background_ratio > config.background_threshold,
    )
    high_second_peak_fraction = _fraction(
        metrics,
        lambda metric: metric.second_peak_ratio > config.second_peak_threshold,
    )
    low_quality_fraction = _quality_fraction(metrics, config.quality_threshold)
    weak_signal_fraction = _fraction(metrics, lambda metric: metric.signal < config.signal_threshold)
    classification = _classify_window(
        median_snr=median_snr,
        low_end_snr=low_end_snr,
        median_purity=median_purity,
        low_snr_fraction=low_snr_fraction,
        low_purity_fraction=low_purity_fraction,
        high_background_fraction=high_background_fraction,
        high_second_peak_fraction=high_second_peak_fraction,
        low_quality_fraction=low_quality_fraction,
        weak_signal_fraction=weak_signal_fraction,
        config=config,
    )
    positions = [metric.reference_position for metric in metrics]
    return WindowResult(
        alignment_start=metrics[0].alignment_index,
        alignment_end=metrics[-1].alignment_index,
        start=min(positions),
        end=max(positions),
        valid_bases=len(metrics),
        median_snr=median_snr,
        low_end_snr=low_end_snr,
        median_purity=median_purity,
        median_background_ratio=median(metric.background_ratio for metric in metrics),
        median_second_peak_ratio=median(metric.second_peak_ratio for metric in metrics),
        median_quality=median(quality_values) if quality_values else None,
        low_snr_fraction=low_snr_fraction,
        low_purity_fraction=low_purity_fraction,
        high_background_fraction=high_background_fraction,
        high_second_peak_fraction=high_second_peak_fraction,
        low_quality_fraction=low_quality_fraction,
        weak_signal_fraction=weak_signal_fraction,
        classification=classification,
    )


def _build_windows(metrics: list[BaseMetric], config: NoiseConfig) -> tuple[WindowResult, ...]:
    """Build overlapping windows with enough valid called bases."""
    windows: list[WindowResult] = []
    for start in range(0, len(metrics), config.window_step):
        window_metrics = metrics[start : start + config.window_size]
        if len(window_metrics) < config.min_valid_bases:
            break
        windows.append(_window_result(window_metrics, config))
    return tuple(windows)


def _range_from_windows(windows: list[WindowResult]) -> NoiseRange:
    """Create one merged range from overlapping likely-noisy windows."""
    return NoiseRange(
        start=min(window.start for window in windows),
        end=max(window.end for window in windows),
        alignment_start=min(window.alignment_start for window in windows),
        alignment_end=max(window.alignment_end for window in windows),
        supporting_windows=len(windows),
        median_snr=median(window.median_snr for window in windows),
        low_end_snr=min(window.low_end_snr for window in windows),
        median_purity=median(window.median_purity for window in windows),
    )


def _merge_noisy_windows(
    windows: tuple[WindowResult, ...],
    min_supporting_windows: int,
) -> tuple[NoiseRange, ...]:
    """Merge noisy windows when a region has enough independent support."""
    noisy = [window for window in windows if window.classification == "likely_noisy"]
    if not noisy:
        return ()

    groups: list[list[WindowResult]] = []
    current = [noisy[0]]
    for window in noisy[1:]:
        if window.alignment_start <= current[-1].alignment_end + 1:
            current.append(window)
        else:
            groups.append(current)
            current = [window]
    groups.append(current)
    return tuple(_range_from_windows(group) for group in groups if len(group) >= min_supporting_windows)


def extract_base_metrics(
    data: dict[str, Any],
    *,
    ref_start: int,
    is_reverse: bool,
) -> list[BaseMetric]:
    """Calculate metrics for called bases in the trimmed alignment."""
    bounds = alignment_bounds(data)
    if bounds.start_index is None or bounds.end_index is None:
        return []

    metrics: list[BaseMetric] = []
    deletions = 0
    for alignment_index in range(bounds.start_index, bounds.end_index + 1):
        ref_base = data["ref1align"][alignment_index]
        alt_base = data["alt1align"][alignment_index]
        peak_index = alignment_index - deletions
        if ref_base != "-" and alt_base == "-":
            deletions += 1
            continue

        peaks = extract_peak_data(data, peak_index)
        if any(value is None for value in peaks.channels):
            continue
        channels = sorted((value for value in peaks.channels if value is not None), reverse=True)
        if any(not isfinite(value) or value < 0 for value in channels):
            continue
        signal = channels[0]
        noise = median(channels[1:])
        total = sum(channels)
        reference_position = alignment_reference_position(
            alignment_index,
            ref_start=ref_start,
            ref_align=data["ref1align"],
            is_reverse=is_reverse,
        )
        metrics.append(
            BaseMetric(
                alignment_index=alignment_index,
                reference_position=reference_position,
                signal=signal,
                noise=noise,
                snr=signal / max(noise, 1.0),
                purity=signal / total if total else 0.0,
                background_ratio=noise / signal if signal else 1.0,
                second_peak_ratio=channels[1] / signal if signal else 1.0,
                quality=peaks.quality
                if peaks.quality is not None and isfinite(peaks.quality) and peaks.quality >= 0
                else None,
            ),
        )
    return metrics


def analyze_trace_data(
    data: dict[str, Any],
    *,
    sample_id: str,
    filename: str,
    config: NoiseConfig | None = None,
) -> TraceQCResult:
    """Analyze one parsed Tracy JSON mapping without calling variants."""
    config = config or NoiseConfig.from_settings()
    required = {"ref1pos", "ref1align", "alt1align", "basecallPos", "peakA", "peakC", "peakG", "peakT"}
    missing = sorted(required.difference(data))
    if missing:
        return unavailable_result(sample_id, f"Missing Tracy fields: {', '.join(missing)}", filename)

    trace_data = dict(data)
    try:
        ref_start = normalize_ref_positions(trace_data)
        primer_info = detect_primer_type(filename)
        metrics = extract_base_metrics(trace_data, ref_start=ref_start, is_reverse=primer_info.is_reverse)
    except (IndexError, KeyError, TypeError, ValueError) as error:
        return error_result(sample_id, str(error), filename)

    primer = next((name for name in ("HV1F", "HV1R", "HV2F", "HV2R", "HV3F", "HV3R") if name in filename), None)
    strand = (
        "reverse" if primer_info.is_reverse else "forward" if primer_info.is_forward or primer_info.is_hv2f else None
    )
    if not metrics:
        result = unavailable_result(sample_id, "No valid called bases were available for QC", filename)
        return TraceQCResult(
            result.sample_id,
            result.filename,
            result.status,
            primer,
            strand,
            result.bases_evaluated,
            reason=result.reason,
        )

    windows = _build_windows(metrics, config)
    if not windows:
        return TraceQCResult(
            sample_id,
            filename,
            "unavailable",
            primer,
            strand,
            len(metrics),
            reason=f"Fewer than {config.min_valid_bases} valid bases were available",
        )
    ranges = _merge_noisy_windows(windows, config.min_supporting_windows)
    status: TraceStatus
    if ranges:
        status = "likely_noisy"
    elif any(window.classification != "clean" for window in windows):
        status = "suspicious"
    else:
        status = "clean"
    return TraceQCResult(sample_id, filename, status, primer, strand, len(metrics), windows, ranges)


def analyze_trace_file(
    path: Path,
    *,
    sample_id: str,
    config: NoiseConfig | None = None,
) -> TraceQCResult:
    """Load and analyze one Tracy decompose JSON file."""
    try:
        with path.open(encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        return error_result(sample_id, str(error), path.name)
    if not isinstance(data, dict):
        return error_result(sample_id, "Tracy JSON root must be an object", path.name)
    return analyze_trace_data(data, sample_id=sample_id, filename=path.name, config=config)
