"""Shared Tracy QC evidence, decisions, and trace result models."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

from src.config import get_settings

TraceStatus = Literal["clean", "suspicious", "likely_noisy", "unavailable", "error"]
WindowStatus = Literal["clean", "suspicious", "likely_noisy"]
QCDecision = Literal["accept", "review", "exclude"]


@dataclass(frozen=True)
class QCFinding:
    """Trace-local finding over an inclusive canonical reference interval."""

    code: str
    region: str | None
    start: int
    end: int
    traces: tuple[str, ...]
    decision: QCDecision
    metrics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready finding with source-trace provenance."""
        return {**asdict(self), "traces": list(self.traces)}


@dataclass(frozen=True)
class NoiseConfig:
    """Settings for noise masking, base metrics, and overlapping windows."""

    mask_enabled: bool
    window_size: int
    window_step: int
    min_valid_bases: int
    min_supporting_windows: int
    snr_threshold: float
    low_snr_threshold: float
    purity_threshold: float
    background_threshold: float
    second_peak_threshold: float
    quality_threshold: float
    signal_threshold: float
    bad_base_fraction: float
    suspicious_base_fraction: float

    @classmethod
    def from_settings(cls) -> NoiseConfig:
        """Build QC configuration from Tracy application settings."""
        settings = get_settings().tracy
        return cls(
            mask_enabled=settings.noise_mask_enabled,
            window_size=settings.noise_window_size,
            window_step=settings.noise_window_step,
            min_valid_bases=settings.noise_min_valid_bases,
            min_supporting_windows=settings.noise_min_supporting_windows,
            snr_threshold=settings.noise_snr_threshold,
            low_snr_threshold=settings.noise_low_snr_threshold,
            purity_threshold=settings.noise_purity_threshold,
            background_threshold=settings.noise_background_threshold,
            second_peak_threshold=settings.noise_second_peak_threshold,
            quality_threshold=settings.noise_quality_threshold,
            signal_threshold=settings.noise_signal_threshold,
            bad_base_fraction=settings.noise_bad_base_fraction,
            suspicious_base_fraction=settings.noise_suspicious_base_fraction,
        )


@dataclass(frozen=True)
class QCConfig:
    """Independent switches for the optional Tracy QC assessors."""

    enabled: bool
    signal_noise_enabled: bool
    polyc_enabled: bool
    read_edge_enabled: bool
    noise: NoiseConfig

    @classmethod
    def from_settings(cls) -> QCConfig:
        """Build the QC engine configuration from Tracy settings."""
        settings = get_settings().tracy
        return cls(
            enabled=settings.qc_enabled,
            signal_noise_enabled=settings.signal_noise_enabled,
            polyc_enabled=settings.polyc_enabled,
            read_edge_enabled=settings.read_edge_enabled,
            noise=NoiseConfig.from_settings(),
        )


@dataclass(frozen=True)
class BaseMetric:
    """Peak-separation measurements for one called base."""

    alignment_index: int
    reference_position: int
    signal: float
    noise: float
    snr: float
    purity: float
    background_ratio: float
    second_peak_ratio: float
    quality: float | None


@dataclass(frozen=True)
class WindowResult:
    """Aggregated measurements and classification for one trace window."""

    alignment_start: int
    alignment_end: int
    start: int
    end: int
    valid_bases: int
    median_snr: float
    low_end_snr: float
    median_purity: float
    median_background_ratio: float
    median_second_peak_ratio: float
    median_quality: float | None
    low_snr_fraction: float
    low_purity_fraction: float
    high_background_fraction: float
    high_second_peak_fraction: float
    low_quality_fraction: float
    weak_signal_fraction: float
    classification: WindowStatus

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready window record."""
        return asdict(self)


@dataclass(frozen=True)
class NoiseRange:
    """Merged reference and alignment range supported by noisy windows."""

    start: int
    end: int
    alignment_start: int
    alignment_end: int
    supporting_windows: int
    median_snr: float
    low_end_snr: float
    median_purity: float
    classification: Literal["likely_noisy"] = "likely_noisy"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready range record."""
        return asdict(self)


@dataclass(frozen=True)
class ExcludedVariant:
    """Variant removed because it overlaps a likely-noisy range."""

    pos: int | str
    ref: str
    seq: str
    noise_range: NoiseRange
    reason: str = "overlaps likely-noisy range"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready exclusion record."""
        return {
            "pos": self.pos,
            "ref": self.ref,
            "seq": self.seq,
            "reason": self.reason,
            "range": self.noise_range.to_dict(),
        }


@dataclass(frozen=True)
class TraceQCResult:
    """QC result for one Tracy trace or one missing trace output."""

    sample_id: str
    filename: str | None
    status: TraceStatus
    primer: str | None
    strand: str | None
    bases_evaluated: int
    windows: tuple[WindowResult, ...] = ()
    ranges: tuple[NoiseRange, ...] = ()
    excluded_variants: tuple[ExcludedVariant, ...] = ()
    reason: str | None = None
    findings: tuple[QCFinding, ...] = ()
    raw_coverage: dict[str, list[list[int]]] | None = None
    trusted_coverage: dict[str, list[list[int]]] | None = None
    excluded_observations: tuple[dict[str, Any], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a compact JSON-ready trace record."""
        record: dict[str, Any] = {
            "sample_id": self.sample_id,
            "filename": self.filename,
            "status": self.status,
            "primer": self.primer,
            "strand": self.strand,
            "bases_evaluated": self.bases_evaluated,
            "windows": [window.to_dict() for window in self.windows],
            "ranges": [noise_range.to_dict() for noise_range in self.ranges],
            "excluded_variants": [variant.to_dict() for variant in self.excluded_variants],
        }
        if self.reason is not None:
            record["reason"] = self.reason
        if self.findings:
            record["findings"] = [finding.to_dict() for finding in self.findings]
        if self.raw_coverage is not None:
            record["raw_coverage"] = self.raw_coverage
        if self.trusted_coverage is not None:
            record["trusted_coverage"] = self.trusted_coverage
        if self.excluded_observations:
            record["excluded_observations"] = list(self.excluded_observations)
        return record


def unavailable_result(sample_id: str, reason: str, filename: str | None = None) -> TraceQCResult:
    """Create a result for a trace that could not provide QC measurements."""
    return TraceQCResult(sample_id, filename, "unavailable", None, None, 0, reason=reason)


def error_result(sample_id: str, reason: str, filename: str | None = None) -> TraceQCResult:
    """Create a result for an invalid or unreadable trace."""
    return TraceQCResult(sample_id, filename, "error", None, None, 0, reason=reason)
