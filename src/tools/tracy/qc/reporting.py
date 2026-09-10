"""Serialize Tracy QC results into per-sample reports."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from src.tools.tracy.qc.models import NoiseConfig, TraceQCResult


def build_qc_report(results: list[TraceQCResult], config: NoiseConfig | None = None) -> dict[str, Any]:
    """Build the batch-level Tracy QC report."""
    config = config or NoiseConfig.from_settings()
    ordered = sorted(results, key=lambda result: (result.sample_id, result.filename or ""))
    counts = {
        status: sum(result.status == status for result in ordered)
        for status in (
            "clean",
            "suspicious",
            "likely_noisy",
            "unavailable",
            "error",
        )
    }
    return {
        "tool": "tracy",
        "total": len(ordered),
        "analyzed": counts["clean"] + counts["suspicious"] + counts["likely_noisy"],
        **counts,
        "settings": asdict(config),
        "traces": [result.to_dict() for result in ordered],
    }
