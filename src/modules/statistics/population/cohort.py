"""Read an explicitly approved export, never infer PASS or regenerate missing calls."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from loguru import logger
from pydantic import TypeAdapter

from src.core.models import Sample
from src.modules.statistics.population.models import REGIONS, PopulationCohort

_SAMPLES_ADAPTER = TypeAdapter(list[Sample])


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            msg = "Duplicate JSON object key in approved export"
            raise ValueError(msg)
        result[key] = value
    return result


def _validate_coverage_fields(raw: dict[str, Any]) -> None:
    intervals = raw.get("intervals")
    if intervals is None:
        return
    if not isinstance(intervals, dict) or not set(intervals) <= set(REGIONS):
        msg = "Intervals must be keyed by HV1, HV2 or HV3"
        raise ValueError(msg)
    for ranges in intervals.values():
        if not isinstance(ranges, list):
            msg = "Region intervals must be lists of endpoint pairs"
            raise TypeError(msg)
        for interval in ranges:
            if not isinstance(interval, list) or len(interval) != len(("start", "end")):
                msg = "Intervals require exactly two integer endpoints"
                raise ValueError(msg)
            if any(type(endpoint) is not int for endpoint in interval):
                msg = "Interval endpoints must be integers without coercion"
                raise TypeError(msg)


def _validated_entry(sample_id: str, raw: Any) -> dict[str, Any]:  # noqa: ANN401
    if not sample_id or sample_id != sample_id.strip() or not isinstance(raw, dict):
        msg = "Expected nonempty sample IDs mapped to canonical Sample entries"
        raise ValueError(msg)
    if "sample_id" in raw and raw["sample_id"] != sample_id:
        msg = "Embedded sample ID disagrees with export key"
        raise ValueError(msg)
    _validate_coverage_fields(raw)
    variants = raw.get("variants")
    if isinstance(variants, dict):
        keys = set(variants)
        if not (keys <= {"snps", "insertions", "deletions"} or keys <= set(REGIONS)):
            msg = "Unsupported grouped variant schema"
            raise ValueError(msg)
        if any(not isinstance(group, list) for group in variants.values()):
            msg = "Variant groups must be lists"
            raise ValueError(msg)
        variants = [variant for group in variants.values() for variant in group]
    if not isinstance(variants, list):
        msg = "Each export entry must declare its normalized variants"
        raise TypeError(msg)
    # Trace evidence and personal metadata are neither needed nor loaded into analysis.
    canonical = []
    for variant in variants:
        if not isinstance(variant, dict) or not {"pos", "ref", "seq"} <= variant.keys():
            msg = "Expected canonical pos/ref/seq variant entries"
            raise ValueError(msg)
        if isinstance(variant["pos"], bool) or not isinstance(variant["pos"], int | str):
            msg = "Variant positions must be integer or string, never float/bool"
            raise TypeError(msg)
        canonical.append({key: variant[key] for key in ("pos", "ref", "seq")})
    # Trace evidence, HV strings and personal metadata are not loaded: population
    # statistics consume declared intervals and normalized variant alleles only.
    return {
        "variants": canonical,
        **({"intervals": raw["intervals"]} if "intervals" in raw else {}),
    }


def load_cohort(path: str | Path) -> PopulationCohort:
    """Reload an approved snapshot, preserving decimal position tokens losslessly.

    JSON decimals are decoded as their source strings, never binary floats, so
    insertion indices such as 309.10 retain their exact digits. Interval endpoints
    still require JSON integers and reject these decimal strings. Validated
    entries are model-checked in one bulk pass; sample construction performs the
    same pydantic validation as per-sample loading.
    """
    content = Path(path).read_bytes()
    data = json.loads(content, object_pairs_hook=_unique_object, parse_float=str)
    if not isinstance(data, dict):
        msg = "Approved cohort must be a JSON object keyed by sample ID"
        raise TypeError(msg)
    payload = [
        {"sample_id": sample_id, "source_tool": "unknown", **_validated_entry(sample_id, raw)}
        for sample_id, raw in data.items()
    ]
    del data  # free the raw export before model construction
    input_sha256 = hashlib.sha256(content).hexdigest()
    samples = _SAMPLES_ADAPTER.validate_python(payload)
    logger.info("Population cohort: loaded {} approved samples from {} (sha256 {})", len(samples), path, input_sha256)
    return PopulationCohort(samples=samples, input_sha256=input_sha256)
