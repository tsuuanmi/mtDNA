"""Neutral, trace-local ETL evidence contract; no dependency on optional QC."""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Literal

Disposition = Literal["normalized", "untrusted", "out_of_scope", "accepted", "invalid"]


@dataclass
class TraceEvidence:
    """Mutable output accumulator, reset for each ``etl.process`` invocation.

    Raw candidates are deep snapshots after genomic positioning, before
    nomenclature transforms. IDs are trace-local stable alignment-column IDs.
    ``raw_intervals`` is region-keyed aligned coverage, before synthetic extension.
    Observations retain original/source coordinates even when nomenclature moves
    a call. Derived records have no measured quality/peaks in this contract;
    their legacy placeholder values in canonical Variants are not QC evidence.
    """

    raw_intervals: dict[str, list[list[int]]] = field(default_factory=dict)
    filename: str = ""
    has_16189_t_c: bool = False
    raw_candidates: list[dict[str, Any]] = field(default_factory=list)
    observations: list[dict[str, Any]] = field(default_factory=list)


def observe(evidence: TraceEvidence | None, variant: dict[str, Any], disposition: Disposition, reason: str) -> None:
    """Append an independent audit event without promoting synthetic metrics."""
    if evidence is None:
        return
    record = deepcopy(variant)
    record.setdefault("source_positions", [record["pos"]])
    record.update(disposition=disposition, reason=reason)
    if record.get("evidence_kind") == "derived":
        record.pop("quality", None)
        record.pop("peaks", None)
    evidence.observations.append(record)
