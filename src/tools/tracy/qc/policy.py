"""Explicit action policy for Tracy evidence findings."""

from dataclasses import replace

from src.tools.tracy.qc.models import QCDecision, QCFinding

_EXCLUSION_CODES = frozenset({"SIGNAL_NOISE", "POLYC_DIRECTIONAL_UNCERTAINTY", "UNTRUSTED_OBSERVATION"})


def decide(finding: QCFinding, *, noise_mask_enabled: bool = True) -> QCFinding:
    """Promote only validated exclusions; exploratory findings remain REVIEW."""
    decision: QCDecision = "review"
    if finding.code in _EXCLUSION_CODES:
        decision = "exclude"
    if finding.code == "SIGNAL_NOISE" and not noise_mask_enabled:
        decision = "review"
    return replace(finding, decision=decision)
