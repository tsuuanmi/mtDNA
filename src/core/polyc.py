"""Shared directional polyC policy helpers.

Tools call variants before using these helpers to classify directional polyC
artifacts and coverage exclusions. The helpers do not trim sequence data or mutate
variants, so callers can log or retain raw candidate evidence separately from
the trusted consensus input.
"""

from typing import Literal

from src.core.models import Position
from src.core.variants import pos_base

DirectionalTrace = Literal["forward", "reverse"]

HV1_POLYC_TRIGGER_POSITION = 16189
HV1_DIRECTIONAL_POST_POLYC_REASON = "HV1 polyC region variant with 16189T>C present"
_HV1_DIRECTIONAL_TRACE_LABELS: dict[DirectionalTrace, str] = {
    "forward": "HV1F primer - forward",
    "reverse": "HV1R primer - reverse",
}

HV2_POLYC_REGION: tuple[int, int] = (303, 315)
HV2_DIRECTIONAL_EXCLUSION_REASON = "HV2 directional polyC variant"

__all__ = [
    "HV1_DIRECTIONAL_POST_POLYC_REASON",
    "HV1_POLYC_TRIGGER_POSITION",
    "HV2_DIRECTIONAL_EXCLUSION_REASON",
    "HV2_POLYC_REGION",
    "DirectionalTrace",
    "directional_hv1_polyc_suppression_reason",
    "directional_hv2_polyc_exclusion_reason",
    "is_hv1_polyc_created",
    "is_hv2_directionally_excluded",
]


def is_hv1_polyc_created(position: Position, seq: str) -> bool:
    """Return whether a called C at 16189 creates the HV1 polyC condition.

    This deliberately preserves the caller's existing anchor-plus-called-base
    condition without imposing an additional reference-base requirement.
    """
    return pos_base(position) == HV1_POLYC_TRIGGER_POSITION and seq == "C"


def directional_hv1_polyc_suppression_reason(
    position: Position,
    direction: DirectionalTrace | None,
    *,
    polyc_created: bool,
) -> str | None:
    """Return the shared reason for an HV1 directional polyC artifact.

    HV1F candidates after, and HV1R candidates before, the 16189 C trigger
    are suppressed. The trigger position itself and unknown provenance remain
    usable.
    """
    if not polyc_created or direction is None:
        return None

    pos_int = pos_base(position)
    is_post_polyc = (direction == "forward" and pos_int > HV1_POLYC_TRIGGER_POSITION) or (
        direction == "reverse" and pos_int < HV1_POLYC_TRIGGER_POSITION
    )
    if not is_post_polyc:
        return None

    trace_label = _HV1_DIRECTIONAL_TRACE_LABELS[direction]
    return f"{HV1_DIRECTIONAL_POST_POLYC_REASON} ({trace_label} sequencing artifact)"


def is_hv2_directionally_excluded(
    position: Position,
    direction: DirectionalTrace,
) -> bool:
    """Return whether ``position`` is excluded from an HV2 trace by direction.

    Forward traces exclude every anchor from 304 onward; reverse traces exclude
    every anchor through 315. Decimal insertion positions use their integer
    anchor.
    """
    pos_int = pos_base(position)
    return (direction == "forward" and pos_int > HV2_POLYC_REGION[0]) or (
        direction == "reverse" and pos_int <= HV2_POLYC_REGION[1]
    )


def directional_hv2_polyc_exclusion_reason(
    position: Position,
    direction: DirectionalTrace | None,
) -> str | None:
    """Return the shared reason for an HV2 directional exclusion.

    ``None`` direction represents unknown provenance and is intentionally not
    excluded. The policy decision and its reason stay centralized here so
    callers do not duplicate directional coverage logic.
    """
    if direction is None or not is_hv2_directionally_excluded(position, direction):
        return None
    return f"{HV2_DIRECTIONAL_EXCLUSION_REASON} ({direction} trace)"
