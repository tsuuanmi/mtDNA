"""Mutation Surveyor flagging with per-variant level classification.

MSFlag carries the flag level (SAMPLE vs VARIANT) so downstream ETL can
route flags into sample_flags vs variant_flags without guesswork.

The code that creates a flag owns its level.  Downstream consumers
use MSFlag.level to route correctly.
"""

from enum import StrEnum

from pydantic import BaseModel, Field

from src.core.models import Position


class FlagLevel(StrEnum):
    """Classification level for a pipeline flag.

    SAMPLE: flags about the sample as a whole (coverage gaps, pair status, autopass/review).
    Het is now detected via IUPAC codes in Variant.seq by
    core flag_variants().
    """

    SAMPLE = "sample"
    VARIANT = "variant"


# ---------------------------------------------------------------------------
# Region status flag constants
# ---------------------------------------------------------------------------

AUTOPASS_HV1 = "Autopass HV1"
"""Sample-level flag: HV1 region passed auto-review (no review-triggering flags)."""

AUTOPASS_HV2_AND_HV3 = "Autopass HV2 and HV3"
"""Sample-level flag: HV2-3 region passed auto-review (no review-triggering flags)."""

REVIEW_HV1 = "Review HV1"
"""Sample-level flag: HV1 region has review-triggering flags and needs manual review."""

REVIEW_HV2_AND_HV3 = "Review HV2 and HV3"
"""Sample-level flag: HV2-3 region has review-triggering flags and needs manual review."""

# ---------------------------------------------------------------------------
# Flag name sets for review-trigger classification
# ---------------------------------------------------------------------------

# Flag names that do NOT trigger review, even though they may be SAMPLE or
# VARIANT level.  Informational flags (e.g. 'HV1 polyC variant') are
# variant-level but do not trigger review.  Autopass status flags are
# SAMPLE-level but record a passing outcome, not a problem.
_NON_REVIEW_TRIGGER_NAMES: frozenset[str] = frozenset(
    {
        # Informational variant-level flags (e.g. polyC) and
        # autopass status flags do not trigger review classification.
        "HV1 polyC variant",
        AUTOPASS_HV1,
        AUTOPASS_HV2_AND_HV3,
    },
)


class MSFlag(BaseModel):
    """Structured flag with level classification for the Mutation Surveyor pipeline.

    Per-reason: one MSFlag instance per flag reason, with a level annotation
    for routing into sample_flags vs variant_flags. For Excel output, flags
    are flattened to semicolon strings; the level annotation is preserved in
    the structured pipeline only.
    """

    name: str = Field(..., description="Flag name, e.g. 'No full HV1 region', 'PHP'")
    level: FlagLevel = Field(
        ...,
        description="Whether this flag is sample-level or variant-level",
    )
    pos: Position | None = None
    ref: str | None = None
    seq: str | None = None

    @property
    def is_review_trigger(self) -> bool:
        """Whether this flag triggers a Review classification.

        Informational flags (e.g. 'HV1 polyC variant') and autopass status
        flags (e.g. 'Autopass HV1') do not trigger review. All other
        sample-level flags and most variant-level flags do trigger review.
        """
        return self.name not in _NON_REVIEW_TRIGGER_NAMES

    def to_text(self) -> str:
        """Return the flat flag string for Excel columns."""
        return self.name


def sample_flags_from(flags: list[MSFlag]) -> list[str]:
    """Extract sample-level flag names from a mixed flag list."""
    return [f.name for f in flags if f.level == FlagLevel.SAMPLE]


def variant_flags_from(flags: list[MSFlag]) -> dict[str, list[str]]:
    """Convert variant-level flags to {"pos|ref|alt": [flag_names]} dict.

    Only includes flags with a known position. Flags without a position
    are omitted (they are typically informational).
    """
    out: dict[str, list[str]] = {}
    for f in flags:
        if f.level == FlagLevel.VARIANT and f.pos is not None:
            pos_str = str(f.pos)
            ref = f.ref or "?"
            seq = f.seq or "?"
            key = f"{pos_str}|{ref}|{seq}"
            out.setdefault(key, []).append(f.name)
    return out


# ---------------------------------------------------------------------------
# Region status flag building
# ---------------------------------------------------------------------------


def build_region_status_flags(
    hv1_review_class: str,
    hv23_review_class: str,
) -> list[MSFlag]:
    """Build region-level autopass/review status flags from review class strings.

    Maps HV1 and HV2-3 review class strings to explicit status flags that
    downstream tasks can use for dataflow routing without parsing combined
    review strings.

    Args:
        hv1_review_class: HV1 Review Class string (e.g. "Auto - no review-triggering flag", "Review").
        hv23_review_class: HV2-3 Review Class string (e.g. "Auto - no flag", "Review").

    Returns:
        List of MSFlag objects with FlagLevel.SAMPLE, one per region.
    """
    flags: list[MSFlag] = []

    # HV1 region status
    if hv1_review_class and "Auto" in hv1_review_class and "Review" not in hv1_review_class:
        flags.append(MSFlag(name=AUTOPASS_HV1, level=FlagLevel.SAMPLE))
    elif hv1_review_class and hv1_review_class.strip():
        flags.append(MSFlag(name=REVIEW_HV1, level=FlagLevel.SAMPLE))

    # HV2-3 region status
    if hv23_review_class and "Auto" in hv23_review_class and "Review" not in hv23_review_class:
        flags.append(MSFlag(name=AUTOPASS_HV2_AND_HV3, level=FlagLevel.SAMPLE))
    elif hv23_review_class and hv23_review_class.strip():
        flags.append(MSFlag(name=REVIEW_HV2_AND_HV3, level=FlagLevel.SAMPLE))

    return flags
