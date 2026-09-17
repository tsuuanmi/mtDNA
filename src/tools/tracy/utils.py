"""Tracy-specific utilities — primer detection, heteroplasmy, peak validation.

Pure functions only: no file I/O, no side effects, no Sample/Variant dependency.
This module is the single source of truth for Tracy filename parsing, primer type
detection, and variant filtering helpers.
"""

from typing import Any, NamedTuple

from loguru import logger

from src.config import get_settings
from src.core.variants import nucleotide_to_iupac, pos_base

# ---------------------------------------------------------------------------
# mtDNA position constants
# ---------------------------------------------------------------------------

POS_73 = 73
POS_248 = 248
POS_250 = 250
POS_300 = 300
POS_309 = 309
POS_310 = 310
POS_315 = 315
POS_316 = 316
POS_320 = 320
POS_456 = 456
POS_459 = 459
POS_513 = 513
POS_514 = 514
POS_515 = 515
POS_523 = 523
POS_524 = 524
POS_525 = 525
POS_16000 = 16000
POS_16189 = 16189
POS_16193 = 16193

# Number of peak channels (A, C, G, T)
NUM_PEAK_CHANNELS = 4

# ---------------------------------------------------------------------------
# Primer type detection
# ---------------------------------------------------------------------------


class PrimerTypeInfo(NamedTuple):
    """Primer type information detected from a trace filename.

    Adds ``is_forward`` and ``is_hv2f`` fields for convenience.
    """

    is_forward: bool
    is_reverse: bool
    is_hv2f: bool
    is_hv2f_hv3f: bool
    is_hv2r_hv3r: bool
    is_hv1f: bool
    is_hv1r: bool


def detect_primer_type(filename: str) -> PrimerTypeInfo:
    """Determine primer types from a Tracy decompose output filename.

    Returns all primer-type flags in a single ``PrimerTypeInfo`` NamedTuple.
    """
    return PrimerTypeInfo(
        is_forward=any(p in filename for p in ["HV1F", "HV3F"]),
        is_reverse=any(p in filename for p in ["HV1R", "HV2R", "HV3R"]),
        is_hv2f=any(p in filename for p in ["HV2F"]),
        is_hv2f_hv3f=any(p in filename for p in ["HV2F", "HV3F"]),
        is_hv2r_hv3r=any(p in filename for p in ["HV2R", "HV3R"]),
        is_hv1f=any(p in filename for p in ["HV1F"]),
        is_hv1r=any(p in filename for p in ["HV1R"]),
    )


# ---------------------------------------------------------------------------
# Heteroplasmy detection
# ---------------------------------------------------------------------------


def detect_heteroplasmy(
    peaks: list[float | None],
    position: float | str,
    heteroplasmy_threshold: float | None = None,
    sample: str | None = None,
) -> str:
    """Detect heteroplasmy from peak ratios and return IUPAC code for bases above threshold.

    Takes ``heteroplasmy_threshold`` as a parameter (defaults to settings).
    """
    threshold = (
        heteroplasmy_threshold if heteroplasmy_threshold is not None else get_settings().tracy.heteroplasmy_threshold
    )

    # Filter out None values
    peak_values = [p if p is not None else 0 for p in peaks]

    # Map peaks to bases
    base_peaks = {
        "A": peak_values[0],
        "C": peak_values[1],
        "G": peak_values[2],
        "T": peak_values[3],
    }

    # Find maximum peak
    max_peak = max(peak_values)

    if max_peak == 0:
        logger.warning(f"Sample {sample}: Position {position}: All peaks are zero")
        return "N"

    # Identify bases above heteroplasmy threshold
    threshold_value = max_peak * threshold
    detected_bases = [base for base, peak in base_peaks.items() if peak >= threshold_value]

    if len(detected_bases) == 0:
        logger.warning(f"Sample {sample}: Position {position}: No bases detected above threshold")
        return "N"

    # Single base — no heteroplasmy
    if len(detected_bases) == 1:
        return detected_bases[0]

    # Multiple bases — heteroplasmy
    iupac_code = nucleotide_to_iupac(detected_bases)
    logger.info(
        f"Sample {sample}: Heteroplasmy detected at position {position}: "
        f"{detected_bases} -> {iupac_code} "
        f"(peaks: {', '.join([f'{b}={base_peaks[b]:.0f}' for b in detected_bases])})",
    )
    return iupac_code


# ---------------------------------------------------------------------------
# Peak quality validation
# ---------------------------------------------------------------------------


def validate_peak_quality(  # noqa: PLR0911
    variant: dict[str, Any],
    min_peak_value: int | None = None,
    pratio: float | None = None,
    sample: str | None = None,
) -> bool:
    """Validate peak quality for a variant. Heteroplasmic variants auto-pass validation.

    Takes ``min_peak_value`` and ``pratio`` as parameters (default to settings).
    """
    tracy = get_settings().tracy
    min_peak_value = min_peak_value if min_peak_value is not None else tracy.min_peak_value
    pratio = pratio if pratio is not None else tracy.pratio

    # Special case: pos 73 ref='A' seq='G' always passes
    if pos_base(variant.get("pos", 0)) == POS_73 and variant.get("ref") == "A" and variant.get("seq") == "G":
        return True

    # Skip peak validation for deletions
    if variant["seq"] == "-":
        return True

    # Get peak values from the peaks array
    peak_values = variant["peaks"]

    # Filter out None values
    peak_values = [p for p in peak_values if p is not None]

    # Check if we have any peak values
    if not peak_values:
        return False

    # Check minimum peak value threshold
    max_peak = max(peak_values)

    if max_peak < min_peak_value:
        logger.warning(
            f"Sample {sample}: Remove variant: {variant['pos']} "
            f"{variant['ref']}>{variant['seq']} (peak quality: max_peak={max_peak} < threshold={min_peak_value})",
        )
        return False

    # Check if this is a heteroplasmic variant (IUPAC code)
    iupac_codes = "RYSWKMBDHVN"
    is_heteroplasmic = variant["seq"] in iupac_codes

    if is_heteroplasmic:
        # Heteroplasmic variants automatically pass validation
        return True

    # For homoplasmic variants, validate peak ratios
    total_peaks = sum(peak_values)
    if total_peaks == 0:
        return False

    highest_peak_ratio = max_peak / total_peaks
    if highest_peak_ratio <= pratio:
        logger.warning(
            f"Sample {sample}: Remove variant: {variant['pos']} {variant['ref']}>{variant['seq']} "
            f"(peak quality: peak_ratio={highest_peak_ratio:.2f} <= threshold={pratio})"
        )
        return False

    return True


# ---------------------------------------------------------------------------
# Variant condition detection — dispatch table
# ---------------------------------------------------------------------------

# Type alias for condition checker functions
_ConditionChecker = list[tuple[bool, str]]


def _check_deletion_conditions(pos_int: int, seq: str, ref: str) -> _ConditionChecker:
    """Check deletion-specific conditions, returning (flag, reason) pairs."""
    results: _ConditionChecker = []
    if pos_int == POS_248 and seq == "-":
        results.append((True, "has_248_deletion"))
    if pos_int == POS_250 and seq == "-":
        results.append((True, "has_250_deletion"))
    if pos_int == POS_456 and seq == "-":
        results.append((True, "has_456_deletion"))
    if pos_int == POS_513 and seq == "-":
        results.append((True, "has_513_deletion"))
        if ref == "G":
            results.append((True, "has_513_G_deletion"))
    if pos_int == POS_514 and seq == "-":
        results.append((True, "has_514_deletion"))
        if ref == "C":
            results.append((True, "has_514_C_deletion"))
    if pos_int == POS_515 and seq == "-":
        results.append((True, "has_515_deletion"))
    if pos_int == POS_16189 and seq == "-":
        results.append((True, "has_16189_deletion"))
    return results


def _check_snp_conditions(pos_int: int, seq: str, ref: str) -> _ConditionChecker:
    """Check SNP-specific conditions, returning (flag, reason) pairs."""
    results: _ConditionChecker = []
    if pos_int == POS_309 and seq == "T" and ref == "C":
        results.append((True, "has_309_C_T"))
    if pos_int == POS_310 and seq == "C" and ref == "T":
        results.append((True, "has_310_T_C"))
    if pos_int == POS_524 and seq == "T":
        results.append((True, "has_524_T"))
    if pos_int == POS_523 and seq == "C" and ref == "A":
        results.append((True, "has_523_A_C"))
    if pos_int == POS_16189 and seq == "C":
        results.append((True, "has_16189_T_C"))
    return results


def _check_insertion_conditions(pos: float | str, pos_int: int, seq: str, ref: str) -> _ConditionChecker:
    """Check insertion-specific conditions, returning (flag, reason) pairs."""
    results: _ConditionChecker = []
    if pos == "513.1" and seq == "C" and ref == "-":
        results.append((True, "has_513_1_insertion"))
    if pos == "513.2" and seq == "A" and ref == "-":
        results.append((True, "has_513_2_insertion"))
    if pos == "513.3" and seq == "C" and ref == "-":
        results.append((True, "has_513_3_insertion"))
    if pos == "513.4" and seq == "A" and ref == "-":
        results.append((True, "has_513_4_insertion"))
    # Count HV2_polyC insertions
    if POS_309 <= pos_int <= POS_316 and ref == "-":
        results.append((True, "HV2_polyC"))
    return results


def detect_variant_conditions(variants: list[dict[str, Any]]) -> dict[str, Any]:
    """Detect all variant conditions in a single pass.

    Uses a dispatch table pattern (``_check_deletion_conditions``,
    ``_check_snp_conditions``, ``_check_insertion_conditions``).
    """
    conditions: dict[str, Any] = {
        "has_248_deletion": False,
        "has_250_deletion": False,
        "has_309_C_T": False,
        "has_310_T_C": False,
        "has_456_deletion": False,
        "has_513_deletion": False,
        "has_513_1_insertion": False,
        "has_513_2_insertion": False,
        "has_513_3_insertion": False,
        "has_513_4_insertion": False,
        "has_514_deletion": False,
        "has_515_deletion": False,
        "has_524_T": False,
        "has_523_A_C": False,
        "has_16189_T_C": False,
        "has_16189_deletion": False,
        "has_513_G_deletion": False,
        "has_514_C_deletion": False,
        "HV2_polyC": 0,
    }

    for variant in variants:
        pos, seq, ref = variant["pos"], variant["seq"], variant["ref"]
        pos_int = pos_base(pos)

        conditions.update({key: flag for flag, key in _check_deletion_conditions(pos_int, seq, ref)})
        conditions.update({key: flag for flag, key in _check_snp_conditions(pos_int, seq, ref)})
        for flag, key in _check_insertion_conditions(pos, pos_int, seq, ref):
            if key == "HV2_polyC":
                conditions[key] += 1
            else:
                conditions[key] = flag

    return conditions
