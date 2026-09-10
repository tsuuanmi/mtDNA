"""Tracy variant transformations — position-specific, polyC, and strand transforms.

Pure functions that modify variant dicts in-place. These are internal
helpers for the ETL ``process()`` function in ``etl.py``.
"""

from typing import Any, NamedTuple

from Bio.Seq import Seq
from loguru import logger

from src.core.polyc import (
    directional_hv1_polyc_suppression_reason,
    directional_hv2_polyc_exclusion_reason,
)
from src.core.variants import normalize_position, pos_base
from src.tools.tracy.utils import (
    NUM_PEAK_CHANNELS,
    POS_248,
    POS_250,
    POS_309,
    POS_310,
    POS_316,
    POS_456,
    POS_459,
    POS_513,
    POS_514,
    POS_515,
    POS_523,
    POS_524,
    POS_525,
    POS_16189,
    POS_16193,
    PrimerTypeInfo,
    alignment_reference_position,
    detect_heteroplasmy,
    extract_peak_data,
)

# ---------------------------------------------------------------------------
# Parameter types for transforms
# ---------------------------------------------------------------------------


class _CreateParams(NamedTuple):
    """Parameters for creating a variant dict from alignment data."""

    ref_seq: str
    data: dict[str, Any]
    filename: str | None
    sample: str | None
    heteroplasmy_threshold: float


class _PositionParams(NamedTuple):
    """Parameters for position transformation based on strand direction."""

    primers: PrimerTypeInfo
    ref1pos: float
    ref_align: str
    filename: str
    sample: str | None
    heteroplasmy_threshold: float


# ---------------------------------------------------------------------------
# Variant creation helper
# ---------------------------------------------------------------------------


def _find_alignment_index(data: dict[str, Any], pos_float: float) -> int | None:
    """Find the alignment index closest to a target genomic position.

    ``ref1pos`` may be an int (starting position) or a normalized per-column
    position list.
    """
    if isinstance(data["ref1pos"], int):
        ref_start = data["ref1pos"]
        ref_align = data["ref1align"]
        best_idx = 0
        best_dist = float("inf")
        current_pos = ref_start
        for idx, ch in enumerate(ref_align):
            target = current_pos if ch != "-" else current_pos - 1
            dist = abs(target - pos_float)
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
            if ch != "-":
                current_pos += 1
        return None if best_dist == float("inf") else best_idx
    try:
        return min(range(len(data["ref1pos"])), key=lambda i: abs(data["ref1pos"][i] - pos_float))
    except ValueError:
        return None


# Fixed-shape created variants (deletion/insertion): ref from the reference
# sequence, fixed quality=50 and placeholder peaks=[1000,10,10,10], no
# heteroplasmy override.
_CREATED_FIXED_PEAKS: list[float] = [1000, 10, 10, 10]
_CREATED_FIXED_QUALITY = 50


def _create_variant(pos: int | str, variant_type: str, params: _CreateParams) -> dict[str, Any] | None:
    """Create a new variant dict at a given position.

    - deletion: ref from ``ref_seq[pos-1]``, seq "-", quality 50,
      peaks ``[1000,10,10,10]`` (no heteroplasmy override, no alignment lookup).
    - insertion (309-316 polyC): ref "-", seq "C", quality 50,
      peaks ``[1000,10,10,10]``.
    - snp: ref from ``ref_seq[pos-1]``, peaks/quality from the alignment column
      closest to ``pos`` (via ``_find_alignment_index``), seq from heteroplasmy
      detection.

    ``pos`` may be an int (base position) or a decimal string (e.g. "315.1")
    for insertions; reference lookup uses ``pos_base``.
    """
    pos_int = pos_base(pos)
    try:
        ref_base = params.ref_seq[pos_int - 1]
    except IndexError:
        logger.warning("Cannot create variant at position {}: out of reference range", pos)
        return None

    if variant_type == "deletion":
        return {
            "pos": pos,
            "ref": ref_base,
            "seq": "-",
            "peaks": list(_CREATED_FIXED_PEAKS),
            "index": pos_int,
            "peak_index": pos_int,
            "quality": _CREATED_FIXED_QUALITY,
        }

    # HV2 polyC insertion (309.x / 315.x)
    if variant_type == "insertion" and POS_309 <= pos_int <= POS_316:
        return {
            "pos": pos,
            "ref": "-",
            "seq": "C",
            "peaks": list(_CREATED_FIXED_PEAKS),
            "index": pos_int,
            "peak_index": pos_int,
            "quality": _CREATED_FIXED_QUALITY,
        }

    # SNP: look up the alignment column closest to pos for peak data
    data = params.data
    best_idx = _find_alignment_index(data, float(pos_int))
    if best_idx is None:
        return None
    peak_a, peak_c, peak_g, peak_t, quality = extract_peak_data(data, best_idx)
    peaks: list[float | None] = [peak_a, peak_c, peak_g, peak_t]
    if all(pk is None for pk in peaks):
        logger.warning("Cannot create SNP variant at position {}: no valid peak data", pos)
        return None
    detected_base = detect_heteroplasmy(peaks, pos, params.heteroplasmy_threshold, params.sample)
    return {
        "pos": pos,
        "ref": ref_base,
        "seq": detected_base,
        "peaks": peaks,
        "index": best_idx,
        "peak_index": best_idx,
        "quality": quality,
    }


# ---------------------------------------------------------------------------
# Position-specific transforms (dispatch table)
# ---------------------------------------------------------------------------
# The dispatch table maps position integers to transform functions.


def _transform_pos_248(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Convert 248/250 deletions to position 249 (ref='A')."""
    pos_int = pos_base(variant["pos"])
    if pos_int in (POS_248, POS_250) and (conditions.get("has_248_deletion") or conditions.get("has_250_deletion")):
        variant["pos"] = normalize_position(249)
        variant["ref"] = "A"


def _transform_pos_309(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Handle 309 C>T with 310 T>C — mark for removal.

    The 309 deletion variant is created separately in ``_apply_polyc_transforms``.
    """
    if conditions["has_309_C_T"] and conditions["has_310_T_C"]:
        variant["remove"] = True
        variant["reason"] = "309C>T and 310T>C conversion to 309 deletion + 315.1 insertion"


def _transform_pos_456(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Convert 456 deletion to 459."""
    if conditions["has_456_deletion"]:
        variant["pos"] = normalize_position(POS_459)


def _transform_pos_513(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Handle 513 deletions.

    - 513G- + 514C- → convert 513 deletion to G>A SNP (seq="A")
    - 513 deletion with 514 deletion → convert to 523 ref=A

    The 513.x insertion remapping (513.x -> 524.x) is handled separately by
    ``_apply_insertion_remap()`` on the decimal-formatted positions.
    """
    # Special case: 513G- + 514C- → convert deletion to G>A SNP
    if conditions.get("has_513_G_deletion") and conditions.get("has_514_C_deletion"):
        variant["seq"] = "A"
        variant["ref"] = "G"

    # Convert 513/514 deletions to 523/524 (exclude 513G- + 514C- case)
    if (
        conditions["has_513_deletion"]
        and conditions["has_514_deletion"]
        and not (conditions.get("has_513_G_deletion") and conditions.get("has_514_C_deletion"))
    ):
        variant["pos"] = normalize_position(POS_523)
        variant["ref"] = "A"


def _transform_pos_514(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Handle 514 deletions.

    - 514C- with 513G- → remove 514, create 523/524 deletions
    - 513/514 deletions (exclude special case) → convert 514 to 524 ref=C
    - 514/515 deletions → convert 514 to 523 ref=A
    """
    # Special case: 513G- + 514C- → remove 514
    if conditions.get("has_513_G_deletion") and conditions.get("has_514_C_deletion"):
        variant["remove"] = True
        variant["reason"] = "514C- deletion converted to 523A- and 524C- deletions"

    # Convert 513/514 deletions to 523/524 (exclude special case)
    if (
        conditions["has_513_deletion"]
        and conditions["has_514_deletion"]
        and not (conditions.get("has_513_G_deletion") and conditions.get("has_514_C_deletion"))
    ):
        variant["pos"] = normalize_position(POS_524)
        variant["ref"] = "C"

    # Convert 514/515 deletions to 523/524
    if conditions["has_514_deletion"] and conditions["has_515_deletion"]:
        variant["pos"] = normalize_position(POS_523)
        variant["ref"] = "A"


def _transform_pos_515(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Handle 515 deletion (with 514 deletion → 524 ref=C)."""
    if conditions["has_514_deletion"] and conditions["has_515_deletion"]:
        variant["pos"] = normalize_position(POS_524)
        variant["ref"] = "C"


def _transform_pos_523(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Convert 523A>C to 523A-."""
    if conditions.get("has_523_A_C"):
        variant["seq"] = "-"


def _transform_pos_16189(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Handle 16189 deletion (seq → "C")."""
    if conditions["has_16189_deletion"]:
        variant["seq"] = "C"


def _apply_position_specific_transforms(
    variant: dict[str, Any],
    conditions: dict[str, Any],
) -> None:
    """Apply position-specific transformations to a variant in place.

    Args:
        variant: Variant dict to transform.
        conditions: Detected variant conditions.
    """
    pos_int = pos_base(variant["pos"])
    for transform in _POSITION_TRANSFORMS.get(pos_int, []):
        transform(variant, conditions)


# Dispatch table for position-specific transforms (non-polyC)
_POSITION_TRANSFORMS: dict[int, list] = {
    POS_248: [_transform_pos_248],
    POS_250: [_transform_pos_248],  # 250 also maps to _transform_pos_248 (same logic: convert to 249)
    POS_309: [_transform_pos_309],
    POS_456: [_transform_pos_456],
    POS_513: [_transform_pos_513],
    POS_514: [_transform_pos_514],
    POS_515: [_transform_pos_515],
    POS_523: [_transform_pos_523],
    POS_16189: [_transform_pos_16189],
}


# ---------------------------------------------------------------------------
# PolyC transformations
# ---------------------------------------------------------------------------


def _apply_polyc_insertion(variant: dict[str, Any]) -> None:
    """Replace raw HV2 polyC insertions with canonical generated insertions."""
    pos_int = pos_base(variant["pos"])
    if POS_309 <= pos_int <= POS_316 and variant["ref"] == "-" and not variant.get("keep", False):
        variant["remove"] = True
        variant["reason"] = "HV2 polyC insertion normalized to canonical insertion"


def _apply_polyc_primer_filters(
    variant: dict[str, Any],
    primers: PrimerTypeInfo,
    conditions: dict[str, Any],
) -> None:
    """Apply shared HV2 and HV1 directional polyC filters after calling."""
    pos_int = pos_base(variant["pos"])

    # The shared policy owns directional exclusion coverage and reasons.
    # Tracy only maps its primer metadata to canonical trace directions.
    forward_reason = directional_hv2_polyc_exclusion_reason(
        variant["pos"],
        "forward" if primers.is_hv2f_hv3f else None,
    )
    if forward_reason is not None:
        variant["remove"] = True
        variant["reason"] = forward_reason

    reverse_reason = directional_hv2_polyc_exclusion_reason(
        variant["pos"],
        "reverse" if primers.is_hv2r_hv3r else None,
    )
    if reverse_reason is not None:
        variant["remove"] = True
        variant["reason"] = reverse_reason

    # Remove variants at position > 525 for HV2F/HV3F primers
    if pos_int > POS_525 and primers.is_hv2f_hv3f:
        variant["remove"] = True
        variant["reason"] = "Variant in position >= 525 (HV2F/HV3F primers)"

    # The shared policy owns the HV1 trigger, directional boundary, and reason.
    for direction in (
        "forward" if primers.is_hv1f else None,
        "reverse" if primers.is_hv1r else None,
    ):
        hv1_reason = directional_hv1_polyc_suppression_reason(
            variant["pos"],
            direction,
            polyc_created=conditions.get("has_16189_T_C", False),
        )
        if hv1_reason is not None:
            variant["remove"] = True
            variant["reason"] = hv1_reason


def _apply_conversion_deletions(
    variant: dict[str, Any],
    conditions: dict[str, Any],
    params: _CreateParams,
    additional_variants: list[dict[str, Any]],
) -> None:
    """Create deletion variants for the 309/514/524 conversion cases.

    These need alignment/ref-seq lookup (so they live here where ``params`` is
    available, not in the no-params dispatch transforms). A variant is at one
    position, so the cases are mutually exclusive.
    """
    pos_int = pos_base(variant["pos"])
    # 309C>T + 310T>C → kept 309 deletion (309 SNP removal is in _transform_pos_309)
    if pos_int == POS_309 and conditions["has_309_C_T"] and conditions["has_310_T_C"]:
        deletion_309 = _create_variant(POS_309, "deletion", params)
        if deletion_309 is not None:
            deletion_309["keep"] = True
            additional_variants.append(deletion_309)
        return
    # 513G- + 514C- → 523A- and 524C- (514 removal is in _transform_pos_514)
    if pos_int == POS_514 and conditions.get("has_513_G_deletion") and conditions.get("has_514_C_deletion"):
        for position in (POS_523, POS_524):
            created = _create_variant(position, "deletion", params)
            if created is not None:
                additional_variants.append(created)
        return
    # 523- + 524T → 522T + 523- + 524-
    if pos_int == POS_524 and conditions.get("has_524_T"):
        variant["pos"] = normalize_position(522)
        deletion_524 = _create_variant(POS_524, "deletion", params)
        if deletion_524 is not None:
            additional_variants.append(deletion_524)


def _apply_polyc_transforms(
    variant: dict[str, Any],
    conditions: dict[str, Any],
    primers: PrimerTypeInfo,
    params: _CreateParams,
    additional_variants: list[dict[str, Any]],
) -> None:
    """Apply polyC-specific transformations to a variant."""
    # Normalize raw polyC insertions; canonical insertions are created later.
    _apply_polyc_insertion(variant)

    # Primer-specific polyC filtering
    _apply_polyc_primer_filters(variant, primers, conditions)

    # Special case: 16189 deletion with 16193 deletion
    if pos_base(variant["pos"]) == POS_16189 and conditions["has_16189_deletion"]:
        variant["seq"] = "C"
        deletion_16193 = _create_variant(POS_16193, "deletion", params)
        if deletion_16193 is not None:
            additional_variants.append(deletion_16193)

    # Conversion deletion creations (309/514/524); a variant sits at one
    # position so the branches are mutually exclusive.
    _apply_conversion_deletions(variant, conditions, params, additional_variants)

    # 310 T>C with 309 C>T → remove (309 handling is in _transform_pos_309)
    if pos_base(variant["pos"]) == POS_310 and conditions["has_309_C_T"] and conditions["has_310_T_C"]:
        variant["remove"] = True
        variant["reason"] = "309C>T and 310T>C conversion to 309 deletion + 315.1 insertion"


def _apply_insertion_remap(variant: dict[str, Any], conditions: dict[str, Any]) -> None:
    """Remap AC-repeat insertions 513.x → 524.x when paired conditions hold.

    Operates on the decimal-formatted insertion positions (e.g. "513.1") produced
    by ``update_variant_positions()``, so it must run after the position update.
    """
    insertion_map: dict[str, tuple[str, str, bool]] = {
        "513.1": (
            "524.1",
            "A",
            conditions["has_513_1_insertion"] and conditions["has_513_2_insertion"],
        ),
        "513.2": (
            "524.2",
            "C",
            conditions["has_513_1_insertion"] and conditions["has_513_2_insertion"],
        ),
        "513.3": (
            "524.3",
            "A",
            conditions["has_513_3_insertion"] and conditions["has_513_4_insertion"],
        ),
        "513.4": (
            "524.4",
            "C",
            conditions["has_513_3_insertion"] and conditions["has_513_4_insertion"],
        ),
    }
    pos = variant["pos"]
    if pos in insertion_map:
        new_pos, new_seq, should_convert = insertion_map[pos]
        if should_convert:
            variant["pos"] = normalize_position(new_pos)
            variant["seq"] = new_seq


_POLYC_PATTERNS: dict[int, list[str]] = {
    1: ["315.1"],
    2: ["309.1", "315.1"],
    3: ["309.1", "309.2", "315.1"],
}


def _create_polyc_insertions(
    variants: list[dict[str, Any]],
    conditions: dict[str, Any],
) -> None:
    """Create canonical HV2 polyC insertions from the HV2_polyC count.

    Primer-agnostic: created for every file whose count > 0 (HV2F and HV3R both
    produce 315.1). The count is ``conditions["HV2_polyC"] + (1 if 309C>T & 310T>C)``.
    Each created insertion has ref="-", seq="C", quality=50, peaks=[1000,10,10,10],
    with the decimal-string position preserved. Raw polyC insertions are marked
    for canonical normalization before this runs, so this does not duplicate them.
    """
    hv2_polyc = conditions["HV2_polyC"]
    if conditions["has_309_C_T"] and conditions["has_310_T_C"]:
        hv2_polyc += 1
    variants.extend(
        {
            "pos": normalize_position(pos_str),
            "ref": "-",
            "seq": "C",
            "peaks": [1000, 10, 10, 10],
            "index": pos_base(pos_str),
            "peak_index": pos_base(pos_str),
            "quality": 50,
        }
        for pos_str in _POLYC_PATTERNS.get(hv2_polyc, [])
    )


def apply_all_transformations(
    variants: list[dict[str, Any]],
    conditions: dict[str, Any],
    primers: PrimerTypeInfo,
    params: _CreateParams,
) -> None:
    """Apply all variant transformations in-place.

    Order: position-specific dispatch → insertion remap (513.x → 524.x) →
    polyC transforms → additional variant creation. Insertion remap runs on
    decimal positions produced by ``update_variant_positions()``, which runs
    before this function in ``etl.process()``.
    """
    additional_variants: list[dict[str, Any]] = []

    for variant in variants:
        # Position-specific transforms (non-polyC)
        _apply_position_specific_transforms(variant, conditions)

        # AC-repeat insertion remap (513.x → 524.x) — runs on decimal positions
        _apply_insertion_remap(variant, conditions)

        # PolyC transforms
        _apply_polyc_transforms(variant, conditions, primers, params, additional_variants)

    # Append any additional variants created by transforms
    variants.extend(additional_variants)

    # Create canonical HV2 polyC insertions from the HV2_polyC count. Runs after
    # the transform loop (primer removal already handled the polyC variants).
    _create_polyc_insertions(variants, conditions)


# ---------------------------------------------------------------------------
# Position/strand transforms
# ---------------------------------------------------------------------------


def _compute_file_intervals(
    params: _PositionParams,
    start_index: int | None,
    end_index: int | None,
) -> dict[str, list[list[int]]]:
    """Compute file intervals based on strand direction."""
    if start_index is None or end_index is None:
        return {}

    start_pos = alignment_reference_position(
        start_index,
        ref_start=int(params.ref1pos),
        ref_align=params.ref_align,
        is_reverse=params.primers.is_reverse,
    )
    end_pos = alignment_reference_position(
        end_index,
        ref_start=int(params.ref1pos),
        ref_align=params.ref_align,
        is_reverse=params.primers.is_reverse,
    )
    return {params.filename: [[min(start_pos, end_pos), max(start_pos, end_pos)]]}


def _transform_reverse_strand(variant: dict[str, Any], params: _PositionParams) -> None:
    """Transform a variant for reverse strand: position, bases, and peaks."""
    variant["pos"] = normalize_position(
        alignment_reference_position(
            variant["index"],
            ref_start=int(params.ref1pos),
            ref_align=params.ref_align,
            is_reverse=True,
        ),
    )
    # Transform bases to complements for reverse strand
    if variant["ref"] != "-":
        variant["ref"] = str(Seq(variant["ref"]).reverse_complement())
    if variant["seq"] != "-":
        variant["seq"] = str(Seq(variant["seq"]).reverse_complement())
    # Swap peaks: A<->T, G<->C
    peaks = variant["peaks"]
    if len(peaks) >= NUM_PEAK_CHANNELS:
        if peaks[0] is not None and peaks[3] is not None:
            peaks[0], peaks[3] = peaks[3], peaks[0]  # A<->T
        if peaks[1] is not None and peaks[2] is not None:
            peaks[1], peaks[2] = peaks[2], peaks[1]  # G<->C
        variant["peaks"] = peaks


def _redetect_heteroplasmy(variant: dict[str, Any], params: _PositionParams) -> None:
    """Re-detect heteroplasmy after strand transformation (SNPs only)."""
    if variant["ref"] == "-" or variant["seq"] == "-":
        return
    peaks = variant["peaks"]
    if len(peaks) < NUM_PEAK_CHANNELS or not any(p is not None for p in peaks):
        return
    detected_base = detect_heteroplasmy(peaks, variant["pos"], params.heteroplasmy_threshold, params.sample)
    variant["seq"] = detected_base
    if variant["ref"] == variant["seq"]:
        variant["remove"] = True
        variant["reason"] = "Reference equals sequence after heteroplasmy detection"


def _transform_variant_position(
    variant: dict[str, Any],
    params: _PositionParams,
    insertion_count_at_pos: dict[int, int],
) -> None:
    """Transform a single variant's position and bases based on strand direction.

    Consecutive insertions at the same base position are numbered with a running
    counter (e.g., 315.1, 315.2, 315.3 and 513.1, 513.2, 513.3, 513.4), as
    required by the AC-repeat insertion remap (513.x -> 524.x).
    """
    if params.primers.is_forward or params.primers.is_hv2f:
        variant["pos"] = normalize_position(
            alignment_reference_position(
                variant["index"],
                ref_start=int(params.ref1pos),
                ref_align=params.ref_align,
                is_reverse=False,
            ),
        )
    elif params.primers.is_reverse:
        _transform_reverse_strand(variant, params)

    # Re-detect heteroplasmy after strand transformation (SNPs only)
    _redetect_heteroplasmy(variant, params)

    # Format insertion positions (e.g., 316 insert C -> 315.1); consecutive
    # insertions at the same base position get .1, .2, .3, ... suffixes.
    if variant["ref"] == "-":
        base_pos = pos_base(variant["pos"]) - 1
        if base_pos not in insertion_count_at_pos:
            insertion_count_at_pos[base_pos] = 0
        insertion_count_at_pos[base_pos] += 1
        variant["pos"] = f"{base_pos}.{insertion_count_at_pos[base_pos]}"


def update_variant_positions(
    variants: list[dict[str, Any]],
    position_params: _PositionParams,
    start_index: int | None,
    end_index: int | None,
) -> dict[str, list[list[int]]]:
    """Update variant positions based on strand direction and compute intervals.

    Consecutive insertions at the same base position are numbered via a shared
    ``insertion_count_at_pos`` counter, as needed by the AC-repeat insertion
    remap (513.x -> 524.x).
    """
    file_intervals = _compute_file_intervals(position_params, start_index, end_index)

    # Process variants by strand direction
    variants_to_process = variants
    if position_params.primers.is_reverse:
        variants_to_process = sorted(variants, key=lambda x: x["index"], reverse=True)

    insertion_count_at_pos: dict[int, int] = {}
    for variant in variants_to_process:
        _transform_variant_position(variant, position_params, insertion_count_at_pos)

    return file_intervals
