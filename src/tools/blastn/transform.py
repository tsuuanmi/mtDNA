"""BLASTn variant transforms — indel normalization, conversion rules, polyC metadata.

Extracted from ``etl.py`` to keep the ETL orchestration module focused on
variant calling, consensus, intervals, and ``process()`` orchestration. This
module mirrors the ``transforms.py`` convention used by ``src/tools/tracy/``.

Contents:
- Indel right-normalization: ``_normalize_single_indel``, ``_normalize_indel_block``,
  ``_normalize_indel_blocks``, ``_right_norm_indel``
- YAML conversion rules: ``_load_conversion_rules``, ``_check_must_have``,
  ``_check_must_not_have``, ``_apply_rule_changes``, ``_parse_variant_rule``,
  ``_apply_conversion_rules``
- polyC post-hoc removal: ``_polyc_stretch_detected``, ``_handle_polyc_removal``,
  ``_apply_in_polyc_rule``
- region validation: ``_regions_from_files``, ``_apply_region_validation_rules``
- polyC warning flags: ``_polyc_warning_flags``
"""

from pathlib import Path
from typing import Any

import yaml
from loguru import logger

from src.config import get_settings
from src.core.models import Position
from src.core.variants import get_hv_region_for_position, normalize_position, pos_base, pos_sort_key

__all__ = [
    "_apply_conversion_rules",
    "_apply_region_validation_rules",
    "_handle_polyc_removal",
    "_polyc_warning_flags",
    "_right_norm_indel",
]

# ---------------------------------------------------------------------------
# Conversion-rule loader (YAML)
# ---------------------------------------------------------------------------


def _load_conversion_rules(convert_rules_dir: str | None = None) -> list[dict[str, Any]]:
    """Load YAML conversion rules from the convert rules directory.

    Conversion rules allow adding or removing specific variants based on
    conditions (MustHave, MustNotHave).

    Args:
        convert_rules_dir: Path to conversion rules directory.
            If None, uses settings.directories.convert_rules.

    Returns:
        List of parsed YAML conversion rule dicts.
    """
    rules: list[dict[str, Any]] = []

    convert_rules_dir = str(get_settings().directories.convert_rules)

    rules_dir = Path(convert_rules_dir)
    if not rules_dir.is_dir():
        logger.warning("Conversion rules directory not found: {}", convert_rules_dir)
        return rules

    for rule_file in sorted(rules_dir.glob("*.yaml")) + sorted(rules_dir.glob("*.yml")):
        try:
            with rule_file.open() as f:
                rule = yaml.safe_load(f)
            if rule:
                rules.append(rule)
        except (OSError, yaml.YAMLError) as e:
            logger.error("Error loading conversion rule {}: {}", rule_file, e)

    return rules


# ---------------------------------------------------------------------------
# PolyC region helpers
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Right-normalization of indels
# ---------------------------------------------------------------------------


def _normalize_single_indel(variant: dict[str, Any], rcrs_seq: str) -> dict[str, Any]:
    """Normalize a single indel variant."""
    if variant["ref"] == "-":  # Insertion
        new_pos = pos_base(variant["pos"]) + 1
        ref_seq = ""
        seq_seq = variant["seq"]
    else:  # Deletion
        new_pos = pos_base(variant["pos"]) + 1
        ref_seq = variant["ref"]
        seq_seq = ""

    can_right_norm = True
    while can_right_norm and new_pos < len(rcrs_seq):
        base_at_new_pos = rcrs_seq[new_pos - 1 : new_pos]  # 0-based indexing
        ref_seq = ref_seq + base_at_new_pos
        seq_seq = seq_seq + base_at_new_pos

        if len(ref_seq) > 0 and len(seq_seq) > 0 and ref_seq[0] == seq_seq[0]:
            ref_seq = ref_seq[1:]
            seq_seq = seq_seq[1:]
            new_pos += 1
        else:
            can_right_norm = False

    new_pos = normalize_position(f"{new_pos - 1}.1") if variant["ref"] == "-" else new_pos - 1

    return {"pos": new_pos, "ref": variant["ref"], "seq": variant["seq"]}


def _normalize_indel_block(block: list[dict[str, Any]], rcrs_seq: str, indel_type: str) -> list[dict[str, Any]]:
    """Normalize consecutive indels as a single unit."""
    start_pos = min(pos_base(var["pos"]) for var in block)

    if indel_type == "deletion":
        deleted_bases = ""
        for var in sorted(block, key=lambda x: pos_sort_key(x["pos"])):
            deleted_bases += var["ref"]

        normalized_deletions: list[dict[str, Any]] = []
        for i, base in enumerate(deleted_bases):
            normalized_deletions.append({"pos": normalize_position(start_pos + i), "ref": base, "seq": "-"})
        return normalized_deletions

    # insertion block — normalize insertions individually
    normalized_insertions: list[dict[str, Any]] = []
    for i, variant in enumerate(sorted(block, key=lambda x: pos_sort_key(x["pos"]))):
        normalized_var = _normalize_single_indel(variant, rcrs_seq)
        if i > 0:
            base = pos_base(normalized_var["pos"])
            normalized_var["pos"] = normalize_position(f"{base}.{i + 1}")
        normalized_insertions.append(normalized_var)
    return normalized_insertions


def _normalize_indel_blocks(variants: list[dict[str, Any]], rcrs_seq: str, indel_type: str) -> list[dict[str, Any]]:
    """Group consecutive indels into blocks and normalize each block."""
    if not variants:
        return []

    blocks: list[list[dict[str, Any]]] = []
    current_block = [variants[0]]

    for i in range(1, len(variants)):
        prev_pos = pos_base(variants[i - 1]["pos"])
        curr_pos = pos_base(variants[i]["pos"])

        if curr_pos == prev_pos + 1:
            current_block.append(variants[i])
        else:
            blocks.append(current_block)
            current_block = [variants[i]]

    blocks.append(current_block)

    normalized_variants: list[dict[str, Any]] = []
    for block in blocks:
        if len(block) == 1:
            normalized_variants.append(_normalize_single_indel(block[0], rcrs_seq))
        else:
            normalized_variants.extend(_normalize_indel_block(block, rcrs_seq, indel_type))

    return normalized_variants


def _right_norm_indel(variants: list[dict[str, Any]], rcrs_seq: str) -> list[dict[str, Any]]:
    """Right-normalize indels by grouping consecutive indels.

    Separates insertions and deletions, sorts each by position, groups
    consecutive indels into blocks, and normalizes each block.
    """
    if not variants:
        return []

    insertions = [v for v in variants if v["ref"] == "-"]
    deletions = [v for v in variants if v["seq"] == "-"]

    insertions = sorted(insertions, key=lambda x: pos_sort_key(x["pos"]))
    deletions = sorted(deletions, key=lambda x: pos_sort_key(x["pos"]))

    norm_variants: list[dict[str, Any]] = []
    norm_variants.extend(_normalize_indel_blocks(insertions, rcrs_seq, "insertion"))
    norm_variants.extend(_normalize_indel_blocks(deletions, rcrs_seq, "deletion"))

    return norm_variants


# ---------------------------------------------------------------------------
# Collect variants by position
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Conversion rules
# ---------------------------------------------------------------------------


def _check_must_have(must_have: list[str], variants_dict: dict[Any, dict[str, Any]]) -> bool:
    """Check if all MustHave variants are present.

    Args:
        must_have: List of variant rule strings.
        variants_dict: Position-keyed variant data.

    Returns:
        True if all MustHave conditions are satisfied.
    """
    for variant_str in must_have:
        pos, ref, seq = _parse_variant_rule(variant_str)
        pos_norm = normalize_position(pos)
        if pos_norm not in variants_dict:
            return False
        v = variants_dict[pos_norm]
        if v.get("ref") != ref and ref != "*":
            return False
        if v.get("seq") != seq and seq != "*":
            return False
    return True


def _check_must_not_have(must_not_have: list[str], variants_dict: dict[Any, dict[str, Any]]) -> bool:
    """Check if any MustNotHave variant is present.

    Args:
        must_not_have: List of variant rule strings.
        variants_dict: Position-keyed variant data.

    Returns:
        True if any MustNotHave variant is present (meaning rule should NOT apply).
    """
    for variant_str in must_not_have:
        pos, ref, seq = _parse_variant_rule(variant_str)
        pos_norm = normalize_position(pos)
        if pos_norm in variants_dict:
            v = variants_dict[pos_norm]
            if (v.get("ref") == ref or ref == "*") and (v.get("seq") == seq or seq == "*"):
                return True
    return False


def _apply_rule_changes(changes: dict[str, Any], variants_dict: dict[Any, dict[str, Any]]) -> None:
    """Apply Add/Remove changes from a conversion rule to the variants dict.

    Mutates variants_dict in place.

    Args:
        changes: Dict with 'Remove' and 'Add' lists.
        variants_dict: Position-keyed variant data.
    """
    for variant_str in changes.get("Remove", []):
        pos, _, _ = _parse_variant_rule(variant_str)
        pos_norm = normalize_position(pos)
        variants_dict.pop(pos_norm, None)
        logger.info("Conversion rule: Remove variant at {}", pos_norm)

    for variant_str in changes.get("Add", []):
        pos, ref, seq = _parse_variant_rule(variant_str)
        pos_norm = normalize_position(pos)
        variants_dict[pos_norm] = {"pos": pos_norm, "ref": ref, "seq": seq, "file": ["conversion_rule"]}
        logger.info("Conversion rule: Add variant at {} {} {}", pos_norm, ref, seq)


def _apply_conversion_rules(
    variants: list[dict[str, Any]],
    convert_rules_dir: str | None = None,
) -> list[dict[str, Any]]:
    """Apply YAML conversion rules to variants.

    Conversion rules can add or remove specific variants based on
    MustHave/MustNotHave conditions.

    Args:
        variants: List of variant dicts.
        convert_rules_dir: Directory containing YAML conversion rules.

    Returns:
        List of variants after applying conversion rules.
    """
    rules = _load_conversion_rules(convert_rules_dir)
    if not rules:
        return variants

    # Convert to position-keyed dict for rule matching
    variants_dict: dict[Any, dict[str, Any]] = {}
    for v in variants:
        key = v.get("pos")
        if key is not None:
            variants_dict[key] = v

    for rule in rules:
        # Check MustHave conditions
        must_have = rule.get("MustHave", [])
        if must_have and not _check_must_have(must_have, variants_dict):
            continue

        # Check MustNotHave conditions
        must_not_have = rule.get("MustNotHave", [])
        if must_not_have and _check_must_not_have(must_not_have, variants_dict):
            continue

        # Apply changes
        changes = rule.get("Change", {})
        _apply_rule_changes(changes, variants_dict)

    return list(variants_dict.values())


def _parse_variant_rule(variant_str: str) -> tuple[Position, str, str]:
    """Parse a variant rule string like '309.1TC' or '16193CT'.

    Args:
        variant_str: Variant string in format 'positionREF+ALT'.

    Returns:
        Tuple of (position, ref, seq).
    """
    pos_str = variant_str[:-2]
    ref = variant_str[-2]
    seq = variant_str[-1]
    return normalize_position(pos_str), ref, seq


# ---------------------------------------------------------------------------
# polyC region handling
# ---------------------------------------------------------------------------

# Max per-trace Phred score required to keep a low-confidence variant that
# falls in / after a polyC region or in a non-native HV region.
_QUALITY_KEEP_THRESHOLD: int = 40

# HV1 position whose T>C transition creates a polyC stretch.
_HV1_POLYC_CREATING_POS: int = 16189

# HV2 is considered to have a polyC stretch when the reference block contains
# at least this many C bases.
_HV2_C_COUNT_THRESHOLD: int = 7


def _polyc_stretch_detected(
    region_name: str,
    variants_by_pos: dict[Any, dict[str, Any]],
    ref_seq: str,
) -> bool:
    """Check whether a polyC stretch is detected for ``region_name``.

    Looks at the called variants (16189T>C for HV1; 309.1/309.2/315.1 or a
    C-rich reference block for HV2) to decide whether a region is "affected"
    and therefore subject to special polyC handling.

    Args:
        region_name: Region key in ``regions.POLYC_REGIONS``.
        variants_by_pos: Position-keyed variant dicts (consensus output).
        ref_seq: Reference sequence string (rCRS).

    Returns:
        True if a polyC stretch is detected for the region.
    """
    polyc_regions = get_settings().regions.POLYC_REGIONS

    if region_name == "HV1":
        v16189 = variants_by_pos.get(_HV1_POLYC_CREATING_POS)
        return v16189 is not None and v16189.get("seq") == "C"

    if region_name == "HV2":
        if any(p in variants_by_pos for p in ("309.1", "309.2", "315.1")):
            return True
        region_range = polyc_regions.get("HV2")
        if not region_range:
            return False
        seq_section = ref_seq[region_range[0] - 1 : region_range[1]]
        return seq_section.count("C") >= _HV2_C_COUNT_THRESHOLD

    return False


def _max_quality(variant: dict[str, Any]) -> int:
    """Return the max per-trace quality score, or -1 when unavailable.

    Returns -1 for empty quality so the keep-threshold fails (treated as
    low-quality / removable).
    """
    quality = variant.get("quality") or []
    if isinstance(quality, list):
        return max((int(q) for q in quality if q is not None), default=-1)
    if quality is not None:
        try:
            return int(quality)
        except (TypeError, ValueError):
            return -1
    return -1


def _keep_by_quality(variant: dict[str, Any]) -> bool:
    """True when a variant's max per-trace quality meets the keep threshold."""
    return _max_quality(variant) >= _QUALITY_KEEP_THRESHOLD


def _region_for_pos(
    pos_int: int,
    regions: dict[str, list[int]],
) -> str | None:
    """Return the region name whose range contains ``pos_int``, or None."""
    for rname, rrange in regions.items():
        if rrange[0] <= pos_int <= rrange[1]:
            return rname
    return None


def _apply_in_polyc_rule(
    variant: dict[str, Any],
    region: str,
) -> list[dict[str, Any]]:
    """Apply special handling for a variant inside an affected polyC region.

    Returns the list of retained variants (0 or 1 element) for this variant.
    """
    pos = variant.get("pos", 0)
    pos_int = pos_base(pos)

    if region == "HV2":
        # Keep 309.1C/309.2C; retain all other HV2 polyC variants.
        if (
            isinstance(pos, str)
            and "." in pos
            and pos_base(pos) == 309  # noqa: PLR2004
            and pos in ("309.1", "309.2")
            and variant.get("seq") == "C"
        ):
            variant["polyc_region"] = True
        return [variant]

    if region == "HV1":
        # Keep 16189T>C; retain all other HV1 polyC variants.
        if pos_int == _HV1_POLYC_CREATING_POS and variant.get("ref") == "T" and variant.get("seq") == "C":
            variant["polyc_region"] = True
        return [variant]

    # HV3 and any other non-HV1/HV2 polyC region: quality-gated removal.
    if _keep_by_quality(variant):
        variant["polyc_region"] = True
        variant["kept_by_quality"] = True
        return [variant]
    logger.info(
        "Remove {} {} {} filter_polyc_{}: quality=[] or low",
        pos,
        variant.get("ref"),
        variant.get("seq"),
        region,
    )
    return []


def _handle_polyc_removal(
    variants: list[dict[str, Any]],
    ref_seq: str,
) -> list[dict[str, Any]]:
    """Handle variants in/after polyC regions, removing low-quality calls.

    This is the post-hoc REMOVAL pass (distinct from the consensus-level
    ``flagged_for_filtering`` filter already applied in
    ``_compute_consensus_variants``):

    - Variants in a ``SKIP_AFTER_POLYC`` range are removed unless their max
      per-trace quality is >= 40 (kept variants are marked ``kept_by_quality``).
    - Variants inside a ``POLYC_REGIONS`` region where a polyC stretch is
      detected are specially handled: 309.1C/309.2C (HV2) and 16189T>C (HV1)
      are always kept; HV3 / other non-HV1-HV2 polyC variants are removed
      unless max quality >= 40. Other in-region variants are retained.

    Variants kept by quality or by a special-case rule are tagged with
    ``polyc_region`` (and ``kept_by_quality`` where applicable) so the
    downstream warning-flags pass can annotate them.

    Args:
        variants: List of consensus variant dicts (pos/ref/seq/quality/file).
        ref_seq: Reference sequence string (rCRS).

    Returns:
        Filtered variant list with polyC metadata on retained variants.
    """
    polyc_regions = get_settings().regions.POLYC_REGIONS
    skip_after = get_settings().regions.SKIP_AFTER_POLYC

    variants_by_pos: dict[Any, dict[str, Any]] = {v.get("pos"): v for v in variants}
    affected_regions = {region for region in polyc_regions if _polyc_stretch_detected(region, variants_by_pos, ref_seq)}

    kept: list[dict[str, Any]] = []
    for v in variants:
        pos_int = pos_base(v.get("pos", 0))

        after_region = _region_for_pos(pos_int, skip_after)
        if after_region is not None:
            if _keep_by_quality(v):
                v["polyc_region"] = True
                v["kept_by_quality"] = True
                kept.append(v)
            else:
                logger.info(
                    "Remove {} {} {} after_polyc_region_{}: quality=[] or low",
                    v.get("pos"),
                    v.get("ref"),
                    v.get("seq"),
                    after_region,
                )
            continue

        in_region = _region_for_pos(pos_int, polyc_regions)
        if in_region is not None and in_region in affected_regions:
            kept.extend(_apply_in_polyc_rule(v, in_region))
            continue

        kept.append(v)

    return kept


# ---------------------------------------------------------------------------
# Region validation
# ---------------------------------------------------------------------------


def _regions_from_files(files: list[str] | None) -> list[str]:
    """Extract HV region names from a variant's source file names.

    A file name containing "HV1"/"HV2"/"HV3" contributes that region label.
    """
    if not files:
        return []
    regions: list[str] = []
    for f in files:
        if not isinstance(f, str):
            continue
        if "HV1" in f:
            regions.append("HV1")
        elif "HV2" in f:
            regions.append("HV2")
        elif "HV3" in f:
            regions.append("HV3")
    return list(set(regions))


def _apply_region_validation_rules(
    variants: list[dict[str, Any]],
    no_filter_positions: set[str] | None = None,
) -> list[dict[str, Any]]:
    """Reject variants seen only in non-native HV regions with low quality.

    For each variant: if it is detected only in regions other than its native
    region (derived from its position via ``get_hv_region_for_position``) and
    its max per-trace quality is < 40, it is removed. Positions listed in
    ``no_filter_positions`` are always kept.

    Args:
        variants: List of variant dicts with ``file`` (source names) and
            ``quality`` fields.
        no_filter_positions: Position strings that should never be removed.
            If None, loaded from rules/filter/no_filter.txt.

    Returns:
        Filtered variant list.
    """
    if no_filter_positions is None:
        from src.tools.blastn.etl import _load_no_filter_positions  # noqa: PLC0415

        no_filter_positions = _load_no_filter_positions()

    kept: list[dict[str, Any]] = []
    for v in variants:
        pos = v.get("pos", 0)
        pos_str = str(normalize_position(pos))

        if pos_str in no_filter_positions:
            kept.append(v)
            continue

        native_region = get_hv_region_for_position(pos)
        if not native_region:
            kept.append(v)
            continue

        detected_regions = _regions_from_files(v.get("file"))
        if not detected_regions or native_region in detected_regions:
            kept.append(v)
            continue

        # Detected only in non-native regions: quality-gated removal.
        if _keep_by_quality(v):
            logger.info(
                "Keep {} {} {} region_validation: kept due to high quality",
                pos,
                v.get("ref"),
                v.get("seq"),
            )
            kept.append(v)
        else:
            logger.warning(
                "Remove {} {} {} region_validation: only in {} (native {})",
                pos,
                v.get("ref"),
                v.get("seq"),
                ",".join(detected_regions),
                native_region,
            )
    return kept


# ---------------------------------------------------------------------------
# polyC warning flags
# ---------------------------------------------------------------------------


def _variant_flag_key(v: dict[str, Any]) -> str:
    """Build the canonical ``pos|ref|alt`` flag key for a variant dict."""
    return f"{normalize_position(v.get('pos', 0))}|{v.get('ref', '')}|{v.get('seq', '')}"


def _polyc_warning_flags(variants: list[dict[str, Any]]) -> dict[str, list[str]]:
    """Produce per-variant polyC warning flag strings.

    Warnings are emitted as ``variant_flags`` strings (the flagging layer is
    the sole authority); the standardized model has no metadata fields.

    Adds, keyed by ``pos|ref|alt``:
    - ``"polyC region {region}: potential length heteroplasmy"`` for variants
      inside a ``POLYC_REGIONS`` region.
    - ``"near polyC region {region}: poor signal quality"`` for variants in a
      ``SKIP_AFTER_POLYC`` range.
    - ``"16189T>C creates polyC stretch"`` for the 16189 T>C special case.

    Args:
        variants: List of variant dicts (carrying polyC metadata set by
            ``_handle_polyc_removal`` where applicable).

    Returns:
        Dict mapping ``pos|ref|alt`` -> list[str] of warning flags.
    """
    polyc_regions = get_settings().regions.POLYC_REGIONS
    skip_after = get_settings().regions.SKIP_AFTER_POLYC

    flags: dict[str, list[str]] = {}
    for v in variants:
        pos_int = pos_base(v.get("pos", 0))
        reasons: list[str] = []

        in_region = _region_for_pos(pos_int, polyc_regions)
        if in_region is not None:
            reasons.append(f"polyC region {in_region}: potential length heteroplasmy")

        after_region = _region_for_pos(pos_int, skip_after)
        if after_region is not None:
            reasons.append(f"near polyC region {after_region}: poor signal quality")

        if pos_int == _HV1_POLYC_CREATING_POS and v.get("ref") == "T" and v.get("seq") == "C":
            reasons.append("16189T>C creates polyC stretch")

        if reasons:
            flags[_variant_flag_key(v)] = reasons

    return flags
