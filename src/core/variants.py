"""
Variant and position domain utilities for mtDNA analysis.

This module centralizes all variant-related and position-related functions
used across the mtDNA pipeline. Other modules should import from here
for variant/position operations.

Position contract:
- Base positions are int (e.g., 309, 16569)
- Insertion positions are str (e.g., "309.1", "217.10")
- normalize_position() returns int for bases, str for insertions
- Use pos_base() to extract integer base position for comparisons
- Use pos_sort_key() for sorting positions
"""

from typing import Any

from loguru import logger

from src.config import get_settings

# ---------------------------------------------------------------------------
# Position type and helpers
# ---------------------------------------------------------------------------
# Position type: canonical source is src.core.models.
# Re-exported here so callers can import Position/validate_position from
# src.core.variants; declared in __all__ to mark them as intentional exports.
from src.core.models import Position, Variant, validate_position

__all__ = ["Position", "validate_position"]


def pos_base(pos: float | str) -> int:
    """Extract the integer base position from any Position value.

    Base positions are int, insertion positions are str.
    Use this function whenever you need to compare positions numerically.

    Args:
        pos: A Position value (int for base, str for insertion like "309.1").
             Also accepts int/float for backward compatibility.

    Returns:
        Integer base position (e.g., pos_base(309) -> 309, pos_base("309.1") -> 309)
    """
    if isinstance(pos, int):
        return pos
    if isinstance(pos, float):
        return int(pos)
    return int(pos.split(".")[0])


def normalize_position(pos: float | str) -> int | str:
    """Normalize a genomic position to its canonical type.

    Returns int for base positions and str for insertion positions.
    Base positions (no decimal) become int (e.g., 309, 16569).
    Insertions (with decimal) stay str to preserve precision (e.g., "309.1", "217.10").

    Args:
        pos: Position as int, float, or string. Strings may include commas
             and whitespace (e.g., "1,234.10").

    Returns:
        int for base positions, str for insertion positions.

    Examples:
        >>> normalize_position(309)
        309
        >>> normalize_position(309.0)
        309
        >>> normalize_position(309.1)
        '309.1'
        >>> normalize_position("217.10")
        '217.10'
        >>> normalize_position("1,234.10")
        '1234.10'
        >>> normalize_position("309")
        309
        >>> normalize_position(" 309 ")
        309
        >>> normalize_position("309.0")
        309
        >>> normalize_position(0)
        0
        >>> normalize_position(-1)
        -1
    """
    # Handle integer input - already canonical
    if isinstance(pos, int):
        return pos

    # Normalize floats by converting to their string form and falling through to
    # the string handler below. This keeps one code path for decimal positions
    # while preserving the Position contract (int for bases, str for insertions).
    # Note: float(217.10) == 217.1, so the ".10" suffix cannot be recovered from
    # a float; we only get "217.1" via str().
    if isinstance(pos, float):
        pos = str(pos)

    if not isinstance(pos, str):
        msg = f"Cannot normalize position of type {type(pos).__name__}: {pos!r}"
        raise TypeError(msg)

    # Strip whitespace and commas
    cleaned = pos.strip().replace(",", "")

    # Try integer - base position -> return int
    try:
        return int(cleaned)
    except ValueError:
        pass

    # Try float to check if it is a numeric value with a decimal point
    try:
        float_val = float(cleaned)
    except ValueError:
        # Not a number - return the cleaned string as-is
        return cleaned

    # If float is integer-valued (e.g., 309.0), return as int
    if float_val.is_integer():
        return int(float_val)

    # Decimal position (insertion) - return str to preserve precision.
    # "217.10" stays "217.10", "309.1" stays "309.1".
    return cleaned


def pos_sort_key(pos: float | str) -> tuple[int, int]:
    """Generate a sort key for genomic positions.

    Sorts positions by base position first, then by insertion index.
    Base positions sort before their insertions (base -1 < insertion 1).

    Args:
        pos: Position as int, float, or string

    Returns:
        Tuple of (base_position, insertion_index) where insertion_index is
        -1 for base positions and >= 1 for insertion positions.

    Examples:
        >>> pos_sort_key(309)
        (309, -1)
        >>> pos_sort_key("309.1")
        (309, 1)
        >>> pos_sort_key("217.10")
        (217, 10)
        >>> pos_sort_key("217.9")
        (217, 9)
        >>> pos_sort_key("1,234.5")
        (1234, 5)
        >>> pos_sort_key("  309  ")
        (309, -1)
    """
    # Normalize to Position type first
    normalized = normalize_position(pos)

    # Extract base and insertion index
    base = pos_base(normalized)

    if isinstance(normalized, int):
        return (base, -1)
    # String insertion position like "309.1" or "217.10"
    _left, sep, right = str(normalized).partition(".")
    if sep:
        try:
            insertion_index = int(right)
        except ValueError:
            # Non-numeric insertion - use string comparison as fallback
            return (base, -1)
        return (base, insertion_index)
    # No decimal point - base position
    return (base, -1)


# ---------------------------------------------------------------------------
# Variant classification and statistics
# ---------------------------------------------------------------------------


def split_variant_types(
    variants: list[dict],
) -> tuple[list[dict], list[dict], list[dict]]:
    """
    Split variants by type and normalize position keys.

    Normalizes each variant's ``pos`` key to canonical form
    (int for base positions, str for insertions) before splitting,
    so callers always receive position-normalized dicts.

    Args:
        variants: List of variant dictionaries

    Returns:
        Tuple containing sorted lists of SNPs, insertions, and deletions
    """
    snps = []
    insertions = []
    deletions = []
    for variant in variants:
        variant["pos"] = normalize_position(variant["pos"])
        if variant["ref"] == "-" and variant["seq"] != "-":
            insertions.append(variant)
        elif variant["ref"] != "-" and variant["seq"] == "-":
            deletions.append(variant)
        else:
            snps.append(variant)

    sorted_snps = sorted(snps, key=lambda x: pos_sort_key(x["pos"]))
    sorted_insertions = sorted(insertions, key=lambda x: pos_sort_key(x["pos"]))
    sorted_deletions = sorted(deletions, key=lambda x: pos_sort_key(x["pos"]))

    return sorted_snps, sorted_insertions, sorted_deletions


def statistic_variants(variants: dict, hv_seq_refs: dict, hv_seqs: dict) -> tuple[int, int, int, int, int]:
    """
    Calculate statistics for variants.

    Args:
        variants: Dictionary of variants
        hv_seq_refs: Dictionary of HV reference sequences
        hv_seqs: Dictionary of HV sequences (may include insertions, making it
                 longer than hv_seq_refs)

    Returns:
        Tuple containing counts of SNPs, insertions, deletions, identicals, and unreads
    """
    no_snps = 0
    no_ins = 0
    no_dels = 0

    for variant in variants.values():
        if variant["ref"] == "-" and variant["seq"] != "-":
            no_ins += 1
        elif variant["ref"] != "-" and variant["seq"] == "-":
            no_dels += 1
        elif variant["ref"] != variant["seq"]:
            no_snps += 1

    no_identicals = 0
    no_unread = 0

    for hv_name in hv_seq_refs:
        ref_seq = hv_seq_refs[hv_name]["seq"]
        con_seq = hv_seqs[hv_name]["seq"]
        # The consensus sequence may be longer than the reference when insertions
        # are present. Only compare up to the reference length; extra bases in
        # the consensus (insertions) count as differences.
        min_len = min(len(ref_seq), len(con_seq))
        for idx in range(min_len):
            if con_seq[idx] == "N":
                no_unread += 1
            elif con_seq[idx] == ref_seq[idx]:
                no_identicals += 1
            # Bases that differ are already counted in no_snps above
        # Extra consensus bases beyond reference length are insertions — already
        # counted via the variant loop above

    return no_snps, no_ins, no_dels, no_identicals, no_unread


# ---------------------------------------------------------------------------
# IUPAC nucleotide encoding
# ---------------------------------------------------------------------------


def nucleotide_to_iupac(bases: list[str]) -> str:
    """
    Convert nucleotides to IUPAC code.

    Args:
        bases: List of bases

    Returns:
        IUPAC code
    """
    conveted_bases = []
    for base in bases:
        real_bases = iupac_to_bases(base)
        conveted_bases.extend(real_bases)

    iupac_codes = {
        "A": "A",
        "T": "T",
        "C": "C",
        "G": "G",
        "AT": "W",
        "AG": "R",
        "AC": "M",
        "GT": "K",
        "CT": "Y",
        "CG": "S",
        "AGT": "D",
        "ACT": "H",
        "ACG": "V",
        "CGT": "B",
    }
    conveted_bases = sorted(set(conveted_bases))  # remove duplicates and sort
    key = "".join(conveted_bases)
    return iupac_codes.get(key, "N")  # return 'N' if no match


def iupac_to_bases(iupac: str) -> list[str]:
    """
    Convert IUPAC code to list of bases.

    Args:
        iupac: IUPAC code

    Returns:
        List of bases
    """
    iupac_map = {
        "A": ["A"],
        "T": ["T"],
        "C": ["C"],
        "G": ["G"],
        "W": ["A", "T"],
        "R": ["A", "G"],
        "M": ["A", "C"],
        "K": ["T", "G"],
        "Y": ["T", "C"],
        "S": ["C", "G"],
        "D": ["A", "T", "G"],
        "H": ["A", "T", "C"],
        "V": ["A", "C", "G"],
        "B": ["T", "C", "G"],
        "N": ["A", "T", "C", "G"],
    }
    return iupac_map.get(iupac.upper(), [])


def iupac_bases_compatible(left: str, right: str) -> bool:
    """Return whether two valid IUPAC nucleotide calls share a possible base."""
    left_bases = set(iupac_to_bases(left))
    right_bases = set(iupac_to_bases(right))
    return bool(left_bases and right_bases and left_bases & right_bases)


# ---------------------------------------------------------------------------
# Variant formatting
# ---------------------------------------------------------------------------


def format_variants_simplified(variants_list: list[dict[str, Any]]) -> str:
    """
    Format variant list in simplified format for display.

    Args:
        variants_list: List of variant dictionaries with 'pos', 'ref', and 'seq' keys

    Returns:
        Space-separated string of formatted variants (e.g., "73G 263G 315.1C")
        Returns "None" if the list is empty

    Examples:
        >>> variants = [{'pos': '73', 'ref': 'A', 'seq': 'G'}, {'pos': '263', 'ref': 'A', 'seq': 'G'}]
        >>> format_variants_simplified(variants)
        '73G 263G'
    """
    if not variants_list:
        return "None"

    variant_strs = []
    sorted_variants = sorted(
        variants_list,
        key=lambda x: pos_sort_key(x.get("pos", 0)),
    )

    for v in sorted_variants:
        ref = v.get("ref", "")
        seq = v.get("seq", "")
        pos = v.get("pos", 0)

        if ref == "-":  # Insertion
            variant_strs.append(f"{pos}{seq}")
        elif seq == "-":  # Deletion
            variant_strs.append(f"{pos}DEL")
        else:  # SNP
            variant_strs.append(f"{pos}{seq}")

    return " ".join(variant_strs)


def format_variants_list(variants_list: list[str], separator: str = " ", empty_value: str = "None") -> str:
    """
    Format variant list for display in TSV or other outputs.

    Args:
        variants_list: List of variant strings
        separator: Separator to use between variants ("space" or "comma")
        empty_value: Value to return if list is empty

    Returns:
        Formatted string of variants

    Examples:
        >>> format_variants_list(["73G", "263G", "315.1C"])
        '73G 263G 315.1C'
        >>> format_variants_list(["73G", "263G"], separator="comma")
        '73G, 263G'
        >>> format_variants_list([])
        'None'
    """
    if not variants_list:
        return empty_value

    if separator == "space":
        return " ".join(variants_list)
    if separator == "comma":
        return ", ".join(variants_list)
    return separator.join(variants_list)


def parse_variant_position_allele(variant: str) -> tuple[str, str]:
    """
    Parse variant string to extract position and allele.

    Args:
        variant: Variant string like "152C", "315.1C", "524DEL"

    Returns:
        Tuple of (position, allele)

    Examples:
        >>> parse_variant_position_allele("152C")
        ('152', 'C')
        >>> parse_variant_position_allele("315.1C")
        ('315.1', 'C')
        >>> parse_variant_position_allele("524DEL")
        ('524', 'DEL')
    """
    # Handle deletion format
    if variant.endswith("DEL"):
        pos = variant[:-3]
        allele = "DEL"
    else:
        # Find where the position ends and allele begins
        pos = ""
        allele = ""
        for i, char in enumerate(variant):
            if char.isalpha():
                pos = variant[:i]
                allele = variant[i:]
                break

        if not pos:  # fallback
            pos = variant
            allele = ""

    return pos, allele


# ---------------------------------------------------------------------------
# Variant/position region classification
# ---------------------------------------------------------------------------


def get_hv_region_for_position(pos: int | str, regions: dict | None = None) -> str | None:
    """
    Get HV region for a given position.

    Args:
        pos: Genomic position (int or str for insertion positions)
        regions: Optional custom regions dictionary. If not provided, uses get_settings().regions.REGIONS.
                Can be in format:
                - {'HV1': [start, end], 'HV2': [start, end]} (flat format)
                - {'HV1': [[start, end]], 'HV2': [[start, end]]} (nested format)

    Returns:
        HV region name or None if position is not in any region
    """
    # Normalize position for comparison
    pos_int = pos_base(pos)

    # Use default regions if none provided
    if regions is None:
        regions = get_settings().regions.REGIONS

    for hv_region, region_ranges in regions.items():
        # Handle both formats: [start, end] and [[start, end]]
        if isinstance(region_ranges[0], list):
            # Nested format: [[start, end], [start2, end2], ...]
            for start, end in region_ranges:
                if start <= pos_int <= end:
                    return hv_region
        else:
            # Flat format: [start, end]
            start, end = region_ranges[0], region_ranges[1]
            if start <= pos_int <= end:
                return hv_region
    return None


# ---------------------------------------------------------------------------
# Interval math (variant-region context)
# ---------------------------------------------------------------------------


def get_overlap(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    """
    Get overlap between two lists of intervals.

    Args:
        a: First list of intervals
        b: Second list of intervals

    Returns:
        List of overlapping intervals
    """
    return [
        [max(first[0], second[0]), min(first[1], second[1])]
        for first in a
        for second in b
        if max(first[0], second[0]) <= min(first[1], second[1])
    ]


def is_position_in_intervals(pos: int | str, intervals: list[tuple[int, int]]) -> bool:
    """
    Check if a variant position is within any of the provided intervals.

    Args:
        pos: Position of the variant (int or str, may include decimal like "573.1")
        intervals: List of inclusive (start, end) reference-base intervals

    Returns:
        True if the position is within any interval, False otherwise. An insertion
        ``n.x`` is after reference base ``n`` and therefore requires ``n < end``.

    Examples:
        >>> is_position_in_intervals("152", [(73, 340), (16024, 16365)])
        True
        >>> is_position_in_intervals("16193.1", [(16024, 16193)])
        False
        >>> is_position_in_intervals("500", [(73, 340), (16024, 16365)])
        False
    """
    if not intervals:
        return False

    try:
        normalized = normalize_position(pos)
        base_pos = pos_base(normalized)
        is_insertion = isinstance(normalized, str) and "." in normalized

        # Intervals describe covered reference bases. An insertion n.x is
        # located after base n, so an interval ending at n does not cover it.
        if is_insertion:
            return any(start <= base_pos < end for start, end in intervals)
        return any(start <= base_pos <= end for start, end in intervals)
    except (ValueError, AttributeError):
        # If position cannot be parsed, assume it's invalid
        logger.warning(f"Could not parse position: {pos}")
        return False


SPECIAL_POSITIONS = ["455", "463", "573", "309"]


def is_special_position(
    pos: str | float,
    ref: str | None = None,
    seq: str | None = None,
) -> bool:
    """Check if position should be ignored in concordance determination."""
    pos_str = str(pos)

    # Check if this is a special position or its variant (e.g., 573.1, 309.1)
    is_special = pos_str in SPECIAL_POSITIONS or any(pos_str.startswith(f"{sp}.") for sp in SPECIAL_POSITIONS)

    if not is_special:
        return False

    # If we have ref/seq information, only filter INDELs
    if ref is not None and seq is not None:
        # This is an INDEL if either ref or seq is "-" or if they're different length
        return ref == "-" or seq == "-" or len(ref) != len(seq)

    return True


def variant_to_dict(variant: Variant | dict[str, Any]) -> dict[str, Any] | None:
    """Convert a Variant model object to canonical simple dict format.

    The canonical format uses simple keys matching the statistic_fullbatch.json
    contract: {pos, ref, seq, file, quality, peaks}.

    Args:
        variant: A Variant model object or dict with position/ref/alt or pos/ref/seq keys.

    Returns:
        Dict with keys: pos, ref, seq, file, quality, peaks (if present)
    """

    if isinstance(variant, Variant):
        result: dict[str, Any] = {
            "pos": normalize_position(variant.pos),
            "ref": variant.ref,
            "seq": variant.seq,
            "file": list(variant.files) if variant.files else [],
        }
        if variant.peaks is not None:
            result["peaks"] = variant.peaks

        result["quality"] = list(variant.quality) if variant.quality else []

        return result

    # Already a dict — normalize pos and ensure correct keys
    d = dict(variant)
    raw_pos = d.get("pos", d.get("position"))

    if raw_pos is None:
        logger.warning(f"Skipping variant with None position: {variant}")
        return None
    d["pos"] = normalize_position(raw_pos)

    # Map alt → seq if needed
    if "alt" in d and "seq" not in d:
        d["seq"] = d.pop("alt")

    # Ensure required keys
    d.setdefault("file", [])
    d.setdefault("quality", [])

    # Keep only canonical keys (plus peaks and quality if present).
    result = {
        "pos": d["pos"],
        "ref": d["ref"],
        "seq": d["seq"],
        "file": d["file"],
    }
    if "peaks" in d and d["peaks"] is not None:
        result["peaks"] = d["peaks"]

    result["quality"] = d.get("quality", [])

    return result
