#!/usr/bin/env python
"""Low-level utilities for the TNLS/HCLS verification pipeline.

These functions operate on plain data structures (dicts, lists, intervals)
and carry no pipeline state — they are pure and stateless.
"""

import json
from pathlib import Path
from typing import Any

from src.core.variants import format_variants_simplified, iupac_bases_compatible, normalize_position, pos_base

# ── I/O ─────────────────────────────────────────────────────────────────────────


def load_json(path: Path | str) -> dict[str, Any]:
    """Load and return the contents of a JSON file."""
    with Path(path).open(encoding="utf-8") as f:
        return json.load(f)


def export_sample_lists(
    data_dir: Path | str,
    *,
    results_name: str = "results",
) -> list[tuple[Path, Path, int]]:
    """Scan all .json files under ``data_dir`` (recursively, skipping any
    directory whose name matches ``results_name``) and write a sidecar
    ``<name>.json.txt`` alongside each one containing its top-level keys.

    Args:
        data_dir: Root directory to scan.
        results_name: Sub-directory name to exclude (default "results").

    Returns:
        List of ``(json_path, txt_path, key_count)`` tuples — useful for
        logging what was written.

    Example::

        results = export_sample_lists("/data", results_name="results")
        # writes  /data/foo.json  ->  /data/foo.json.txt
        # writes  /data/sub/bar.json  ->  /data/sub/bar.json.txt
    """
    data_dir = Path(data_dir)
    rows: list[tuple[Path, Path, int]] = []

    for json_path in sorted(data_dir.rglob("*.json")):
        if results_name in json_path.parts:
            continue
        txt_path = Path(str(json_path) + ".txt")
        keys = list(load_json(json_path).keys())
        with txt_path.open("w", encoding="utf-8") as out:
            out.write("\n".join(keys) + "\n")
        rows.append((json_path, txt_path, len(keys)))

    return rows


# ── Interval arithmetic ───────────────────────────────────────────────────────


def check_number_in_intervals(intervals: list[list[int]], number: int | str) -> bool:
    """Return True if ``number`` falls inside any of the given intervals.

    Used to test whether a genomic position (e.g. a variant at bp 16129)
    lies within the overlapping region between two samples.
    """
    return any(s <= pos_base(number) <= e for s, e in intervals)


# ── Variant extraction ────────────────────────────────────────────────────────


type VariantPosition = int | str


def get_variant_positions(
    variants: dict[str, list[dict[str, Any]]], overlap_intervals: list[list[int]]
) -> set[VariantPosition]:
    """Return canonical variant positions inside the overlap intervals."""
    return {
        normalize_position(record["pos"])
        for records in variants.values()
        for record in records
        if check_number_in_intervals(overlap_intervals, record["pos"])
    }


def get_hcls_unknown_snp_positions(
    variants: dict[str, list[dict[str, Any]]], overlap_intervals: list[list[int]]
) -> set[VariantPosition]:
    """Return overlapping SNP positions where HCLS reports an unknown ``N`` base."""
    return {
        normalize_position(record["pos"])
        for record in variants.get("snps", [])
        if str(record.get("seq", "")).strip().upper() == "N"
        and check_number_in_intervals(overlap_intervals, record["pos"])
    }


def count_position_mismatches(
    left: set[VariantPosition], right: set[VariantPosition], ignored: set[VariantPosition]
) -> int:
    """Count positions present on only one side after removing ignored positions."""
    return len((left - ignored) ^ (right - ignored))


# ── Variant comparison ─────────────────────────────────────────────────────────


def combined_variants(data: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Flatten all variant types into a single typed-key dict.

    The key prefixes ``snp_``, ``ins_``, ``del_`` prevent silent collisions
    when different variant types share the same genomic position
    (e.g. a SNP and an insertion at the same bp).

    Keys are returned sorted for deterministic output order.

    Example key set::

        snp_16129  -> {"ref": "T", "seq": "C"}
        ins_309    -> {"ref": "-", "seq": "CC"}
        del_249    -> {"ref": "A", "seq": "-"}

    .. warning::
        Keys are typed strings, not bare integers.
        Use ``combined[f"snp_{pos}"]``, not ``combined[pos]``.
    """
    prefixes = {"snps": "snp", "insertions": "ins", "deletions": "del"}
    combined: dict[str, dict[str, str]] = {}
    for variant_type, prefix in prefixes.items():
        for variant in data["variants"].get(variant_type, []):
            position = normalize_position(variant["pos"])
            combined[f"{prefix}_{position}"] = {"ref": variant["ref"], "seq": variant["seq"]}
    return dict(sorted(combined.items()))


def _variant_in_overlap(key: str, ranges: list[list[int]]) -> bool:
    """Return whether a typed variant key falls inside an overlap range."""
    try:
        position = normalize_position(key.split("_", 1)[1])
    except (ValueError, IndexError):
        return False
    return any(start <= pos_base(position) <= end for start, end in ranges)


def _variants_match(key: str, hcls: dict[str, str], tnls: dict[str, str]) -> bool:
    """Return whether two variants have compatible reference and called alleles."""
    if hcls["ref"] != tnls["ref"]:
        return False
    if key.startswith("snp_"):
        return iupac_bases_compatible(hcls["seq"], tnls["seq"])
    return hcls["seq"] == tnls["seq"]


def _is_unknown_hcls_snp(key: str, variant: dict[str, str]) -> bool:
    """Return whether a typed HCLS call is an unresolved SNP base."""
    return key.startswith("snp_") and variant["seq"].strip().upper() == "N"


def compare_variants(
    hcls_variants: dict[str, dict[str, str]],
    tnls_variants: dict[str, dict[str, str]],
    overlap_ranges: list[list[int]],
) -> tuple[int, int, dict[str, Any], dict[str, Any]]:
    """Classify every variant key present in either sample as match or mismatch.

    Variant keys are typed: ``snp_<pos>``, ``ins_<pos>``, ``del_<pos>``.
    Only keys whose genomic position falls within ``overlap_ranges`` are
    evaluated — variants outside the overlapping region are ignored.

    HCLS SNP calls with sequence ``N`` are ignored because they contain no
    resolved base information. For other SNPs, IUPAC calls match when their
    possible-base sets overlap. Insertions and deletions must match exactly.

    For each remaining key in the union of both samples' keys:

        Both samples have a compatible key and allele  ->  MATCH
        Both samples have incompatible alleles         ->  MISMATCH
        Only one sample has the key                    ->  MISMATCH (one-sided)

    Parameters
    ----------
    hcls_variants, tnls_variants : dict
        Output of ``combined_variants()``: mapping of typed variant keys
        to ``{"ref": ..., "seq": ...}``.
    overlap_ranges : list of [start, end]
        The intersected genomic intervals between the two samples.

    Returns
    -------
    tuple
        (n_matches, n_mismatches, match_dict, mismatch_dict)

    Example mismatch_dict entry for a one-sided variant::

        {"base_sample": {"ref": "T", "seq": "C"}, "target_sample": None}
    """
    matches, mismatches = 0, 0
    match_dict: dict[str, Any] = {}
    mismatch_dict: dict[str, Any] = {}

    for key in sorted(hcls_variants.keys() | tnls_variants.keys()):
        if not _variant_in_overlap(key, overlap_ranges):
            continue

        hcls = hcls_variants.get(key)
        tnls = tnls_variants.get(key)
        if hcls is not None and _is_unknown_hcls_snp(key, hcls):
            continue
        if hcls is None or tnls is None:
            mismatches += 1
            mismatch_dict[key] = {"base_sample": hcls, "target_sample": tnls}
        elif _variants_match(key, hcls, tnls):
            matches += 1
            match_dict[key] = hcls
        else:
            mismatches += 1
            mismatch_dict[key] = {"base_sample": hcls, "target_sample": tnls}

    return matches, mismatches, match_dict, mismatch_dict


def get_range_details(ranges: list[list[int]]) -> list[str]:
    """Convert numeric [start, end] intervals to human-readable HV-region labels.

    Only intervals that are *fully contained* within a known HV region are
    labelled.  Partial overlaps are silently dropped.

    Returns a list of strings, one per labelled interval::

        [[16024, 16365]]  ->  ["HV1: 16024-16365"]

    Used by ``compare_TNLS_HCLS_stats`` to produce the ``overlap_details``
    column in the output TSV.
    """
    hv_ranges = {"HV1": (16024, 16365), "HV2": (73, 340), "HV3": (438, 576)}
    labelled = []
    for start, end in ranges:
        for label, (rs, re) in hv_ranges.items():
            if rs <= start <= end <= re:
                labelled.append(f"{label}: {start}-{end}")
                break
    return labelled


# ── Variant formatting (simplified display) ────────────────────────────────────


def _typed_key_to_pos(key: str) -> str:
    """Strip the ``snp_``/``ins_``/``del_`` prefix from a typed variant key."""
    return key.split("_", 1)[1] if "_" in key else key


def variant_match_to_list(variant_match: dict[str, Any]) -> list[dict[str, Any]]:
    """Convert a typed-key ``variant_match`` dict to a flat ``{pos, ref, seq}`` list.

    Keys are typed (``snp_<pos>`` / ``ins_<pos>`` / ``del_<pos>``) mapped to
    ``{ref, seq}``; this rebuilds the ``{pos, ref, seq}`` shape expected by
    :func:`src.core.variants.format_variants_simplified`.
    """
    return [
        {"pos": _typed_key_to_pos(key), "ref": val.get("ref", ""), "seq": val.get("seq", "")}
        for key, val in variant_match.items()
        if isinstance(val, dict)
    ]


def format_matched_variants(variant_match: dict[str, Any]) -> str:
    """Format a ``variant_match`` dict (typed keys) into the simplified display form.

    Returns ``"None"`` when there are no matched variants.
    """
    if not variant_match:
        return "None"
    return format_variants_simplified(variant_match_to_list(variant_match))


def format_sample_variants(sample: dict[str, Any]) -> str:
    """Format a sample's full variant profile (SNPs + insertions + deletions).

    Combines every variant type from ``sample["variants"]`` into the flat
    ``{pos, ref, seq}`` list and formats it via
    :func:`src.core.variants.format_variants_simplified` so the complete profile
    (all regions, not only the overlap) is shown (e.g. ``73G 249DEL 309.1C
    16093C``). Returns ``"None"`` when the sample has no variants.
    """
    variants = sample.get("variants", {})
    combined = [
        *variants.get("snps", []),
        *variants.get("insertions", []),
        *variants.get("deletions", []),
    ]
    return format_variants_simplified(combined)


def format_mismatch_variants(variant_mismatch: dict[str, Any]) -> str:
    """Format a ``variant_mismatch`` dict (typed keys) into a simplified display form.

    Each mismatched site is rendered as ``<pos><base>><target>`` where each
    allele is the seq char, ``DEL`` for a deletion, or ``·`` when that side is
    absent (one-sided mismatch). Entries are space-separated and sorted by
    numeric position. Returns ``"None"`` when there are no mismatches.
    """
    if not variant_mismatch:
        return "None"

    def _allele(side: dict[str, str] | None) -> str:
        if not isinstance(side, dict):
            return "·"
        seq = side.get("seq", "")
        return "DEL" if seq == "-" else str(seq)

    records: list[tuple[int, str]] = []
    for key, val in variant_mismatch.items():
        if not isinstance(val, dict):
            continue
        pos = _typed_key_to_pos(key)
        base = _allele(val.get("base_sample"))  # type: ignore[arg-type]
        target = _allele(val.get("target_sample"))  # type: ignore[arg-type]
        records.append((pos_base(normalize_position(pos)), f"{pos}{base}>{target}"))

    records.sort(key=lambda x: x[0])
    return " ".join(text for _, text in records) if records else "None"
