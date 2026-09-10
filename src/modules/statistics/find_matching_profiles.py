#!/usr/bin/env python3
"""Find samples sharing the same mtDNA profile.

Given a sample ID or a TXT file containing one sample ID per line and a
merged JSON file (produced by the merge step, e.g. ``merge_json_files.py``),
return every sample ID whose mtDNA profile matches each queried sample.

A "profile" is defined canonically as the sorted set of variants (SNPs,
insertions, and deletions) expressed as ``(pos, ref, seq)`` tuples — the
forensic haplotype. Two samples whose variant sets are identical share the
same mtDNA profile, even if their HV1/HV2/HV3 sequence strings differ due to
interval coverage or read gaps.

Example usage::

    # CLI
    python -m src.modules.statistics.find_matching_profiles sample_ids.txt \
        -i "/mnt/nas/bca/STR/mtDNA/science/results/merged/merged_statistics.json" \
        -o matching_profiles.tsv \
        --no-self \
        --show-profile

    # Programmatic
    from src.modules.statistics.find_matching_profiles import find_matching_samples

    matching = find_matching_samples("IF0040", merged_data)
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, TextIO

from loguru import logger

from src.core.variants import format_variants_simplified, normalize_position, pos_sort_key


def load_merged_json(file_path: str | Path) -> dict[str, Any]:
    """Load a merged sample JSON file keyed by sample ID.

    Args:
        file_path: Path to the merged JSON file.

    Returns:
        Dictionary mapping sample ID to per-sample data.

    Raises:
        FileNotFoundError: If the file does not exist.
        json.JSONDecodeError: If the file contains invalid JSON.
    """
    path = Path(file_path)
    try:
        with path.open(encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        logger.error(f"Merged JSON file not found: {path}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in {path}: {e}")
        raise


def load_sample_ids(sample_id_or_txt: str | Path) -> list[str]:
    """Load query sample IDs from a literal ID or a TXT file.

    TXT input is interpreted as one sample ID per non-empty line.
    """
    value = Path(sample_id_or_txt)
    if value.is_file():
        with value.open(encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]
    if value.suffix.lower() == ".txt":
        raise FileNotFoundError(value)
    return [str(sample_id_or_txt)]


def _variant_entries(sample_data: dict[str, Any]) -> list[dict[str, Any]]:
    """Collect all variant entries (snps + insertions + deletions) from a sample."""
    variants = sample_data.get("variants", {})
    entries: list[dict[str, Any]] = []
    for key in ("snps", "insertions", "deletions"):
        entries.extend(variants.get(key, []))
    return entries


def _format_intervals(sample_data: dict[str, Any]) -> str:
    """Flatten and format the sample intervals as a JSON-like list string."""
    intervals = sample_data.get("intervals", {})
    flat = sorted(
        (iv for ivs in intervals.values() for iv in ivs),
        key=lambda iv: (iv[0], iv[1]),
    )
    return str(flat)


def build_variant_signature(sample_data: dict[str, Any]) -> tuple[tuple[Any, ...], ...]:
    """Build a canonical, hashable signature from a sample's variant set.

    Each variant is reduced to ``(pos, ref, seq)`` and the full set is sorted
    by position (insertion positions sort after their base position) using
    ``pos_sort_key``. Positions are normalized via ``normalize_position`` so
    the mixed ``"309.1"`` (str) and ``309.1`` (float) forms present in the
    merged JSON compare equal.

    Args:
        sample_data: Per-sample dictionary from the merged JSON.

    Returns:
        Sorted tuple of ``(pos, ref, seq)`` tuples. An empty variant set
        yields an empty tuple (still a valid, comparable signature).
    """
    entries = _variant_entries(sample_data)
    reduced = [(normalize_position(e["pos"]) if "pos" in e else 0, e.get("ref"), e.get("seq")) for e in entries]
    reduced.sort(key=lambda v: pos_sort_key(v[0]) if v[0] is not None else (0, 0))
    return tuple(reduced)


def build_profile_index(data: dict[str, Any]) -> dict[tuple[Any, ...], list[str]]:
    """Group all sample IDs by their mtDNA profile signature.

    Args:
        data: Merged JSON keyed by sample ID.

    Returns:
        Dictionary mapping profile signature to a sorted list of sample IDs.
    """
    index: dict[tuple[Any, ...], list[str]] = {}
    for sample_id, sample_data in data.items():
        signature = build_variant_signature(sample_data)
        index.setdefault(signature, []).append(sample_id)

    for sample_ids in index.values():
        sample_ids.sort()
    return index


def find_matching_samples(
    sample_id: str,
    data: dict[str, Any],
    *,
    include_self: bool = True,
) -> list[str]:
    """Return all sample IDs sharing the queried sample's mtDNA profile.

    Args:
        sample_id: The query sample ID.
        data: Merged JSON keyed by sample ID.
        include_self: If True (default), the query sample ID is included in
            the result. If False, it is omitted.

    Returns:
        Sorted list of matching sample IDs. Always non-empty when
        ``include_self`` is True and the sample exists.

    Raises:
        KeyError: If ``sample_id`` is not present in ``data``.
    """
    if sample_id not in data:
        msg = f"Sample ID {sample_id!r} not found in merged data ({len(data)} samples)"
        raise KeyError(msg)

    index = build_profile_index(data)
    signature = build_variant_signature(data[sample_id])

    matches = list(index.get(signature, []))
    if not include_self:
        matches = [sid for sid in matches if sid != sample_id]
    matches.sort()
    return matches


TSV_FIELDS = ["Sample ID", "Intervals", "Profile", "Matched Sample ID"]


def find_matching_profile_rows(
    sample_ids: list[str],
    data: dict[str, Any],
    *,
    include_self: bool = True,
) -> list[dict[str, str]]:
    """Return one TSV-ready row per queried sample with comma-separated matches."""
    index = build_profile_index(data)
    rows: list[dict[str, str]] = []

    for sample_id in sample_ids:
        if sample_id not in data:
            msg = f"Sample ID {sample_id!r} not found in merged data ({len(data)} samples)"
            raise KeyError(msg)

        signature = build_variant_signature(data[sample_id])
        matches = list(index.get(signature, []))
        if not include_self:
            matches = [sid for sid in matches if sid != sample_id]
        matches.sort()

        profile = format_variants_simplified(_variant_entries(data[sample_id]))
        rows.append(
            {
                "Sample ID": sample_id,
                "Intervals": _format_intervals(data[sample_id]),
                "Profile": profile,
                "Matched Sample ID": ", ".join(matches),
            },
        )

    return rows


def write_matching_profile_tsv(rows: list[dict[str, str]], output: str | Path | TextIO | None = None) -> None:
    """Write matching profile rows as TSV to a path or stdout."""
    close_fh = False
    if output is None:
        fh = sys.stdout
    elif isinstance(output, str | Path):
        fh = Path(output).open("w", newline="", encoding="utf-8")
        close_fh = True
    else:
        fh = output

    try:
        writer = csv.DictWriter(fh, fieldnames=TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if close_fh:
            fh.close()


def _arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Find all sample IDs sharing the same mtDNA profile as the queried sample.",
    )
    parser.add_argument("sample_id", help="Query sample ID or TXT file with one sample ID per line")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to merged JSON file keyed by sample ID",
    )
    parser.add_argument(
        "--no-self",
        action="store_true",
        help="Exclude the queried sample ID from the result",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV path. Defaults to stdout.",
    )
    parser.add_argument(
        "--show-profile",
        action="store_true",
        help="Log each queried sample's profile in addition to writing TSV output",
    )
    return parser


def main() -> None:
    """CLI entry point: write matching sample IDs for one or more query samples."""
    parser = _arg_parser()
    args = parser.parse_args()

    try:
        data = load_merged_json(args.input)
        sample_ids = load_sample_ids(args.sample_id)
    except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
        logger.error(f"Failed to load input: {e}")
        sys.exit(1)

    logger.info(f"Loaded {len(data)} samples from {args.input}")
    logger.info(f"Loaded {len(sample_ids)} query sample IDs")

    try:
        rows = find_matching_profile_rows(
            sample_ids,
            data,
            include_self=not args.no_self,
        )
    except KeyError as e:
        logger.error(str(e))
        sys.exit(1)

    if args.show_profile:
        for sample_id in sample_ids:
            profile = format_variants_simplified(_variant_entries(data[sample_id]))
            logger.info(f"{sample_id} profile: {profile}")

    write_matching_profile_tsv(rows, args.output)
    if args.output:
        logger.success(f"Wrote {len(rows)} matching profile rows to {args.output}")


if __name__ == "__main__":
    main()
