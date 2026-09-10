"""Calculate population-style mtDNA variant frequencies from merged sample profiles."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, TextIO

from loguru import logger

from src.core.variants import format_variants_simplified, normalize_position, pos_sort_key

DEFAULT_INPUT_JSON = Path("/mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json")
DEFAULT_OUTPUT_TSV = Path("results/modules/statistics/variant_frequencies.tsv")
VARIANT_GROUPS = ("snps", "insertions", "deletions")
TSV_FIELDS = ["Variant", "Type", "Pos", "Ref", "Seq", "Sample_Count", "Total_Samples", "Frequency"]
VariantKey = tuple[int | str, str, str]


def load_merged_json(file_path: str | Path) -> dict[str, Any]:
    """Load a merged sample JSON file keyed by sample ID."""
    path = Path(file_path)
    try:
        with path.open(encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        logger.error("Merged JSON file not found: {}", path)
        raise
    except json.JSONDecodeError as e:
        logger.error("Invalid JSON in {}: {}", path, e)
        raise

    if not isinstance(data, dict):
        msg = f"Expected merged JSON object keyed by sample ID, got {type(data).__name__}"
        raise TypeError(msg)
    return data


def variant_type(ref: str, seq: str) -> str:
    """Return the variant class from reference and observed alleles."""
    if ref == "-" and seq != "-":
        return "insertion"
    if ref != "-" and seq == "-":
        return "deletion"
    return "snp"


def iter_variant_keys(sample_data: dict[str, Any]) -> list[VariantKey]:
    """Return canonical variant keys for one sample."""
    variants = sample_data.get("variants", {})
    keys: list[VariantKey] = []
    for group in VARIANT_GROUPS:
        keys.extend(
            (
                normalize_position(variant["pos"]),
                str(variant.get("ref", "")),
                str(variant.get("seq", "")),
            )
            for variant in variants.get(group, [])
        )
    return keys


def count_variant_frequencies(data: dict[str, Any]) -> Counter[VariantKey]:
    """Count how many samples contain each unique variant."""
    counts: Counter[VariantKey] = Counter()
    for sample_data in data.values():
        for key in set(iter_variant_keys(sample_data)):
            counts[key] += 1
    return counts


def build_frequency_rows(data: dict[str, Any]) -> list[dict[str, str]]:
    """Build TSV-ready variant frequency rows from merged sample data."""
    total_samples = len(data)
    counts = count_variant_frequencies(data)
    rows: list[dict[str, str]] = []

    for pos, ref, seq in sorted(
        counts,
        key=lambda key: (-counts[key], pos_sort_key(key[0]), key[1], key[2]),
    ):
        sample_count = counts[(pos, ref, seq)]
        frequency = sample_count / total_samples if total_samples else 0.0
        variant = {"pos": pos, "ref": ref, "seq": seq}
        rows.append(
            {
                "Variant": format_variants_simplified([variant]),
                "Type": variant_type(ref, seq),
                "Pos": str(pos),
                "Ref": ref,
                "Seq": seq,
                "Sample_Count": str(sample_count),
                "Total_Samples": str(total_samples),
                "Frequency": f"{frequency:.8f}",
            },
        )

    return rows


def _write_frequency_tsv_to_file(rows: list[dict[str, str]], fh: TextIO) -> None:
    """Write variant frequency rows to an open file handle."""
    writer = csv.DictWriter(fh, fieldnames=TSV_FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)


def write_frequency_tsv(rows: list[dict[str, str]], output: str | Path | TextIO) -> None:
    """Write variant frequency rows as TSV."""
    if isinstance(output, str | Path):
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as fh:
            _write_frequency_tsv_to_file(rows, fh)
        return

    _write_frequency_tsv_to_file(rows, output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=DEFAULT_INPUT_JSON,
        help=f"Merged statistics JSON keyed by sample ID (default: {DEFAULT_INPUT_JSON})",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_TSV,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT_TSV})",
    )
    return parser.parse_args()


def main() -> None:
    """CLI entry point."""
    args = parse_args()
    try:
        data = load_merged_json(args.input)
        rows = build_frequency_rows(data)
        write_frequency_tsv(rows, args.output)
    except (FileNotFoundError, json.JSONDecodeError, OSError, TypeError, KeyError) as e:
        logger.error("Failed to calculate variant frequencies: {}", e)
        sys.exit(1)

    logger.success(
        "Wrote {} variant frequency rows for {} samples to {}",
        len(rows),
        len(data),
        args.output,
    )


if __name__ == "__main__":
    main()
