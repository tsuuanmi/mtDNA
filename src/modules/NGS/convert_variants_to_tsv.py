#!/usr/bin/env python
"""
Convert variant data from JSON format to simplified TSV format.

This script reads variant data from statistic_fullbatch.json and converts it to a
tab-separated format with Sample_ID and Genotype columns, where genotypes are
formatted as simplified variant strings (e.g., "73G 263G 315.1C").
"""

import json
import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import Any

from loguru import logger

from src.core.variants import format_variants_simplified, pos_sort_key


def load_json_data(json_path: Path) -> dict[str, Any]:
    """
    Load JSON data from file.

    Args:
        json_path: Path to JSON file

    Returns:
        Dictionary containing the JSON data

    Raises:
        FileNotFoundError: If the file doesn't exist
        json.JSONDecodeError: If the file is not valid JSON
    """
    if not json_path.exists():
        msg = f"JSON file not found: {json_path}"
        raise FileNotFoundError(msg)

    try:
        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse JSON file: {e}")
        raise
    else:
        logger.info(f"Successfully loaded JSON data from {json_path}")
        return data


def extract_variants_from_sample(sample_data: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Extract and combine all variants (SNPs, insertions, deletions) from sample data.

    Args:
        sample_data: Dictionary containing sample variant information

    Returns:
        List of variant dictionaries sorted by position
    """

    # Get all variant types
    variant_data = sample_data.get("variants", {})
    snps = variant_data.get("snps", [])
    insertions = variant_data.get("insertions", [])
    deletions = variant_data.get("deletions", [])

    # Combine all variants
    all_variants = snps + insertions + deletions

    # Sort by position (handle decimal positions like 315.1)
    return sorted(all_variants, key=lambda x: pos_sort_key(x.get("pos", 0)))


def convert_json_to_tsv(json_path: Path, output_path: Path) -> None:
    """
    Convert JSON variant data to simplified TSV format.

    Args:
        json_path: Path to input JSON file (statistic_fullbatch.json)
        output_path: Path to output TSV file
    """
    # Load JSON data
    data = load_json_data(json_path)

    # Prepare output data
    output_lines = ["Sample_ID\tGenotype"]

    # Process each sample
    for sample_id, sample_data in sorted(data.items()):
        # Extract all variants
        variants = extract_variants_from_sample(sample_data)

        # Format variants using the simplified format
        genotype = format_variants_simplified(variants)

        # Add to output
        output_lines.append(f"{sample_id}\t{genotype}")

        logger.debug(f"Processed {sample_id}: {len(variants)} variants")

    # Write to TSV file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(output_lines))

    logger.success(f"Successfully converted {len(data)} samples to TSV format: {output_path}")
    logger.info(f"Total samples processed: {len(data)}")


def arg_parser() -> ArgumentParser:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Convert variant JSON data to simplified TSV format")
    parser.add_argument("json_path", type=str, help="Path to input JSON file (statistic_fullbatch.json)")
    parser.add_argument("output_path", type=str, help="Path to output TSV file")
    return parser


def main() -> None:
    """Main entry point for the script."""
    parser = arg_parser()
    args = parser.parse_args()

    json_path = Path(args.json_path)
    output_path = Path(args.output_path)

    try:
        convert_json_to_tsv(json_path, output_path)
    except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
        logger.error(f"Failed to convert JSON to TSV: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
