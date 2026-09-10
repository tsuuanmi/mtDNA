#!/usr/bin/env python3
"""
Variant Frequency Analyzer for mtDNA Data

This script analyzes the merged_statistics.json file to count variant frequencies
across all samples and generates a TSV database of all variants with their frequencies.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

_ZERO = 0
_HUNDRED = 100
_FIFTY = 50
_TEN = 10
_ONE = 1


def load_json_data(file_path: str | Path) -> dict[str, Any]:
    """Load the merged statistics JSON file."""
    try:
        with Path(file_path).open() as f:
            data = json.load(f)
    except FileNotFoundError:
        logger.error(f"File {file_path} not found")
        sys.exit(1)
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON format - {e}")
        sys.exit(1)
    else:
        logger.info(f"Successfully loaded data for {len(data)} samples")
        return data


def extract_variants(data: dict[str, Any]) -> tuple[defaultdict[str, int], int]:
    """Extract all variants from the data and count their frequencies."""
    variant_counts = defaultdict(int)
    total_samples = 0
    samples_with_variants = 0

    for sample_data in data.values():
        total_samples += 1

        # Skip samples that only have batch information
        if "variants" not in sample_data:
            continue

        variants = sample_data.get("variants", {})
        sample_has_variants = False

        # Process SNPs
        for snp in variants.get("snps", []):
            variant_key = f"SNP_{snp['pos']}_{snp['ref']}>{snp['seq']}"
            variant_counts[variant_key] += 1
            sample_has_variants = True

        # Process insertions
        for insertion in variants.get("insertions", []):
            variant_key = f"INS_{insertion['pos']}_{insertion['ref']}>{insertion['seq']}"
            variant_counts[variant_key] += 1
            sample_has_variants = True

        # Process deletions
        for deletion in variants.get("deletions", []):
            variant_key = f"DEL_{deletion['pos']}_{deletion['ref']}>{deletion['seq']}"
            variant_counts[variant_key] += 1
            sample_has_variants = True

        if sample_has_variants:
            samples_with_variants += 1

    logger.info(f"Total samples processed: {total_samples}")
    logger.info(f"Samples with variant data: {samples_with_variants}")
    logger.info(f"Unique variants found: {len(variant_counts)}")

    return variant_counts, samples_with_variants


def create_variant_database(variant_counts: defaultdict[str, int], total_samples_with_variants: int) -> pd.DataFrame:
    """Create a pandas DataFrame with variant information and frequencies."""
    variants_data = []

    for variant_key, count in variant_counts.items():
        # Parse variant key
        parts = variant_key.split("_")
        variant_type = parts[0]
        position = parts[1]
        ref_alt = parts[2].split(">")
        ref_allele = ref_alt[0]
        alt_allele = ref_alt[1]

        # Calculate frequency
        frequency = count / total_samples_with_variants if total_samples_with_variants > _ZERO else 0
        percentage = frequency * 100

        variants_data.append(
            {
                "Variant_Type": variant_type,
                "Position": position,
                "Reference_Allele": ref_allele,
                "Alternative_Allele": alt_allele,
                "Count": count,
                "Frequency": round(frequency, 6),
                "Percentage": round(percentage, 2),
            },
        )

    # Create DataFrame and sort by frequency (descending)
    df = pd.DataFrame(variants_data)
    df = df.sort_values(["Frequency", "Position"], ascending=[False, True])  # type: ignore[reportCallIssue]
    return df.reset_index(drop=True)


def save_results(df: pd.DataFrame, output_file: str | Path) -> None:
    """Save the variant database to a TSV file."""
    try:
        df.to_csv(output_file, sep="\t", index=False)
        logger.success(f"Variant frequency database saved to: {output_file}")
        logger.info(f"Total variants in database: {len(df)}")
    except OSError as e:
        logger.error(f"Error saving file: {e}")
        sys.exit(1)


def print_summary_statistics(df: pd.DataFrame) -> None:
    """Print summary statistics about the variants."""
    logger.info("=" * 60)
    logger.info("VARIANT FREQUENCY ANALYSIS SUMMARY")
    logger.info("=" * 60)

    # Overall statistics
    total_variants = len(df)
    snp_count = len(df[df["Variant_Type"] == "SNP"])
    ins_count = len(df[df["Variant_Type"] == "INS"])
    del_count = len(df[df["Variant_Type"] == "DEL"])

    logger.info(f"Total unique variants: {total_variants}")
    logger.info(f"  - SNPs: {snp_count} ({snp_count / total_variants * _HUNDRED:.1f}%)")
    logger.info(f"  - Insertions: {ins_count} ({ins_count / total_variants * _HUNDRED:.1f}%)")
    logger.info(f"  - Deletions: {del_count} ({del_count / total_variants * _HUNDRED:.1f}%)")

    # Frequency statistics
    logger.info("Frequency Distribution:")
    top = df.iloc[0]
    variant_id = f"{top['Variant_Type']}_{top['Position']}_{top['Reference_Allele']}>{top['Alternative_Allele']}"
    logger.info(f"  - Most common variant: {variant_id} ({df.iloc[0]['Percentage']:.1f}%)")
    logger.info(f"  - Median frequency: {df['Percentage'].median():.2f}%")
    logger.info(f"  - Variants >50%: {len(df[df['Percentage'] > _FIFTY])}")
    logger.info(f"  - Variants >10%: {len(df[df['Percentage'] > _TEN])}")
    logger.info(f"  - Variants >1%: {len(df[df['Percentage'] > _ONE])}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze variant frequencies from merged_statistics.json")
    parser.add_argument(
        "--input",
        "-i",
        default="merged_statistics.json",
        help="Input JSON file (default: merged_statistics.json)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="variant_frequency_database.tsv",
        help="Output TSV file (default: variant_frequency_database.tsv)",
    )

    args = parser.parse_args()

    logger.info("mtDNA Variant Frequency Analyzer")
    logger.info("=" * 40)

    # Load data
    logger.info(f"Loading data from: {args.input}")
    data = load_json_data(args.input)

    # Extract and count variants
    logger.info("Extracting variants")
    variant_counts, total_samples_with_variants = extract_variants(data)

    if not variant_counts:
        logger.warning("No variants found in the data!")
        sys.exit(1)

    # Create variant database
    logger.info("Creating variant frequency database")
    df = create_variant_database(variant_counts, total_samples_with_variants)

    # Save results
    save_results(df, args.output)

    # Print summary
    print_summary_statistics(df)

    logger.success(f"Analysis complete! Results saved to {args.output}")


if __name__ == "__main__":
    main()
