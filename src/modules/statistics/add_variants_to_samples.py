#!/usr/bin/env python
"""
Module for adding variant information to blind sample data.

This script processes a sample mapping TSV file and adds variant information
from JSON files to create enhanced sample metadata.
"""

import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from src.config import get_settings
from src.core.variants import (
    format_variants_simplified,
    is_position_in_intervals,
    normalize_position,
    pos_sort_key,
)
from src.modules.statistics.add_variants_to_metadata import get_intervals_for_sample

REGIONS = get_settings().regions.REGIONS


def load_json_data(file_path: Path) -> dict[str, Any]:
    """
    Load data from a JSON file.

    Args:
        file_path: Path to the JSON file

    Returns:
        Dictionary containing the loaded JSON data

    Raises:
        FileNotFoundError: If the file doesn't exist
        json.JSONDecodeError: If the file contains invalid JSON
    """
    try:
        with file_path.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in {file_path}: {e}")
        raise


def load_sample_data(file_path: Path) -> pd.DataFrame:
    """
    Load sample mapping data from a TSV file.

    Args:
        file_path: Path to the TSV file

    Returns:
        DataFrame containing the sample mapping data with proper column names

    Raises:
        FileNotFoundError: If the file doesn't exist
        pd.errors.EmptyDataError: If the file is empty
        pd.errors.ParserError: If the file format is invalid
    """
    try:
        # Load identifiers as strings so they match the string keys in JSON.
        df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["Sample_ID", "Barcode"],
            dtype=str,
        )
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        logger.error(f"Error parsing {file_path}: {e}")
        raise
    else:
        logger.info(f"Loaded {len(df)} samples from {file_path}")
        return df


def create_df_from_json(variants_data: dict[str, Any]) -> pd.DataFrame:
    """
    Create a DataFrame from JSON variants data when no input TSV file is provided.

    Args:
        variants_data: Dictionary with variants data keyed by sample ID

    Returns:
        DataFrame with Sample_ID column and placeholder Barcode column
    """
    sample_ids = list(variants_data.keys())
    df = pd.DataFrame(
        {
            "Sample_ID": sample_ids,
            "Barcode": ["N/A"] * len(sample_ids),  # Placeholder since no barcode mapping available
        },
    )
    logger.info(f"Created DataFrame from JSON data with {len(df)} samples")
    return df


def save_dataframe(df: pd.DataFrame, file_path: Path, *, excel: bool = False) -> None:
    """
    Save DataFrame to a file.

    Args:
        df: DataFrame to save
        file_path: Path to save the file
        excel: Whether to save as Excel format (True) or TSV (False)

    Raises:
        PermissionError: If the file can't be written due to permissions
        OSError: If another I/O error occurs
    """
    try:
        if excel:
            df.to_excel(file_path, index=False)
        else:
            df.to_csv(file_path, sep="\t", index=False)
        logger.success(f"Successfully saved data to {file_path}")
    except PermissionError:
        logger.error(f"Permission denied when saving to {file_path}")
        raise
    except OSError as e:
        logger.error(f"Error saving to {file_path}: {e}")
        raise


def format_variants_standard(variants_list: list[dict[str, Any]]) -> str:
    """
    Format variant list in standard format (space-separated, position+base).

    Delegates the per-variant formatting to ``src.core.variants.format_variants_simplified``
    (canonical pipeline formatter) while preserving the module-specific empty-list
    sentinel ``"No variants"`` (the core helper returns ``"None"``).

    Args:
        variants_list: List of variant dictionaries

    Returns:
        Formatted string of variants in standard format
    """
    if not variants_list:
        return "No variants"
    return format_variants_simplified(variants_list)


def get_variants_for_region(
    variants_list: list[dict[str, Any]],
    region_start: int,
    region_end: int,
) -> list[dict[str, Any]]:
    """
    Filter variants that fall within a specific genomic region.

    Args:
        variants_list: List of variant dictionaries
        region_start: Start of the region (inclusive)
        region_end: End of the region (inclusive)

    Returns:
        List of variants within the specified region
    """
    region_variants = []
    for variant in variants_list:
        if is_position_in_intervals(variant["pos"], [(region_start, region_end)]):
            variant_pos = normalize_position(variant["pos"])
            region_variants.append({"pos": variant_pos, "ref": variant["ref"], "seq": variant["seq"]})
    return region_variants


def get_variants_for_sample(sample_id: str, variants_data: dict[str, dict]) -> str:
    """
    Extract variant information for a sample from variants data.

    Args:
        sample_id: Sample ID to look up
        variants_data: Dictionary with variants data keyed by sample ID

    Returns:
        Formatted string of variants for the sample
    """
    if not sample_id or pd.isna(sample_id):
        return "No data"

    # Use sample_id directly as it should match the JSON keys
    if sample_id not in variants_data:
        logger.warning(f"Sample {sample_id} not found in variants data")
        return "No data"

    sample_data = variants_data[sample_id]
    variants = []

    # Combine SNPs, insertions, and deletions
    if "variants" in sample_data:
        for variant_type in ["snps", "insertions", "deletions"]:
            if variant_type in sample_data["variants"]:
                variants.extend(sample_data["variants"][variant_type])

    # Format variants for easier reading
    formatted_variants = []
    for variant in sorted(variants, key=lambda x: pos_sort_key(x["pos"])):
        variant_pos = normalize_position(variant["pos"])
        formatted_variants.append({"pos": variant_pos, "ref": variant["ref"], "seq": variant["seq"]})

    return format_variants_standard(formatted_variants)


def get_region_variants_for_sample(sample_id: str, variants_data: dict[str, dict], region_name: str) -> str:
    """
    Extract variant information for a sample from variants data within a specific region.

    Args:
        sample_id: Sample ID to look up
        variants_data: Dictionary with variants data keyed by sample ID
        region_name: Name of the region (HV1, HV2, or HV3)

    Returns:
        Formatted string of variants for the sample within the specified region
    """
    if not sample_id or pd.isna(sample_id):
        return "No data"

    if sample_id not in variants_data:
        return "No data"

    if region_name not in REGIONS:
        return "No data"

    region_start, region_end = REGIONS[region_name]
    sample_data = variants_data[sample_id]
    variants = []

    # Combine SNPs, insertions, and deletions
    if "variants" in sample_data:
        for variant_type in ["snps", "insertions", "deletions"]:
            if variant_type in sample_data["variants"]:
                variants.extend(sample_data["variants"][variant_type])

    # Filter variants for the specific region
    region_variants = get_variants_for_region(variants, region_start, region_end)

    return format_variants_standard(region_variants)


def get_sample_statistics(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate statistics about the samples.

    Args:
        df: DataFrame with sample information

    Returns:
        Dictionary with sample statistics
    """
    stats = {"total_samples": len(df), "samples_with_variants": 0, "samples_without_variants": 0}

    # Count samples with and without variants
    if "Variants" in df.columns:
        stats["samples_with_variants"] = len(df[df["Variants"] != "No data"])
        stats["samples_without_variants"] = len(df[df["Variants"] == "No data"])

    # Add statistics for each region
    for region_name in REGIONS:
        if region_name in df.columns:
            no_variants = df[region_name].isin(["No data", "No variants"])
            stats[f"{region_name}_with_variants"] = len(df[~no_variants])
            stats[f"{region_name}_without_variants"] = len(df[no_variants])

    return stats


def log_sample_statistics(stats: dict[str, Any]) -> None:
    """
    Log sample statistics to the terminal.

    Args:
        stats: Dictionary with statistics
    """
    logger.info("Sample Statistics:")
    logger.info(f"  Total samples: {stats['total_samples']}")
    logger.info(f"  Samples with variants: {stats['samples_with_variants']}")
    logger.info(f"  Samples without variants: {stats['samples_without_variants']}")

    # Log region-specific statistics
    for region_name in REGIONS:
        with_key = f"{region_name}_with_variants"
        without_key = f"{region_name}_without_variants"
        if with_key in stats and without_key in stats:
            logger.info(f"  {region_name} with variants: {stats[with_key]}")
            logger.info(f"  {region_name} without variants: {stats[without_key]}")


def process_sample_data(sample_df: pd.DataFrame, variants_data: dict[str, Any]) -> pd.DataFrame:
    """
    Process sample data by adding variant information.

    Args:
        sample_df: DataFrame with sample mapping data
        variants_data: Dictionary with variants data

    Returns:
        Processed DataFrame with added variant information
    """
    # Make a copy to avoid modifying the original
    df = sample_df.copy()

    # Add variants column
    logger.info("Adding variants data to samples")
    df["Variants"] = df["Sample_ID"].apply(lambda x: get_variants_for_sample(x, variants_data))

    # Add region-specific variant columns
    logger.info("Adding region-specific variant columns")
    for region_name in REGIONS:
        logger.info(f"Processing {region_name} variants")
        df[region_name] = df["Sample_ID"].apply(
            lambda x, region_name=region_name: get_region_variants_for_sample(x, variants_data, region_name)
        )

    # Add analyzed intervals column (reuses the shared TNLS interval formatter)
    logger.info("Adding analyzed intervals data to samples")
    df["Analyzed_Range"] = df["Sample_ID"].apply(lambda x: get_intervals_for_sample(x, variants_data))

    return df


def arg_parser() -> ArgumentParser:
    """
    Set up command line argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    parser = ArgumentParser(description="Add variant information to sample mapping TSV or process JSON-only data")
    parser.add_argument(
        "-i",
        "--input-file",
        type=Path,
        required=False,
        default=None,
        help="Path to input sample mapping TSV file (optional - if not provided, will process JSON-only)",
    )
    parser.add_argument(
        "-j",
        "--json-file",
        type=Path,
        required=True,
        help="Path to JSON file with variants data (merged.json)",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=Path,
        required=True,
        help="Path to output TSV file with added variants column (Excel is written alongside with .xlsx suffix)",
    )
    return parser


def main() -> None:
    """
    Main function to process sample data and add variant information.
    """
    # Parse command line arguments
    args = arg_parser().parse_args()

    try:
        # Load variants data first
        logger.info(f"Loading variants data from {args.json_file}")
        variants_data = load_json_data(args.json_file)

        # Load or create sample data based on whether input file is provided
        if args.input_file:
            logger.info(f"Loading sample data from {args.input_file}")
            sample_df = load_sample_data(args.input_file)
        else:
            logger.info("No input TSV file provided, creating DataFrame from JSON data")
            sample_df = create_df_from_json(variants_data)

        # Process the data
        logger.info("Processing sample data and adding variants")
        processed_df = process_sample_data(sample_df, variants_data)

        # Generate and log statistics
        stats = get_sample_statistics(processed_df)
        log_sample_statistics(stats)

        # Save the processed data
        logger.info(f"Saving processed data to {args.output_file}")
        save_dataframe(processed_df, args.output_file)

        # Always save an Excel copy alongside the TSV (same directory, same stem, .xlsx)
        excel_path = args.output_file.with_suffix(".xlsx")
        logger.info(f"Saving Excel file to {excel_path}")
        save_dataframe(processed_df, excel_path, excel=True)

        logger.success("Processing completed successfully!")

    except Exception as e:
        logger.error(f"An error occurred during processing: {e}")
        logger.exception("Traceback:")
        raise


if __name__ == "__main__":
    main()
