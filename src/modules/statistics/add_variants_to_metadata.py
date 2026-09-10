#!/usr/bin/env python
import json
from argparse import ArgumentParser
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from src.core.variants import format_variants_simplified, normalize_position, pos_base, pos_sort_key

SPECIAL_POSITIONS: set[int] = {16193, 455, 463, 573, 309}
_INTERVAL_PAIR_LEN = 2
_TWO_MEMBERS = 2


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


def load_metadata(file_path: Path) -> pd.DataFrame:
    """
    Load metadata from a TSV file.

    Args:
        file_path: Path to the TSV file

    Returns:
        DataFrame containing the metadata

    Raises:
        FileNotFoundError: If the file doesn't exist
        pd.errors.EmptyDataError: If the file is empty
        pd.errors.ParserError: If the file format is invalid
    """
    try:
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        logger.error(f"Error parsing {file_path}: {e}")
        raise


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


def save_statistics(stats: dict[str, Any], file_path: Path) -> None:
    """
    Save batch statistics to a text file.

    Args:
        stats: Dictionary with batch statistics
        file_path: Path to save the statistics

    Raises:
        PermissionError: If the file can't be written due to permissions
        OSError: If another I/O error occurs
    """
    try:
        with file_path.open("w") as f:
            # Write family statistics if available
            if "family" in stats:
                f.write("Family Statistics:\n")
                f.write(f"  Total Families: {stats['family']['total_families']}\n")
                f.write(f"  Families with 1 member: {stats['family']['single_member']}\n")
                f.write(f"  Families with 2 members: {stats['family']['two_members']}\n")
                f.write(f"  Families with more than 2 members: {stats['family']['more_than_two']}\n")
                f.write("\n")

            f.write("STR Batches:\n")
            f.writelines(f"  {batch}: {count}\n" for batch, count in stats["str"].most_common())

            f.write("\nmtDNA Batches:\n")
            f.writelines(f"  {batch}: {count}\n" for batch, count in stats["mtdna"].most_common())
        logger.success(f"Statistics saved to {file_path}")
    except PermissionError:
        logger.error(f"Permission denied when saving to {file_path}")
        raise
    except OSError as e:
        logger.error(f"Error saving to {file_path}: {e}")
        raise


def log_statistics(stats: dict[str, Any]) -> None:
    """
    Log family statistics to the terminal.

    Args:
        stats: Dictionary with statistics
    """
    # Log family statistics if available
    if "family" in stats:
        logger.info("Family Statistics:")
        logger.info(f"  Total Families: {stats['family']['total_families']}")
        logger.info(f"  Families with 1 member: {stats['family']['single_member']}")
        logger.info(f"  Families with 2 members: {stats['family']['two_members']}")
        logger.info(f"  Families with more than 2 members: {stats['family']['more_than_two']}")


def format_identifier(value: object | None, *, remove_decimal: bool = True) -> str | None:
    """
    Format identifiers consistently by handling different input formats.

    Args:
        value: The identifier value to format
        remove_decimal: Whether to remove decimal part if it's just .0

    Returns:
        Formatted identifier or None if invalid
    """
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return None

    # Convert to string and strip whitespace
    str_val = str(value).strip()

    # If empty string after stripping, return None
    if not str_val:
        return None

    # If we should remove decimal points ending in .0
    if remove_decimal and str_val.endswith(".0"):
        str_val = str_val[:-2]

    # For barcodes, accept both numeric and alphanumeric formats
    # This handles both old numeric barcodes and new TN-prefixed barcodes
    if remove_decimal:
        # Accept if it's all digits OR starts with letters followed by digits (like TN4587612228)
        if str_val.isdigit() or (str_val.isalnum() and any(c.isdigit() for c in str_val)):
            return str_val
        return None

    return str_val


def is_special_position(pos: float | str) -> bool:
    """
    Check if position is among those that should be ignored.

    Uses pos_base() for numeric comparison against SPECIAL_POSITIONS
    and checks insertion variants at position 309.

    Args:
        pos: Position to check (int, float, or str)

    Returns:
        True if this is a special position to be ignored
    """
    base_pos = pos_base(pos)

    # Check if the base position is in our special set
    return base_pos in SPECIAL_POSITIONS


def get_variants_for_sample(sample_id: object | None, variants_data: dict[str, dict]) -> str:
    """
    Extract variant information for a sample from variants data.

    Args:
        sample_id: Sample ID to look up
        variants_data: Dictionary with variants data keyed by sample ID

    Returns:
        Formatted string of variants for the sample
    """
    if not sample_id or pd.isna(sample_id):  # type: ignore[reportGeneralTypeIssues]
        return "No data"

    # Format the sample ID properly before lookup
    formatted_id = str(sample_id).split(".")[0]

    if formatted_id not in variants_data:
        return "No data"

    sample_data = variants_data[formatted_id]
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

        # Skip variants at special positions
        if is_special_position(variant_pos):
            continue

        formatted_variants.append({"pos": variant_pos, "ref": variant["ref"], "seq": variant["seq"]})

    return format_variants_simplified(formatted_variants)


def get_intervals_for_sample(sample_id: object | None, variants_data: dict[str, dict]) -> str:
    """
    Extract analyzed intervals information for a sample from variants data.

    Args:
        sample_id: Sample ID to look up
        variants_data: Dictionary with variants data keyed by sample ID

    Returns:
        Formatted string of analyzed intervals for the sample
    """
    if not sample_id or pd.isna(sample_id):  # type: ignore[reportGeneralTypeIssues]
        return "No data"

    # Format the sample ID properly before lookup
    formatted_id = str(sample_id).split(".")[0]

    if formatted_id not in variants_data:
        return "No data"

    sample_data = variants_data[formatted_id]

    # Check if intervals data exists
    if "intervals" not in sample_data:
        return "No data"

    intervals_data = sample_data["intervals"]
    interval_tuples = []

    # Process each HV region in a consistent order
    for region in ["HV1", "HV2", "HV3"]:
        if intervals_data.get(region):
            for interval in intervals_data[region]:
                if len(interval) == _INTERVAL_PAIR_LEN:
                    start, end = interval
                    interval_tuples.append((start, end))

    # If no intervals found, return "No data"
    if not interval_tuples:
        return "No data"

    # Sort intervals by start position and format them
    interval_tuples.sort(key=lambda x: x[0])
    interval_ranges = [f"{start}-{end}" for start, end in interval_tuples]

    # Join all intervals with spaces
    return " ".join(interval_ranges)


def get_batch_statistics(df: pd.DataFrame) -> dict[str, Counter]:
    """
    Generate statistics about the number of samples in each batch.

    Args:
        df: DataFrame with batch information

    Returns:
        Dictionary with statistics by batch type
    """
    stats = {"str": Counter(), "mtdna": Counter()}

    # Count STR batches
    if "STR_batch" in df.columns:
        str_counts = df["STR_batch"].dropna().value_counts()
        for batch, count in str_counts.items():
            stats["str"][batch] = count

    # Count mtDNA batches
    if "mtDNA_batch" in df.columns:
        mtdna_counts = df["mtDNA_batch"].dropna().value_counts()
        for batch, count in mtdna_counts.items():
            stats["mtdna"][batch] = count

    return stats


def get_family_statistics(df: pd.DataFrame) -> dict[str, int]:
    """
    Generate statistics about the families in the dataset.

    Args:
        df: DataFrame with family information

    Returns:
        Dictionary with family statistics
    """
    stats = {"total_families": 0, "single_member": 0, "two_members": 0, "more_than_two": 0}

    # Check if Family_ID column exists
    if "Family_ID" not in df.columns:
        logger.warning("Family_ID column not found in data. Unable to calculate family statistics.")
        return stats

    # Count members in each family
    family_counts = df["Family_ID"].dropna().value_counts()

    # Calculate statistics
    stats["total_families"] = len(family_counts)
    stats["single_member"] = sum(1 for count in family_counts if count == 1)
    stats["two_members"] = sum(1 for count in family_counts if count == _TWO_MEMBERS)
    stats["more_than_two"] = sum(1 for count in family_counts if count > _TWO_MEMBERS)

    return stats


def process_metadata(metadata_df: pd.DataFrame, variants_data: dict[str, Any]) -> pd.DataFrame:
    """
    Process metadata by formatting columns and adding variant and intervals information.

    Args:
        metadata_df: DataFrame with metadata
        variants_data: Dictionary with variants data

    Returns:
        Processed DataFrame with added variant and intervals information
    """
    # Make a copy to avoid modifying the original
    df = metadata_df.copy()

    # Format barcodes
    for col in ["mtDNA_barcode", "STR_barcode"]:
        if col in df.columns:
            df[col] = df[col].apply(lambda x: format_identifier(x, remove_decimal=True))

    # Format IDs
    for col in ["mtDNA_ID", "STR_ID"]:
        if col in df.columns:
            df[col] = df[col].apply(lambda x: format_identifier(x, remove_decimal=False))

    # Add variants column
    logger.info("Adding variants data to metadata")
    df["Variants"] = df["mtDNA_ID"].apply(
        lambda x: get_variants_for_sample(x, variants_data) if pd.notna(x) else "No ID",
    )

    # Add analyzed intervals column
    logger.info("Adding analyzed intervals data to metadata")
    df["Analyzed_Range"] = df["mtDNA_ID"].apply(
        lambda x: get_intervals_for_sample(x, variants_data) if pd.notna(x) else "No ID",
    )

    # Remove Soldier_Name column if it exists (for privacy)
    if "Soldier_Name" in df.columns:
        logger.info("Removing Soldier_Name column for privacy")
        df = df.drop(columns=["Soldier_Name"])

    # Sort by Family_ID if it exists
    if "Family_ID" in df.columns:
        logger.info("Sorting by Family_ID")
        df = df.sort_values(by="Family_ID")

    return df


def arg_parser() -> ArgumentParser:
    """
    Set up command line argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    parser = ArgumentParser(description="Add variant and analyzed intervals information to metadata TSV file")
    parser.add_argument("-i", "--input-file", type=Path, required=True, help="Path to input TSV file with metadata")
    parser.add_argument("-j", "--json-file", type=Path, required=True, help="Path to JSON file with variants data")
    parser.add_argument(
        "-o",
        "--output-file",
        type=Path,
        required=True,
        help="Path to output TSV file (Excel written alongside as .xlsx)",
    )
    parser.add_argument(
        "-s",
        "--stats-file",
        type=Path,
        default=None,
        help="Path to output stats file with batch statistics (optional)",
    )
    return parser


def main() -> None:
    """
    Main function to process metadata and add variant and intervals information.
    """
    # Parse command line arguments
    args = arg_parser().parse_args()

    try:
        # Load variants data
        logger.info(f"Loading variants data from {args.json_file}")
        variants_data = load_json_data(args.json_file)

        # Load metadata
        logger.info(f"Loading metadata from {args.input_file}")
        metadata_df = load_metadata(args.input_file)

        # Process metadata and add variants
        enriched_df = process_metadata(metadata_df, variants_data)

        # Save results as TSV
        logger.info(f"Saving result to {args.output_file}")
        save_dataframe(enriched_df, args.output_file)

        # Always save an Excel copy alongside the TSV (same directory, same stem, .xlsx)
        excel_path = args.output_file.with_suffix(".xlsx")
        logger.info(f"Saving Excel result to {excel_path}")
        save_dataframe(enriched_df, excel_path, excel=True)

        # Calculate statistics
        logger.info("Calculating family statistics")
        stats = {"str": Counter(), "mtdna": Counter(), "family": {}}

        # Get batch statistics for file output only
        if args.stats_file:
            batch_stats = get_batch_statistics(enriched_df)
            stats["str"] = batch_stats["str"]
            stats["mtdna"] = batch_stats["mtdna"]

        # Get family statistics
        stats["family"] = get_family_statistics(enriched_df)

        # Log family statistics to terminal
        log_statistics(stats)

        # Save all statistics to file if requested
        if args.stats_file:
            logger.info(f"Saving statistics to {args.stats_file}")
            save_statistics(stats, args.stats_file)

        logger.success("Variants and intervals added successfully")

    except (FileNotFoundError, json.JSONDecodeError, OSError, ValueError, KeyError, pd.errors.ParserError) as e:
        logger.exception(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
