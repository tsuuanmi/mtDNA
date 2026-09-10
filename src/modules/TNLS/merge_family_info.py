#!/usr/bin/env python3
"""
This script merges information from two TSV files based on mtDNA_barcode or STR_barcode.
It can also process TNLS_with_soldier_names.tsv to create Family and Family_ID columns.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from loguru import logger

_BARCODE_LEN = 12


def is_valid_barcode(value: object | None) -> bool:
    """Check if barcode is valid (all digits and correct length)."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return False
    str_val = str(value).rstrip(".0")
    # Check if the string contains only digits and has correct length
    return str_val.isdigit() and len(str_val) == _BARCODE_LEN


def format_numeric_id(value: object | None) -> object | None:
    """Format numeric IDs by preserving both leading and trailing zeros."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return None

    # Convert to string
    str_val = str(value)

    # If it's a float value (has decimal point)
    if "." in str_val:
        # Check if it's a float that represents an integer (e.g., 2427.0)
        integer_part, decimal_part = str_val.split(".")
        if decimal_part == "0":
            # It's a whole number, return the integer part
            return integer_part
        # It has meaningful decimal places, preserve them
        return str_val

    # If it's already a string without decimal point, return as is
    return str_val


def format_barcode(value: object | None) -> str | None:
    """Format barcodes by removing .0 and preserving all digits."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return None
    # Convert to string and remove .0 if present
    str_val = str(value).rstrip(".0")
    # Return None if barcode is invalid
    if not is_valid_barcode(str_val):
        return None
    return str_val


def create_family_id(df: pd.DataFrame, start_num: int = 1, prefix: str = "FAM") -> pd.DataFrame:
    """Create a Family_ID column with format FAMxxxx (sequentially numbered)."""
    # Get unique family names
    unique_families = df["Family"].dropna().unique()

    # Create a dictionary mapping family names to IDs
    family_id_map = {}
    for i, family in enumerate(unique_families, start=start_num):
        family_id_map[family] = f"{prefix}{i:04d}"

    # Create the Family_ID column
    df["Family_ID"] = df["Family"].map(family_id_map)  # type: ignore[reportArgumentType]

    return df


def _format_columns(
    df: pd.DataFrame, numeric_columns: list[str], id_columns: list[str]
) -> pd.DataFrame:
    """Format barcodes and IDs, dropping rows with invalid barcodes."""
    for col in numeric_columns:
        if col in df.columns:
            df[col] = df[col].apply(format_barcode)
            df = df.dropna(subset=[col])

    for col in id_columns:
        if col in df.columns:
            df[col] = df[col].apply(format_numeric_id)

    return df



def main() -> int:  # noqa: C901,PLR0911,PLR0912,PLR0915
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Merge TNLS files and create/update Family and Family_ID columns")
    parser.add_argument(
        "--ids-file",
        default="data/TNLS/STR_mtDNA_with_ids.tsv",
        help="Path to the STR_mtDNA_with_ids.tsv file",
    )
    parser.add_argument(
        "--family-file",
        required=True,
        help="Path to the file containing family information (e.g., TNLS_with_soldier_names.tsv)",
    )
    parser.add_argument(
        "--output-file",
        default="data/TNLS/STR_mtDNA_with_ids_family.tsv",
        help="Path to the output merged file",
    )
    parser.add_argument("--start-num", type=int, default=1, help="Starting number for family IDs (default: 1)")
    args = parser.parse_args()

    # Check if files exist
    if not Path(args.ids_file).exists():
        logger.error(f"File not found: {args.ids_file}")
        return 1

    if not Path(args.family_file).exists():
        logger.error(f"File not found: {args.family_file}")
        return 1

    try:
        # Read the input files
        logger.info(f"Reading IDs file: {args.ids_file}")
        ids_df = pd.read_csv(args.ids_file, sep="\t")

        logger.info(f"Reading family info file: {args.family_file}")
        family_df = pd.read_csv(args.family_file, sep="\t")

        # Check if this is the TNLS_with_soldier_names.tsv format
        if "Soldier_Name" in family_df.columns and "Family" not in family_df.columns:
            logger.info("Processing TNLS_with_soldier_names.tsv format")
            # Rename Soldier_Name to Family
            family_df = family_df.rename(columns={"Soldier_Name": "Family"})

            # Create Family_ID column if it doesn't exist
            if "Family_ID" not in family_df.columns:
                family_df = create_family_id(family_df, start_num=args.start_num)

        # Format numeric columns in both dataframes
        numeric_columns = ["mtDNA_barcode", "STR_barcode"]
        id_columns = ["mtDNA_ID", "STR_ID"]

        ids_df = _format_columns(ids_df, numeric_columns, id_columns)
        family_df = _format_columns(family_df, numeric_columns, id_columns)

        # Check if required columns are present in the family file
        barcode_columns = ["mtDNA_barcode", "STR_barcode"]

        # At least one of the barcode columns must be present
        if not any(col in family_df.columns for col in barcode_columns):
            logger.error(f"The family file must contain at least one of these columns: {', '.join(barcode_columns)}")
            return 1

        # Check if Family and Family_ID columns are present
        if "Family" not in family_df.columns:
            logger.error("Family column is missing from the family file")
            return 1

        if "Family_ID" not in family_df.columns:
            logger.warning("Family_ID column is missing from the family file. Creating it")
            family_df = create_family_id(family_df, start_num=args.start_num)

        # Create a copy of the IDs dataframe that we'll merge into
        merged_df = ids_df.copy()

        # Add Family and Family_ID columns if they don't exist
        if "Family" not in merged_df.columns:
            merged_df["Family"] = None
        if "Family_ID" not in merged_df.columns:
            merged_df["Family_ID"] = None

        # First try to merge based on mtDNA_barcode
        if "mtDNA_barcode" in family_df.columns and "mtDNA_barcode" in ids_df.columns:
            logger.info("Merging files based on mtDNA_barcode")
            mtdna_family_df = family_df[["mtDNA_barcode", "Family", "Family_ID"]].drop_duplicates()
            merged_mtdna = merged_df.merge(mtdna_family_df, on="mtDNA_barcode", how="left", suffixes=("", "_mtdna"))

            # Update Family and Family_ID where they are not already set
            for col in ["Family", "Family_ID"]:
                mask = merged_df[col].isna() & merged_mtdna[col].notna()
                merged_df.loc[mask, col] = merged_mtdna.loc[mask, col]

        # Next, try to merge based on STR_barcode for records that haven't been matched yet
        if "STR_barcode" in family_df.columns and "STR_barcode" in ids_df.columns:
            logger.info("Merging remaining records based on STR_barcode")
            str_family_df = family_df[["STR_barcode", "Family", "Family_ID"]].drop_duplicates()

            # Create a mask for rows that still don't have family info
            unmatched_mask = merged_df["Family"].isna()

            # Only proceed if there are unmatched records
            if unmatched_mask.any():  # type: ignore[reportGeneralTypeIssues]
                # Get unmatched records
                unmatched_df = merged_df[unmatched_mask].copy()

                # Perform the merge on unmatched records
                str_merged = unmatched_df.merge(str_family_df, on="STR_barcode", how="left", suffixes=("", "_str"))

                # Update the unmatched records with any STR matches
                for col in ["Family", "Family_ID"]:
                    merge_mask = str_merged[col].notna()
                    if merge_mask.any():  # type: ignore[reportGeneralTypeIssues]
                        # Get indices from original dataframe for the matched records
                        matched_indices = unmatched_df.index[str_merged[merge_mask].index]
                        # Update values in the merged dataframe
                        merged_df.loc[matched_indices, col] = str_merged.loc[merge_mask, col].to_numpy()

        # Save the merged file
        logger.info(f"Saving merged data to: {args.output_file}")
        merged_df.to_csv(args.output_file, sep="\t", index=False)

        # Print summary
        match_count = merged_df["Family"].notna().sum()
        total_count = len(merged_df)
        logger.success(
            f"Merge completed successfully! Matched {match_count} out of {total_count} records "
        f"({match_count / total_count * 100:.1f}%)",
        )

    except pd.errors.EmptyDataError:
        logger.error("One of the input files is empty")
        return 1
    except pd.errors.ParserError:
        logger.error("Could not parse one of the TSV files. Please check the file format.")
        return 1
    except (FileNotFoundError, OSError, ValueError, KeyError) as e:
        logger.exception(f"An error occurred: {e!s}")
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main())
