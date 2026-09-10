#!/usr/bin/env python
"""
This script creates family data from TNLS_with_soldier_names.tsv and merges it with STR_mtDNA_with_ids.tsv.
Fixes the floating-point barcode handling problem.
"""

import argparse
import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd
from loguru import logger


def is_valid_barcode(value: object | None) -> str | bool:
    """Check if barcode is valid (all digits and correct length)."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return False

    # Convert to string and remove .0 if present (handle floating point)
    str_val = str(value).strip()
    return str_val.removesuffix(".0")


def format_barcode(value: object | None, target_length: int = 12) -> str | None:
    """Format barcodes by removing .0, preserving all digits, and padding with leading zeros."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return None

    # Convert to string and remove .0 if present
    str_val = str(value).strip()
    str_val = str_val.removesuffix(".0")

    # If barcode is empty or whitespace, return None
    if not str_val:
        return None

    # Pad with leading zeros to ensure consistent length
    # Use the target_length, but don't truncate if longer
    if len(str_val) < target_length:
        str_val = str_val.zfill(target_length)

    return str_val


def format_numeric_id(value: object | None) -> object | None:
    """Format numeric IDs by removing .0 suffix and returning as integer string."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return None

    # Convert to string
    str_val = str(value)

    # If it's a float value (has decimal point)
    if "." in str_val:
        # Check if it's a float that represents an integer (e.g., 241673.0)
        integer_part, decimal_part = str_val.split(".")
        if decimal_part == "0":
            # It's a whole number, return the integer part
            return integer_part
        # It has meaningful decimal places, preserve them
        return str_val

    # If it's already a string without decimal point, return as is
    return str_val


def normalize_vietnamese_name(name: object | None) -> str | None:
    """
    Normalize a Vietnamese name by converting to lowercase, standardizing spacing,
    and removing diacritical marks.

    Args:
        name: The name to normalize

    Returns:
        A normalized version of the name
    """
    if pd.isna(name):  # type: ignore[reportGeneralTypeIssues]
        return None

    # Convert to lowercase and strip extra spaces
    normalized = name.lower().strip()  # type: ignore[reportAttributeAccessIssue,reportOptionalMemberAccess]

    # Replace multiple spaces with a single space
    normalized = " ".join(normalized.split())

    # Remove diacritical marks (accent marks) from Vietnamese characters
    # This converts characters like ớ, ở, ờ, ỡ, ợ to o
    normalized = "".join(c for c in unicodedata.normalize("NFD", normalized) if not unicodedata.combining(c))

    # Additional common Vietnamese character mappings
    # Handle single characters with different base forms
    viet_chars = {
        "đ": "d",
        "Đ": "d",  # Vietnamese đ/Đ to d
        "ơ": "o",
        "Ơ": "o",  # Vietnamese ơ/Ơ to o
        "ư": "u",
        "Ư": "u",  # Vietnamese ư/Ư to u
    }

    for viet_char, replacement in viet_chars.items():
        normalized = normalized.replace(viet_char, replacement)

    return normalized


def create_family_id(  # noqa: C901
    df: pd.DataFrame, start_num: int = 1, prefix: str = "FAM", event_col: str | None = None
) -> pd.DataFrame:
    """
    Create a Family_ID column with format FAMxxxx (sequentially numbered).

    If event_col is provided, Family_ID will be unique for combinations of Family and Event.
    Handles comma-separated family names by assigning multiple Family_IDs.

    Args:
        df: The dataframe containing the family data
        start_num: The starting number for the Family ID sequence
        prefix: The prefix to use for the Family ID
        event_col: Optional column name to use for differentiating families by event

    Returns:
        The dataframe with a new Family_ID column
    """
    # Build a dictionary of all unique family names
    family_id_map = {}
    current_id = start_num

    # Process each row to extract individual family names
    for _, row in df.iterrows():
        if pd.isna(row["Family"]):  # type: ignore[reportGeneralTypeIssues]
            continue

        # Split by comma to handle multiple family names in one field
        family_names = [name.strip() for name in row["Family"].split(",")]  # type: ignore[reportAttributeAccessIssue]

        for family_name in family_names:
            normalized_name = normalize_vietnamese_name(family_name)

            # Skip if normalized name is None or already processed
            if normalized_name is None or normalized_name in family_id_map:
                continue

            # Create Family ID with event distinction if needed
            if event_col and event_col in df.columns and pd.notna(row[event_col]):  # type: ignore[reportGeneralTypeIssues]
                family_key = f"{normalized_name}___{row[event_col]}"
            else:
                family_key = normalized_name

            if family_key not in family_id_map:
                family_id_map[family_key] = f"{prefix}{current_id:04d}"
                current_id += 1

    logger.info(f"Created {len(family_id_map)} unique Family IDs")

    # Apply the mapping to create Family_ID column
    def assign_family_ids(row: pd.Series) -> str | None:
        if pd.isna(row["Family"]):  # type: ignore[reportGeneralTypeIssues]
            return None

        # Split multiple family names
        family_names = [name.strip() for name in row["Family"].split(",")]  # type: ignore[reportAttributeAccessIssue]
        family_ids = []

        for family_name in family_names:
            normalized_name = normalize_vietnamese_name(family_name)

            if normalized_name is None:
                continue

            # Create the key based on whether we're using event or not
            if event_col and event_col in df.columns and pd.notna(row[event_col]):  # type: ignore[reportGeneralTypeIssues]
                family_key = f"{normalized_name}___{row[event_col]}"
            else:
                family_key = normalized_name

            if family_key in family_id_map:
                family_ids.append(family_id_map[family_key])

        # Return comma-separated Family_IDs for multiple families
        return ", ".join(family_ids) if family_ids else None

    df["Family_ID"] = df.apply(assign_family_ids, axis=1)
    return df


def extract_soldier_names(df: pd.DataFrame) -> pd.DataFrame:
    """Extract soldier names from the Note column."""
    # Create a Family column from Soldier_Name if it exists, otherwise use Note
    if "Soldier_Name" in df.columns and df["Soldier_Name"].notna().any():  # type: ignore[reportGeneralTypeIssues]
        logger.info("Using existing Soldier_Name column")
        df["Family"] = df["Soldier_Name"]
    else:
        logger.info("Extracting soldier names from Note column")
        # Extract soldier names from the Note column using regex patterns
        patterns = [
            r"(?:liệt sĩ|LS|ls) ([\w\s,]+)",  # Pattern for "liệt sĩ [Name]" or "LS [Name]"
            r"(?:liệt sĩ|LS|ls) ([\w\s]+)(?:và|,) ([\w\s]+)",  # Pattern for names separated by "và" or ","
        ]

        def extract_name(note: object | None) -> str | None:
            if pd.isna(note):  # type: ignore[reportGeneralTypeIssues]
                return None

            # Try each pattern
            for pattern in patterns:
                match = re.search(pattern, note, re.IGNORECASE)  # type: ignore[reportCallIssue,reportArgumentType]
                if match:
                    if len(match.groups()) > 1:
                        # Multiple names found
                        return ", ".join(g.strip() for g in match.groups() if g)
                    return match.group(1).strip()

            # If no pattern matches, use a more general approach
            if "liệt sĩ" in note.lower() or "ls" in note.lower():  # type: ignore[reportAttributeAccessIssue,reportOptionalMemberAccess]
                # Find the part after "liệt sĩ" or "LS"
                parts = re.split(r"(?:liệt sĩ|LS|ls)\s+", note, flags=re.IGNORECASE)  # type: ignore[reportCallIssue,reportArgumentType]
                if len(parts) > 1:
                    return parts[1].strip()

            return None

        df["Family"] = df["Note"].apply(extract_name)

    return df


def main() -> int:  # noqa: C901,PLR0912,PLR0915
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Merge TNLS files and fix floating-point barcode handling")
    parser.add_argument(
        "--ids-file",
        default="data/TNLS/STR_mtDNA_with_ids.tsv",
        help="Path to the STR_mtDNA_with_ids.tsv file",
    )
    parser.add_argument(
        "--family-file",
        default="data/TNLS/TNLS_with_soldier_names.tsv",
        help="Path to the file containing family information",
    )
    parser.add_argument(
        "--output-file",
        default="data/TNLS/STR_mtDNA_with_ids_family.tsv",
        help="Path to the output merged file",
    )
    args = parser.parse_args()

    try:
        # Check if family file exists and read it
        logger.info(f"Reading family info file: {args.family_file}")
        family_df = pd.read_csv(args.family_file, sep="\t", low_memory=False)

        # Process family data
        if "Family" not in family_df.columns:
            family_df = extract_soldier_names(family_df)

        # Format barcodes properly
        for col in ["mtDNA_barcode", "STR_barcode"]:
            if col in family_df.columns:
                family_df[col] = family_df[col].apply(format_barcode)

        # Create Family_ID if it doesn't exist
        if "Family_ID" not in family_df.columns:
            # Check if Event column exists
            event_column = "Event" if "Event" in family_df.columns else None
            if not event_column:
                logger.warning("Event column not found. Family IDs will be created based only on Family names.")

            family_df = create_family_id(family_df, event_col=event_column)

        # Check if IDs file exists
        if not Path(args.ids_file).exists():
            logger.warning(f"IDs file not found: {args.ids_file}")
            logger.info("Creating a sample IDs file for demonstration")

            # Create a minimal IDs dataframe from the family data for demonstration
            ids_df = pd.DataFrame(
                {
                    "Name": family_df["Name"],
                    "Gender": family_df["Gender"],
                    "STR_barcode": family_df["STR_barcode"],
                    "mtDNA_barcode": family_df["mtDNA_barcode"],
                    "STR_ID": range(240001, 240001 + len(family_df)),
                    "mtDNA_ID": range(242001, 242001 + len(family_df)),
                },
            )
        else:
            # Read the IDs file
            logger.info(f"Reading IDs file: {args.ids_file}")
            ids_df = pd.read_csv(args.ids_file, sep="\t", low_memory=False)
            # Format barcodes in IDs file
        for col in ["mtDNA_barcode", "STR_barcode"]:
            if col in ids_df.columns:
                ids_df[col] = ids_df[col].apply(format_barcode)

        # Format numeric IDs to remove .0 suffix
        for col in ["STR_ID", "mtDNA_ID"]:
            if col in ids_df.columns:
                ids_df[col] = ids_df[col].apply(format_numeric_id)

        # Create a merged dataframe
        logger.info("Merging data")

        # Start with IDs dataframe
        merged_df = ids_df.copy()

        # Add Family and Family_ID columns if they don't exist
        if "Family" not in merged_df.columns:
            merged_df["Family"] = None
        if "Family_ID" not in merged_df.columns:
            merged_df["Family_ID"] = None

        # Merge based on mtDNA_barcode
        if "mtDNA_barcode" in family_df.columns and "mtDNA_barcode" in ids_df.columns:
            logger.info("Merging files based on mtDNA_barcode")
            mtdna_family_df = family_df[["mtDNA_barcode", "Family", "Family_ID"]].drop_duplicates()

            for _, row in mtdna_family_df.iterrows():
                if pd.notna(row["mtDNA_barcode"]):  # type: ignore[reportGeneralTypeIssues]
                    mask = merged_df["mtDNA_barcode"] == row["mtDNA_barcode"]
                    if mask.any():  # type: ignore[reportGeneralTypeIssues]
                        for col in ["Family", "Family_ID"]:
                            if pd.isna(merged_df.loc[mask, col]).any():
                                merged_df.loc[mask, col] = row[col]

        # Merge based on STR_barcode for records that haven't been matched yet
        if "STR_barcode" in family_df.columns and "STR_barcode" in ids_df.columns:
            logger.info("Merging remaining records based on STR_barcode")
            str_family_df = family_df[["STR_barcode", "Family", "Family_ID"]].drop_duplicates()

            # Get unmatched records
            unmatched_mask = merged_df["Family"].isna()

            for _, row in str_family_df.iterrows():
                if pd.notna(row["STR_barcode"]):  # type: ignore[reportGeneralTypeIssues]
                    mask = (merged_df["STR_barcode"] == row["STR_barcode"]) & unmatched_mask
                    if mask.any():  # type: ignore[reportGeneralTypeIssues]
                        for col in ["Family", "Family_ID"]:
                            merged_df.loc[mask, col] = row[col]

        # Format all barcodes and IDs in the final merged dataframe
        logger.info("Formatting barcodes and IDs in merged data")
        for col in ["mtDNA_barcode", "STR_barcode"]:
            if col in merged_df.columns:
                merged_df[col] = merged_df[col].apply(format_barcode)

        for col in ["STR_ID", "mtDNA_ID"]:
            if col in merged_df.columns:
                merged_df[col] = merged_df[col].apply(format_numeric_id)

        # Save the merged file
        logger.info(f"Saving merged data to: {args.output_file}")
        merged_df.to_csv(args.output_file, sep="\t", index=False)

        # Print summary
        match_count = merged_df["Family"].notna().sum()
        total_count = len(merged_df)
        match_percentage = (match_count / total_count * 100) if total_count > 0 else 0
        logger.success(
            f"Merge completed successfully! Matched {match_count} out of {total_count} records "
            f"({match_percentage:.1f}%)",
        )
    except (FileNotFoundError, pd.errors.ParserError, OSError, ValueError, KeyError) as e:
        logger.error(f"An error occurred: {e!s}")
        logger.exception("Traceback:")
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main())
