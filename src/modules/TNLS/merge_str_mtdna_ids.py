#!/usr/bin/env python3
"""
Script to merge STR and mtDNA IDs into the STR_mtDNA.tsv file.
This script reads directly from merged metadata files instead of processing individual batch files.
"""

import argparse

import pandas as pd
from loguru import logger

_DUPLICATE_THRESHOLD = 1


def format_numeric_id(value: object | None) -> object | None:
    """Format numeric IDs by preserving both leading and trailing zeros."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return value

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


def format_barcode(value: object | None) -> object | None:
    """Format barcodes by removing .0 and preserving all digits with proper length."""
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return value
    # Convert to string and remove .0 if present
    str_val = str(value).rstrip(".0")
    # Try to convert to float to check if it's numeric
    try:
        float(str_val)
        # Ensure 12 digits with leading zeros
        return str_val.zfill(12)
    except ValueError:
        return value


def check_duplicates(df: pd.DataFrame, column_name: str) -> bool:
    """Check for duplicate values in a dataframe column and log them."""
    if column_name not in df.columns:
        return False

    duplicates = df[df.duplicated(subset=[column_name], keep=False)]
    if len(duplicates) > 0:
        logger.warning(f"Found {len(duplicates)} rows with duplicate values in column '{column_name}'")
        dup_values = duplicates[column_name].value_counts()  # type: ignore[reportAttributeAccessIssue]
        for value, count in dup_values.items():
            if count > _DUPLICATE_THRESHOLD:  # Only show values that appear more than once
                logger.warning(f"Value '{value}' appears {count} times")
        return True
    return False


def create_mapping_dict(df: pd.DataFrame, key_col: str, value_col: str) -> dict[object, object]:
    """Create a mapping dictionary that handles duplicates by using first occurrence."""
    if key_col not in df.columns or value_col not in df.columns:
        return {}

    # Check if there are duplicates in the key column
    has_duplicates = check_duplicates(df, key_col)

    if has_duplicates:
        logger.info(f"Creating mapping using first occurrence of each {key_col}")
        # Keep only the first occurrence of each key
        unique_df = df.drop_duplicates(subset=[key_col], keep="first")
        # Create mapping dictionary
        mapping = dict(zip(unique_df[key_col], unique_df[value_col].apply(format_numeric_id), strict=False))
    else:
        # If no duplicates, proceed as before
        mapping = dict(zip(df[key_col], df[value_col].apply(format_numeric_id), strict=False))

    return mapping


def main() -> None:  # noqa: C901,PLR0912,PLR0915
    parser = argparse.ArgumentParser(description="Merge STR and mtDNA IDs into STR_mtDNA.tsv file")
    parser.add_argument("--str_merged_file", help="Merged STR metadata file")
    parser.add_argument("--mtdna_merged_file", help="Merged mtDNA metadata file")
    parser.add_argument("--input_file", help="Input STR_mtDNA.tsv file")
    parser.add_argument("--output_file", help="Output TSV file with added IDs")

    args = parser.parse_args()

    # Load the STR_mtDNA.tsv file
    try:
        match_data = pd.read_csv(args.input_file, sep="\t", low_memory=False)

        # Format the barcode columns properly
        match_data["STR_barcode"] = match_data["STR_barcode"].apply(format_barcode)
        match_data["mtDNA_barcode"] = match_data["mtDNA_barcode"].apply(format_barcode)

        logger.info(f"Loaded {len(match_data)} records from {args.input_file}")
    except (FileNotFoundError, pd.errors.ParserError, OSError) as e:
        logger.error(f"Error loading STR_mtDNA.tsv file: {e!s}")
        logger.exception("Traceback:")
        return

    # Load the merged STR metadata file
    try:
        str_data = pd.read_csv(args.str_merged_file, sep="\t", low_memory=False)
        # Identify the barcode and ID columns (assuming they exist in the data)
        str_barcode_col = next((col for col in str_data.columns if "barcode" in col.lower()), None)
        str_id_col = next((col for col in str_data.columns if "lid" in col.lower()), str_data.columns[0])
        str_batch_col = next((col for col in str_data.columns if "batch" in col.lower()), None)

        if str_barcode_col is None:
            logger.warning(
                "Could not identify barcode column in STR metadata file."
            " Using first column as ID and second as barcode.",
            )
            str_id_col = str_data.columns[0]
            str_barcode_col = str_data.columns[1] if len(str_data.columns) > 1 else None

        if str_barcode_col:  # type: ignore[reportGeneralTypeIssues]
            # Format barcodes
            str_data[str_barcode_col]  # type: ignore[reportGeneralTypeIssues,reportGeneralTypeIssues] = str_data[str_barcode_col].apply(format_barcode)

            # Create barcode-to-ID mapping using the new function
            str_barcode_to_id = create_mapping_dict(str_data, str_barcode_col, str_id_col)  # type: ignore[reportArgumentType]

            # Create barcode-to-batch mapping if batch column exists
            str_barcode_to_batch = {}
            if str_batch_col:
                str_barcode_to_batch = create_mapping_dict(str_data, str_barcode_col, str_batch_col)  # type: ignore[reportArgumentType]

            logger.info(f"Loaded {len(str_barcode_to_id)} STR barcode-ID mappings from {args.str_merged_file}")
        else:
            logger.error("Could not find barcode column in STR metadata file.")
            str_barcode_to_id = {}
            str_barcode_to_batch = {}
    except (FileNotFoundError, pd.errors.ParserError, OSError) as e:
        logger.error(f"Error loading STR merged file: {e!s}")
        logger.exception("Traceback:")
        str_barcode_to_id = {}
        str_barcode_to_batch = {}

    # Load the merged mtDNA metadata file
    try:
        mtdna_data = pd.read_csv(args.mtdna_merged_file, sep="\t", low_memory=False)
        # Identify the barcode and ID columns (assuming they exist in the data)
        mtdna_barcode_col = next((col for col in mtdna_data.columns if "barcode" in col.lower()), None)
        mtdna_id_col = next(
            (col for col in mtdna_data.columns if "id" in col.lower() or "sample" in col.lower()),
            mtdna_data.columns[0],
        )
        mtdna_batch_col = next((col for col in mtdna_data.columns if "batch" in col.lower()), None)

        if mtdna_barcode_col is None:
            logger.warning(
                "Could not identify barcode column in mtDNA metadata file."
            " Using first column as ID and second as barcode.",
            )
            mtdna_id_col = mtdna_data.columns[0]
            mtdna_barcode_col = mtdna_data.columns[1] if len(mtdna_data.columns) > 1 else None

        if mtdna_barcode_col:  # type: ignore[reportGeneralTypeIssues]
            # Format barcodes
            mtdna_data[mtdna_barcode_col]  # type: ignore[reportGeneralTypeIssues,reportGeneralTypeIssues] = mtdna_data[mtdna_barcode_col].apply(format_barcode)

            # Create barcode-to-ID mapping using the new function
            mtdna_barcode_to_id = create_mapping_dict(mtdna_data, mtdna_barcode_col, mtdna_id_col)  # type: ignore[reportArgumentType]

            # Create barcode-to-batch mapping if batch column exists
            mtdna_barcode_to_batch = {}
            if mtdna_batch_col:
                mtdna_barcode_to_batch = create_mapping_dict(mtdna_data, mtdna_barcode_col, mtdna_batch_col)  # type: ignore[reportArgumentType]

            logger.info(f"Loaded {len(mtdna_barcode_to_id)} mtDNA barcode-ID mappings from {args.mtdna_merged_file}")
        else:
            logger.error("Could not find barcode column in mtDNA metadata file.")
            mtdna_barcode_to_id = {}
            mtdna_barcode_to_batch = {}
    except (FileNotFoundError, pd.errors.ParserError, OSError) as e:
        logger.error(f"Error loading mtDNA merged file: {e!s}")
        logger.exception("Traceback:")
        mtdna_barcode_to_id = {}
        mtdna_barcode_to_batch = {}

    # Add new columns to match_data
    match_data["STR_ID"] = match_data["STR_barcode"].map(str_barcode_to_id)  # type: ignore[reportArgumentType]
    match_data["mtDNA_ID"] = match_data["mtDNA_barcode"].map(mtdna_barcode_to_id)  # type: ignore[reportArgumentType]

    # Add batch information if available
    if str_barcode_to_batch:
        match_data["STR_batch"] = match_data["STR_barcode"].map(str_barcode_to_batch)  # type: ignore[reportArgumentType]
    if mtdna_barcode_to_batch:
        match_data["mtDNA_batch"] = match_data["mtDNA_barcode"].map(mtdna_barcode_to_batch)  # type: ignore[reportArgumentType]

    # Format the ID columns
    match_data["STR_ID"] = match_data["STR_ID"].apply(format_numeric_id)
    match_data["mtDNA_ID"] = match_data["mtDNA_ID"].apply(format_numeric_id)

    # Count matches found
    str_matches = match_data["STR_ID"].notna().sum()
    mtdna_matches = match_data["mtDNA_ID"].notna().sum()

    logger.info(f"Found STR IDs for {str_matches} out of {len(match_data)} records")
    logger.info(f"Found mtDNA IDs for {mtdna_matches} out of {len(match_data)} records")

    # Count batch matches if columns exist
    if "STR_batch" in match_data.columns:
        str_batch_matches = match_data["STR_batch"].notna().sum()
        logger.info(f"Added STR batch info for {str_batch_matches} out of {len(match_data)} records")

    if "mtDNA_batch" in match_data.columns:
        mtdna_batch_matches = match_data["mtDNA_batch"].notna().sum()
        logger.info(f"Added mtDNA batch info for {mtdna_batch_matches} out of {len(match_data)} records")

    # Save the result to TSV
    match_data.to_csv(args.output_file, sep="\t", index=False)
    logger.success(f"Saved merged data to TSV file: {args.output_file}")


if __name__ == "__main__":
    main()
