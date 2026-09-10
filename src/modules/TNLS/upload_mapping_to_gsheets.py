#!/usr/bin/env python
"""
Script to upload STR, mtDNA, and TNLS data to Google Sheets.
This script takes the STR, mtDNA, and TNLS mapping files and uploads them to specified Google Sheets.
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any

import gspread
import numpy as np
import pandas as pd
from loguru import logger
from oauth2client.service_account import ServiceAccountCredentials


# Implement utility functions directly in this file instead of importing
def load_tsv_data(file_path: str) -> pd.DataFrame:
    """
    Load data from a TSV file into a pandas DataFrame.

    Args:
        file_path: Path to the TSV file

    Returns:
        DataFrame containing the TSV data
    """
    return pd.read_csv(file_path, sep="\t", low_memory=False)


def clean_data_for_sheets(values: list[list[Any]]) -> list[list[Any]]:
    """Clean data to ensure all values are JSON-compatible for Google Sheets API."""
    cleaned_values = []

    for row in values:
        cleaned_row = []
        for cell in row:
            # Handle NaN, inf values and convert all non-compatible values to strings
            if pd.isna(cell) or (isinstance(cell, float) and (np.isnan(cell) or np.isinf(cell))):
                cleaned_row.append("")
            elif isinstance(cell, (int, float, bool)):
                cleaned_row.append(cell)
            else:
                cleaned_row.append(str(cell))
        cleaned_values.append(cleaned_row)

    return cleaned_values


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Upload STR, mtDNA, and TNLS data to Google Sheets")
    parser.add_argument("--str_file", required=True, help="Path to the STR mapping file (.tsv)")
    parser.add_argument("--mtdna_file", required=True, help="Path to the mtDNA mapping file (.tsv)")
    parser.add_argument("--tnls_file", help="Path to the TNLS integrated data file (.tsv)")
    parser.add_argument("--tnls_filtered_file", help="Path to the filtered TNLS integrated data file (.tsv)")
    parser.add_argument("--sheet_id", required=True, help="Google Sheet ID to upload data to")
    parser.add_argument("--credentials", required=True, help="Path to Google API credentials JSON file")
    parser.add_argument("--str_worksheet", default="STR", help="Name of the worksheet to upload STR data to")
    parser.add_argument("--mtdna_worksheet", default="mtDNA", help="Name of the worksheet to upload mtDNA data to")
    parser.add_argument("--tnls_worksheet", default="TNLS", help="Name of the worksheet to upload TNLS data to")
    parser.add_argument(
        "--tnls_filtered_worksheet",
        default="TNLS_mtDNA",
        help="Name of the worksheet to upload filtered TNLS data to",
    )
    parser.add_argument("--replace", action="store_true", help="Replace the worksheets if they exist")
    parser.add_argument("--log_file", help="Path to log file (optional)")
    parser.add_argument(
        "--sensitive_sheet_id",
        help="Google Sheet ID to upload complete data (including sensitive data) to",
    )
    parser.add_argument(
        "--sensitive_worksheet",
        help="Name of the worksheet for sensitive data (defaults to same as regular worksheets)",
    )
    return parser.parse_args()


def authenticate_google_sheets(credentials_file: str) -> gspread.Client | None:
    """
    Authenticate with Google Sheets API.

    Args:
        credentials_file: Path to the Google API credentials JSON file

    Returns:
        Authorized gspread client, or None if authentication fails
    """
    # Define the scope for Google Sheets
    scopes = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/drive"]

    try:
        # Setup credentials
        logger.debug(f"Attempting to authenticate using credentials file: {credentials_file}")
        creds = ServiceAccountCredentials.from_json_keyfile_name(credentials_file, scopes)  # type: ignore[reportArgumentType]
        # Authorize the client
        client = gspread.authorize(creds)  # type: ignore[reportArgumentType]
    except (OSError, ValueError, gspread.exceptions.APIError) as e:
        logger.error(f"Authentication failed: {e}")
        logger.exception("Traceback:")
        return None
    else:
        logger.success("Authentication successful!")
        return client


def upload_to_google_sheets(
    client: gspread.Client,
    sheet_id: str,
    worksheet_name: str,
    data_df: pd.DataFrame,
    *,
    replace: bool = False,
) -> bool:
    """
    Upload data to Google Sheets.

    Args:
        client: Authorized gspread client
        sheet_id: Google Sheet ID
        worksheet_name: Name of the worksheet
        data_df: DataFrame with data to upload
        replace: Whether to replace the worksheet if it exists

    Returns:
        True if successful, False otherwise
    """
    try:
        # Open the spreadsheet
        sheet = client.open_by_key(sheet_id)

        # Check if worksheet exists, create if it doesn't
        try:
            worksheet = sheet.worksheet(worksheet_name)
            if replace:
                # Clear the worksheet if replace is True
                worksheet.clear()
            else:
                # If not replacing, create a new worksheet with a timestamp
                timestamp = time.strftime("%Y%m%d_%H%M%S")
                worksheet_name = f"{worksheet_name}_{timestamp}"
                worksheet = sheet.add_worksheet(title=worksheet_name, rows=data_df.shape[0] + 1, cols=data_df.shape[1])
        except gspread.exceptions.WorksheetNotFound:
            # Create worksheet if it doesn't exist
            worksheet = sheet.add_worksheet(title=worksheet_name, rows=data_df.shape[0] + 1, cols=data_df.shape[1])

        # Convert DataFrame to list of lists (including headers)
        raw_values = [*data_df.columns.tolist(), *data_df.to_numpy().tolist()]

        # Clean the data to ensure JSON compatibility
        values = clean_data_for_sheets(raw_values)

        # Update the worksheet
        try:
            worksheet.update(values)
        except gspread.exceptions.APIError as e:
            # If batch update fails, try cell by cell update as fallback
            logger.warning(f"Batch update failed, trying cell by cell update: {e}")
            for r_idx, row in enumerate(values):
                for c_idx, value in enumerate(row):
                    try:
                        worksheet.update_cell(r_idx + 1, c_idx + 1, value)
                    except gspread.exceptions.APIError as cell_err:
                        logger.warning(f"Could not update cell at ({r_idx + 1}, {c_idx + 1}): {cell_err}")
                        # Set to empty string if all else fails
                        worksheet.update_cell(r_idx + 1, c_idx + 1, "")

    except (OSError, ValueError, gspread.exceptions.APIError) as e:
        logger.error(f"Error uploading to Google Sheets: {e}")
        logger.exception("Traceback:")
        return False
    else:
        logger.success(f"Successfully uploaded data to Google Sheets worksheet: {worksheet_name}")
        logger.info(f"Sheet URL: https://docs.google.com/spreadsheets/d/{sheet_id}")
        return True


def count_and_summarize_data(
    str_df: pd.DataFrame | None,
    mtdna_df: pd.DataFrame | None,
    tnls_df: pd.DataFrame | None = None,
    tnls_filtered_df: pd.DataFrame | None = None,
) -> None:
    """
    Count and summarize data from STR, mtDNA, and TNLS mapping files.

    Args:
        str_df: DataFrame with STR data
        mtdna_df: DataFrame with mtDNA data
        tnls_df: DataFrame with TNLS data (optional)
        tnls_filtered_df: DataFrame with filtered TNLS data (optional)
    """
    str_count = len(str_df) if str_df is not None else 0
    mtdna_count = len(mtdna_df) if mtdna_df is not None else 0
    tnls_count = len(tnls_df) if tnls_df is not None else 0
    tnls_filtered_count = len(tnls_filtered_df) if tnls_filtered_df is not None else 0

    logger.info("----- Data Summary -----")
    logger.info(f"STR mappings: {str_count} records")
    logger.info(f"mtDNA mappings: {mtdna_count} records")
    if tnls_df is not None:
        logger.info(f"TNLS integrated data: {tnls_count} records")
    if tnls_filtered_df is not None:
        logger.info(f"TNLS filtered data: {tnls_filtered_count} records")


def remove_sensitive_columns(df: pd.DataFrame, columns_to_exclude: list[str] | None = None) -> pd.DataFrame:
    """
    Remove sensitive columns from the DataFrame before uploading.

    Args:
        df: DataFrame to process
        columns_to_exclude: List of column names to exclude

    Returns:
        DataFrame with sensitive columns removed
    """
    if df is None:
        return None

    if columns_to_exclude is None:
        columns_to_exclude = ["Name", "Gender"]

    removed_columns = [col for col in columns_to_exclude if col in df.columns]

    if removed_columns:
        logger.info(f"Removing sensitive columns: {', '.join(removed_columns)}")
        return df.drop(columns=removed_columns, errors="ignore")

    return df


def verify_sheet_access(client: gspread.Client, sheet_id: str) -> bool:
    """
    Verify that the Google Sheet exists and is accessible.

    Args:
        client: Authorized gspread client
        sheet_id: Google Sheet ID to verify

    Returns:
        True if sheet is accessible, False otherwise
    """
    try:
        logger.debug(f"Verifying access to sheet ID: {sheet_id}")
        sheet = client.open_by_key(sheet_id)
        # Try to get the sheet title to confirm access
        title = sheet.title
        logger.debug(f"Successfully accessed sheet: '{title}' (ID: {sheet_id})")
    except gspread.exceptions.APIError as e:
        if "404" in str(e):
            logger.error(f"Sheet with ID {sheet_id} not found. Please check the ID is correct.")
        elif "403" in str(e):
            logger.error(f"Permission denied for sheet ID {sheet_id}. Please share the sheet with the service account.")
        elif "400" in str(e) and "This operation is not supported for this document" in str(e):
            logger.error(
                f"Sheet with ID {sheet_id} exists but does not support API operations. "
                f"Make sure it's a Google Sheet (not a Doc or other format).",
            )
        else:
            logger.error(f"Error accessing sheet ID {sheet_id}: {e}")
        return False
    except (OSError, ValueError) as e:
        logger.error(f"Unexpected error verifying sheet ID {sheet_id}: {e}")
        logger.exception("Traceback:")
        return False
    else:
        return True
        return False


def upload_complete_data(
    client: gspread.Client,
    sheet_id: str,
    str_data: pd.DataFrame | None,
    mtdna_data: pd.DataFrame | None,
    tnls_data: pd.DataFrame | None,
    tnls_filtered_data: pd.DataFrame | None = None,
    str_worksheet: str = "STR",
    mtdna_worksheet: str = "mtDNA",
    tnls_worksheet: str = "TNLS",
    tnls_filtered_worksheet: str = "TNLS_mtDNA",
    *,
    replace: bool = False,
) -> None:
    """
    Upload complete data including sensitive fields to a separate Google Sheet.

    Args:
        client: Authorized gspread client
        sheet_id: Google Sheet ID to upload complete data to
        str_data: DataFrame with STR data
        mtdna_data: DataFrame with mtDNA data
        tnls_data: DataFrame with TNLS data
        tnls_filtered_data: DataFrame with filtered TNLS data
        str_worksheet: Name of worksheet for STR data
        mtdna_worksheet: Name of worksheet for mtDNA data
        tnls_worksheet: Name of worksheet for TNLS data
        tnls_filtered_worksheet: Name of worksheet for filtered TNLS data
        replace: Whether to replace existing worksheets
    """

    # Upload complete STR data
    if str_data is not None:
        upload_to_google_sheets(client, sheet_id, str_worksheet, str_data, replace=replace)

    # Upload complete mtDNA data
    if mtdna_data is not None:
        upload_to_google_sheets(client, sheet_id, mtdna_worksheet, mtdna_data, replace=replace)

    # Upload complete TNLS data
    if tnls_data is not None:
        upload_to_google_sheets(client, sheet_id, tnls_worksheet, tnls_data, replace=replace)

    # Upload complete filtered TNLS data
    if tnls_filtered_data is not None:
        upload_to_google_sheets(client, sheet_id, tnls_filtered_worksheet, tnls_filtered_data, replace=replace)

    logger.success(f"Successfully uploaded complete data to Google Sheets (ID: {sheet_id})")


def main() -> None:  # noqa: C901,PLR0912,PLR0915
    """Main function."""
    args = parse_args()

    # Add file logger if log_file is provided
    if args.log_file:
        logger.add(sink=args.log_file, rotation="10 MB", retention="1 week", level="DEBUG")

    # Load the data
    try:
        str_data = load_tsv_data(args.str_file) if args.str_file else None
        mtdna_data = load_tsv_data(args.mtdna_file) if args.mtdna_file else None
        tnls_data = load_tsv_data(args.tnls_file) if args.tnls_file else None
        tnls_filtered_data = load_tsv_data(args.tnls_filtered_file) if args.tnls_filtered_file else None
        logger.success("Input files loaded successfully")
    except (FileNotFoundError, pd.errors.ParserError, OSError) as e:
        logger.error(f"Failed to load one or more input files: {e}")
        logger.exception("Load error details:")
        return

    # Provide summary of data
    count_and_summarize_data(str_data, mtdna_data, tnls_data, tnls_filtered_data)

    # Authenticate with Google Sheets
    client = authenticate_google_sheets(args.credentials)
    if not client:
        return

    # Verify main sheet is accessible
    if not verify_sheet_access(client, args.sheet_id):
        logger.error(f"Unable to access main sheet (ID: {args.sheet_id}). Check permissions and try again.")
        return

    # Upload complete data (including sensitive data) to separate sheet if specified
    if args.sensitive_sheet_id:
        # Verify sensitive sheet is accessible
        if not verify_sheet_access(client, args.sensitive_sheet_id):
            logger.warning(
                f"Unable to access sensitive data sheet (ID: {args.sensitive_sheet_id}). "
                "Skipping upload of complete data. Will only upload to main sheet.",
            )

            # Try to extract service account email to help with troubleshooting
            try:
                with Path(args.credentials).open() as f:
                    creds_data = json.load(f)
                    client_email = creds_data.get("client_email")
                    if client_email:
                        logger.info(f"To fix this, share the Google Sheet with: {client_email}")
                        logger.info("Make sure the service account has Editor permissions")
            except (OSError, json.JSONDecodeError, KeyError):
                logger.debug("Could not extract service account email from credentials")
        else:
            # Use sensitive_worksheet if provided, otherwise use the default worksheet names
            str_ws = args.sensitive_worksheet or args.str_worksheet
            mtdna_ws = args.sensitive_worksheet or args.mtdna_worksheet
            tnls_ws = args.sensitive_worksheet or args.tnls_worksheet
            tnls_filtered_ws = args.sensitive_worksheet or args.tnls_filtered_worksheet

            upload_complete_data(
                client,
                args.sensitive_sheet_id,
                str_data,
                mtdna_data,
                tnls_data,
                tnls_filtered_data,
                str_ws,
                mtdna_ws,
                tnls_ws,
                tnls_filtered_ws,
                replace=args.replace,
            )

    # Remove sensitive columns from all datasets before uploading to main sheet
    str_data_clean = remove_sensitive_columns(str_data)  # type: ignore[reportArgumentType]
    mtdna_data_clean = remove_sensitive_columns(mtdna_data)  # type: ignore[reportArgumentType]
    tnls_data_clean = remove_sensitive_columns(tnls_data)  # type: ignore[reportArgumentType]
    tnls_filtered_data_clean = remove_sensitive_columns(tnls_filtered_data)  # type: ignore[reportArgumentType]

    # Upload cleaned STR data
    if str_data_clean is not None:
        upload_to_google_sheets(client, args.sheet_id, args.str_worksheet, str_data_clean, replace=args.replace)

    # Upload cleaned mtDNA data
    if mtdna_data_clean is not None:
        upload_to_google_sheets(client, args.sheet_id, args.mtdna_worksheet, mtdna_data_clean, replace=args.replace)

    # Upload cleaned TNLS data
    if tnls_data_clean is not None:
        upload_to_google_sheets(client, args.sheet_id, args.tnls_worksheet, tnls_data_clean, replace=args.replace)

        sheet_id_temp = "1rbCoUx9UhNoK7EnEH-ClBheucevSyHMljV2WSKG2yvQ"
        upload_to_google_sheets(client, sheet_id_temp, args.tnls_worksheet, tnls_data_clean, replace=args.replace)

    # Upload cleaned filtered TNLS data
    if tnls_filtered_data_clean is not None:
        upload_to_google_sheets(
            client,
            args.sheet_id,
            args.tnls_filtered_worksheet,
            tnls_filtered_data_clean,
            replace=args.replace,
        )

    logger.success("Upload process completed.")


if __name__ == "__main__":
    main()
