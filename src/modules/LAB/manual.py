#!/usr/bin/env python
import warnings
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from loguru import logger

# Suppress openpyxl Data Validation warning (one-time only)
warnings.filterwarnings("ignore", message="Data Validation extension is not supported and will be removed")


DEFAULT_TRACKING_FILE = Path("data/modules/LAB/tracking.xlsx")
DEFAULT_OUTPUT_FILE = Path("data/modules/LAB/tracking.tsv")


def normalize_batch_name(value: object) -> object:
    """Prefix only standalone numeric batch identifiers with ``mtDNA_``."""
    if isinstance(value, (int, float)) or (isinstance(value, str) and value.isdigit()):
        return f"mtDNA_{value}"
    return value


def arg_parser() -> ArgumentParser:
    """Create the tracking-workbook conversion argument parser."""
    parser = ArgumentParser(description="Convert the LAB tracking workbook to TSV")
    parser.add_argument(
        "--tracking-excel",
        default=DEFAULT_TRACKING_FILE,
        type=Path,
        help=f"Tracking workbook path (default: {DEFAULT_TRACKING_FILE})",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_FILE,
        type=Path,
        help=f"Tracking TSV output path (default: {DEFAULT_OUTPUT_FILE})",
    )
    return parser


def main() -> None:
    """Convert the LAB tracking workbook's active batch sheet to TSV."""
    args = arg_parser().parse_args()
    tracking_file: Path = args.tracking_excel
    output_file: Path = args.output

    # Read the 'From Batch 20' sheet into a DataFrame and first row as header, and skip the first two data rows
    tracking_df = pd.read_excel(tracking_file, sheet_name="From Batch 20", header=0, skiprows=[1, 2])

    # Format the "Sequence Released Date" to YYYY-MM-DD
    try:
        tracking_df["Sequence Released Date"] = pd.to_datetime(
            tracking_df["Sequence Released Date"],
            errors="coerce",
        ).dt.strftime("%Y-%m-%d")
    except (ValueError, KeyError, TypeError, pd.errors.OutOfBoundsDatetime) as e:
        logger.error(f"Error formatting 'Sequence Released Date': {e}")
        logger.error("Data that caused the error:")
        logger.error(tracking_df["Sequence Released Date"])

    # Format the "Run1 Analysis Results Released Date" to YYYY-MM-DD
    try:
        tracking_df["Run1 Analysis Results Released Date"] = pd.to_datetime(
            tracking_df["Run1 Analysis Results Released Date"],
            errors="coerce",
        ).dt.strftime("%Y-%m-%d")
    except (ValueError, KeyError, TypeError, pd.errors.OutOfBoundsDatetime) as e:
        logger.error(f"Error formatting 'Run1 Analysis Results Released Date': {e}")
        logger.error("Data that caused the error:")
        logger.error(tracking_df["Run1 Analysis Results Released Date"])

    # Prefix standalone numeric batch identifiers without changing date-based or MS batch IDs.
    tracking_df["Batch mtDNA"] = tracking_df["Batch mtDNA"].apply(normalize_batch_name)

    # Save to TSV file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    tracking_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Tracking TSV saved to {output_file}")


if __name__ == "__main__":
    main()
