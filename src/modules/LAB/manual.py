#!/usr/bin/env python
import warnings

import pandas as pd
from loguru import logger

# Suppress openpyxl Data Validation warning (one-time only)
warnings.filterwarnings("ignore", message="Data Validation extension is not supported and will be removed")


def main() -> None:
    tracking_file = "src/modules/lab/01. mtDNA Analysis Tracking_CNVV.xlsx"

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

    # Update the "Batch mtDNA" column to have 'mtDNA_' prefix for batch start with numbers
    tracking_df["Batch mtDNA"] = tracking_df["Batch mtDNA"].apply(
        lambda x: f"mtDNA_{x}" if isinstance(x, (int, float)) or (isinstance(x, str) and x[0].isdigit()) else x,
    )

    # Save to TSV file
    tracking_df.to_csv("src/modules/lab/tracking.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
