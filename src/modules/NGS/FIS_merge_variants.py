#!/usr/bin/env python3
# noqa: N999
"""
Script to merge variant information from all sample Excel files into a single TSV file.
Each sample has an Excel file with variant data that needs to be consolidated.
"""

import argparse
from pathlib import Path

import pandas as pd
from loguru import logger


def extract_sample_id_from_path(excel_path: str) -> str:
    """
    Extract sample ID from the Excel file path.

    Args:
        excel_path: Path to the Excel file

    Returns:
        Sample ID extracted from the path
    """
    # Extract sample ID from path like: .../2_ng-Rep3_L01_88-8/2-MT/2_ng-Rep3_L01_88-8-HyperVar-MT.xlsx
    path_parts = Path(excel_path).parts

    # Find the sample directory (contains the sample ID)
    for part in reversed(path_parts):
        if part.endswith("-MT"):
            continue
        if part.endswith(".xlsx"):
            continue
        if any(char.isdigit() for char in part) and ("_ng-" in part or "control" in part or "NC_" in part):
            return part

    # Fallback: extract from filename
    filename = Path(excel_path).name
    return filename.replace("-HyperVar-MT.xlsx", "")


def find_all_excel_files(excel_path: str) -> list[str]:
    """
    Find all Excel files containing variant data.

    Args:
        excel_path: Path pattern for Excel files

    Returns:
        List of paths to Excel files
    """
    p = Path(excel_path)
    return sorted(str(f) for f in p.parent.glob(p.name))


def read_variant_data(excel_path: str) -> pd.DataFrame:
    """
    Read variant data from an Excel file.

    Args:
        excel_path: Path to the Excel file

    Returns:
        DataFrame containing variant data
    """
    try:
        # Read the Excel file
        df = pd.read_excel(excel_path)

        # Ensure required columns are present
        required_columns = [
            "Position",
            "Marker",
            "Genotype",
            "AlleleFrequency",
            "Ref(Depth):Alt(Depth)",
            "Genotype(raw)",
            "QC_Info",
        ]

        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.warning(f"Missing columns in {excel_path}: {missing_columns}")
            # Add missing columns with empty values
            for col in missing_columns:
                df[col] = ""

        return df[required_columns]  # type: ignore[reportReturnType]

    except (OSError, ValueError, KeyError, pd.errors.ParserError) as e:
        logger.error(f"Error reading {excel_path}: {e!s}")
        return pd.DataFrame()


def merge_all_variants(output_file: str, excel_path: str) -> None:
    """
    Merge variant information from all samples into a single TSV file.

    Args:
        output_file: Output TSV file path
        excel_path: Path pattern for Excel files
    """
    logger.info("Finding Excel files...")
    excel_files = find_all_excel_files(excel_path)

    if not excel_files:
        logger.warning("No Excel files found!")
        return

    logger.info(f"Found {len(excel_files)} Excel files")

    all_data: list[pd.DataFrame] = []

    for file_path in excel_files:
        logger.info(f"Processing: {file_path}")

        # Extract sample ID
        sample_id = extract_sample_id_from_path(file_path)
        logger.debug(f"  Sample ID: {sample_id}")

        # Read variant data
        variant_df = read_variant_data(file_path)

        if variant_df.empty:
            logger.warning(f"  Warning: No data found in {file_path}")
            continue

        # Add sample ID column
        variant_df["Sample_ID"] = sample_id

        # Reorder columns to put Sample_ID first
        columns = ["Sample_ID"] + [col for col in variant_df.columns if col != "Sample_ID"]
        variant_df = variant_df[columns]

        all_data.append(variant_df)  # type: ignore[reportArgumentType]
        logger.debug(f"  Added {len(variant_df)} variants")

    if not all_data:
        logger.error("No data to merge!")
        return

    # Combine all data
    logger.info("Merging all data...")
    merged_df = pd.concat(all_data, ignore_index=True)

    # Sort by Sample_ID and Position for better organization
    merged_df = merged_df.sort_values(["Sample_ID", "Position"])

    # Save to TSV file
    logger.info(f"Saving to {output_file}...")
    merged_df.to_csv(output_file, sep="\t", index=False)

    logger.success("Merge complete!")
    logger.info(f"Total samples: {merged_df['Sample_ID'].nunique()}")
    logger.info(f"Total variants: {len(merged_df)}")
    logger.info(f"Output file: {output_file}")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Merge variant information from all sample Excel files into a single TSV file.",
    )

    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file path")

    parser.add_argument(
        "-e",
        "--excel-path",
        type=str,
        required=True,
        help="Path pattern for Excel files (e.g., 'data/20250820/*MT.xlsx')",
    )

    return parser.parse_args()


def transform_variants(input_file: str, output_file: str) -> None:
    """
    Transform merged variants file to have one row per sample with all genotypes in a single column.

    Args:
        input_file: Path to input TSV file
        output_file: Path to output transformed TSV file
    """
    logger.info(f"Reading merged variants from {input_file}...")
    df = pd.read_csv(input_file, sep="\t")

    logger.info(f"Original data shape: {df.shape}")
    logger.info(f"Number of unique samples: {df['Sample_ID'].nunique()}")

    # Group by Sample_ID and collect all genotypes
    logger.info("Transforming data...")
    transformed_data: list[dict[str, str]] = []

    for sample_id in df["Sample_ID"].unique():
        sample_data = df[df["Sample_ID"] == sample_id]

        # Get all genotypes for this sample, sorted by position
        # Filter out NaN values and convert to strings
        genotypes = sample_data.sort_values("Position")["Genotype"].dropna().astype(str).tolist()  # type: ignore[reportCallIssue]

        # Split comma-separated variants into individual variants
        expanded_genotypes: list[str] = []
        for genotype in genotypes:
            if "," in genotype:
                # Split by comma and add each variant separately
                expanded_genotypes.extend(genotype.split(","))
            else:
                expanded_genotypes.append(genotype)

        # Join genotypes with spaces
        genotypes_str = " ".join(expanded_genotypes)

        transformed_data.append({"Sample_ID": sample_id, "Genotype": genotypes_str})

    # Create new dataframe
    result_df = pd.DataFrame(transformed_data)

    logger.info(f"Transformed data shape: {result_df.shape}")

    # Save to file
    logger.info(f"Saving transformed data to {output_file}...")
    result_df.to_csv(output_file, sep="\t", index=False)

    logger.success("Transformation completed successfully!")


def main() -> None:
    """Main function to run the merging and transformation process."""
    args = parse_arguments()

    merged_output = args.output
    transformed_output = merged_output.replace(".tsv", "_transformed.tsv")

    # Create results directory if it doesn't exist
    Path(merged_output).parent.mkdir(parents=True, exist_ok=True)

    logger.info("Starting variant data merge...")
    logger.info(f"Excel path pattern: {args.excel_path}")
    logger.info(f"Merged output file: {merged_output}")
    logger.info(f"Transformed output file: {transformed_output}")
    logger.info("-" * 50)

    # Merge variants
    merge_all_variants(merged_output, args.excel_path)

    # Transform merged variants
    logger.info("Starting variant data transformation...")
    transform_variants(merged_output, transformed_output)


if __name__ == "__main__":
    main()
