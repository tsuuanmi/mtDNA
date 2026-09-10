#!/usr/bin/env python3
# noqa: N999
"""
Script to merge FASTA sequence information from all sample Excel files into a single FASTA file.
Each sample has an Excel file with FASTA sequence data that needs to be consolidated.
"""

import argparse
import shutil
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


def find_all_excel_files(base_path: str, excel_path: str | None = None) -> list[str]:
    """
    Find all Excel files containing FASTA sequence data.

    Args:
        base_path: Base directory to search in
        excel_path: Custom path pattern for Excel files (optional)

    Returns:
        List of paths to Excel files
    """
    if excel_path:
        # Use custom path if provided
        excel_pattern = str(Path(excel_path) if Path(excel_path).is_absolute() else Path(base_path) / excel_path)
    else:
        # Use default pattern
        excel_pattern = str(Path(base_path) / "data/AllReport/Sub_web/*/2-MT/*-HyperVar-MT.xlsx")

    p = Path(excel_pattern)
    return sorted(str(f) for f in p.parent.glob(p.name))


def is_sequence_valid(sequence: str) -> bool:
    """
    Check if a sequence is valid (not just N's).

    Args:
        sequence: DNA sequence string

    Returns:
        True if sequence contains non-N characters, False otherwise
    """
    # Remove whitespace and count non-N characters
    cleaned_seq = sequence.strip().upper()
    non_n_count = sum(1 for char in cleaned_seq if char != "N")

    # Consider valid if at least 10% of characters are not N
    return non_n_count > len(cleaned_seq) * 0.1


def read_fasta_data(excel_path: str) -> pd.DataFrame:
    """
    Read FASTA sequence data from an Excel file (Fasta sheet).

    Args:
        excel_path: Path to the Excel file

    Returns:
        DataFrame containing FASTA sequence data
    """
    try:
        # Check if 'Fasta' sheet exists
        excel = pd.ExcelFile(excel_path)
        if "Fasta" not in excel.sheet_names:
            logger.debug(f"Skipping {excel_path}: No 'Fasta' sheet found")
            return pd.DataFrame()

        # Read the Excel file, specifically the 'Fasta' sheet
        df = pd.read_excel(excel, sheet_name="Fasta")

        # Ensure required columns are present
        required_columns = ["SampleCode", "Marker", "Consensus"]

        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.debug(f"Skipping {excel_path}: Missing columns: {missing_columns}")
            return pd.DataFrame()

        # Filter out rows with empty consensus sequences
        df = df.dropna(subset=["Consensus"])
        df = df[df["Consensus"].str.strip() != ""]  # type: ignore[reportAttributeAccessIssue]

        # Filter out sequences that are just N's
        df["is_valid"] = df["Consensus"].apply(is_sequence_valid)  # type: ignore[reportAttributeAccessIssue]
        df = df[df["is_valid"]].drop(columns=["is_valid"])  # type: ignore[reportAttributeAccessIssue,reportAttributeAccessIssue]

        if df.empty:
            logger.debug(f"Skipping {excel_path}: No valid sequences found (empty or all N's)")
            return df  # type: ignore[reportReturnType]

        return df[required_columns]  # type: ignore[reportReturnType]

    except (OSError, ValueError, KeyError, pd.errors.ParserError) as e:
        logger.debug(f"Skipping {excel_path}: {e!s}")
        return pd.DataFrame()


def format_marker_for_fasta_header(marker: str) -> str:
    """
    Format marker string for FASTA header.

    Args:
        marker: Original marker string like 'HVS-II[MT:73-340]'

    Returns:
        Formatted marker for FASTA header
    """
    # Extract the main part (e.g., 'HVS-II' from 'HVS-II[MT:73-340]')
    if "[" in marker:
        main_part = marker.split("[", maxsplit=1)[0]
        range_part = marker.split("[")[1].rstrip("]")
        # Convert HVS-II to HV2, HVS-III to HV3, HVS-I to HV1
        if main_part == "HVS-I":
            return "HV1"
        if main_part == "HVS-II":
            return f"HV2|{range_part.split(':')[1]}"
        if main_part == "HVS-III":
            return "HV3"
        return main_part
    return marker


def write_fasta_sequences(sequences: list[dict[str, str]], output_file: str) -> None:
    """
    Write sequences to a FASTA file.

    Args:
        sequences: List of sequence dictionaries with 'header' and 'sequence' keys
        output_file: Output FASTA file path
    """
    with Path(output_file).open("w") as f:
        for seq_data in sequences:
            f.write(f">{seq_data['header']}\n")
            f.write(f"{seq_data['sequence']}\n")


def write_sample_fasta_file(sample_id: str, sequences: list[dict[str, str]], output_dir: str) -> str:
    """
    Write sequences for a single sample to its own FASTA file.

    Args:
        sample_id: ID of the sample
        sequences: List of sequence dictionaries for this sample
        output_dir: Output directory path

    Returns:
        Path to the created FASTA file
    """
    output_file = Path(output_dir) / f"{sample_id}.fasta"
    with output_file.open("w") as f:
        for seq_data in sequences:
            f.write(f">{seq_data['header']}\n")
            f.write(f"{seq_data['sequence']}\n")
    return str(output_file)


def merge_all_fasta_sequences(base_path: str, output_dir: str, excel_path: str | None = None) -> None:
    """
    Extract FASTA sequence information from all samples and create individual FASTA files.

    Args:
        base_path: Base directory containing the data
        output_dir: Output directory path for individual FASTA files
        excel_path: Custom path pattern for Excel files (optional)
    """
    # Clean up existing output directory
    output_dir_path = Path(output_dir)
    if output_dir_path.exists():
        logger.info(f"Removing existing output directory: {output_dir}")
        shutil.rmtree(output_dir_path)

    logger.info("Finding Excel files...")
    excel_files = find_all_excel_files(base_path, excel_path)

    if not excel_files:
        logger.warning("No Excel files found!")
        return

    logger.info(f"Found {len(excel_files)} Excel files")

    # Create fresh output directory
    logger.info("Creating fresh output directory")
    output_dir_path.mkdir(parents=True)

    sample_count = 0
    total_sequences = 0
    created_files = []
    failed_samples = []

    for file_path in excel_files:
        logger.info(f"Processing: {file_path}")

        # Extract sample ID
        sample_id = extract_sample_id_from_path(file_path)

        # Read FASTA sequence data
        fasta_df = read_fasta_data(file_path)

        if fasta_df.empty:
            failed_samples.append(sample_id)
            continue  # Skip silently - detailed debug logs already provided

        sample_sequences = []
        sequence_count = 0

        # Process each sequence for this sample
        for _, row in fasta_df.iterrows():
            marker = row["Marker"]
            consensus = row["Consensus"]

            # Format the FASTA header
            formatted_marker = format_marker_for_fasta_header(marker)  # type: ignore[reportArgumentType]

            sample_sequences.append({"header": formatted_marker, "sequence": consensus})
            sequence_count += 1

        if sample_sequences:
            # Write individual FASTA file for this sample
            output_file = write_sample_fasta_file(sample_id, sample_sequences, output_dir)
            created_files.append(output_file)
            sample_count += 1
            total_sequences += sequence_count

    if not created_files:
        logger.error("No FASTA files created!")
        return

    logger.success("FASTA extraction complete!")
    logger.info(f"Total samples processed successfully: {sample_count}")
    logger.info(f"Total sequences: {total_sequences}")
    logger.info(f"Failed samples: {len(failed_samples)}")
    if failed_samples:
        logger.debug("Failed sample IDs (no FASTA sheet or only N sequences):")
        for sample_id in failed_samples:
            logger.debug(f"  - {sample_id}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Created {len(created_files)} FASTA files")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequence information from all sample Excel files into individual FASTA files.",
    )

    parser.add_argument(
        "-b",
        "--base-path",
        type=str,
        default=Path(__file__).parent.resolve(),
        help="Base directory containing the data (default: current working directory)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output directory path for individual FASTA files (default: {base_path}/results/fasta_sequences/)",
    )

    parser.add_argument(
        "-e",
        "--excel-path",
        type=str,
        help="Custom path pattern for Excel files (default: data/AllReport/Sub_web/*/2-MT/*-HyperVar-MT.xlsx)",
    )

    return parser.parse_args()


def main() -> None:
    """Main function to run the FASTA extraction process."""
    args = parse_arguments()

    base_path = args.base_path

    # Set default output directory if not provided
    output_dir = args.output or str(Path(base_path) / "results" / "fasta_sequences")

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logger.info("Starting FASTA sequence extraction...")
    logger.info(f"Base path: {base_path}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("-" * 50)

    merge_all_fasta_sequences(base_path, output_dir, args.excel_path)


if __name__ == "__main__":
    main()
