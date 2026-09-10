#!/usr/bin/env python3

"""
This script merges the comparison results from FIS-Sanger and FIS-FASTA comparisons
showing ALL variants detected by each method for comprehensive analysis.

The script takes two TSV files as input:
1. fis_sanger_comparison.tsv - Results from comparing FIS NGS with Sanger pipeline
2. fis_fasta_comparison_direct.tsv - Results from comparing FIS NGS with FASTA pipeline

The merged output shows ALL variants detected by FIS and FASTA methods, not limited
to Sanger analyzed intervals, providing a complete view of variant detection across
all sequencing and analysis approaches.

UPDATED MATCH LOGIC:
- Perfect Match calculation now only considers variants within Sanger analyzed intervals
- Variants outside analyzed intervals do not affect Perfect Match status
- Extra Variants shows only variants within analyzed intervals that differ from Sanger
- This provides fair comparison focusing only on regions Sanger actually analyzed

Example usage:
    python merge_comparison_results.py \
        -s results/fis_sanger_comparison.tsv \
        -f results/FIS_FASTA_comparison.tsv \
        -o results/merged_comparison_results.tsv

Output columns (showing all variants):
    - Sanger_Sample_ID: Sanger sample identifier (reference)
    - Sanger_Variants: Variants found by Sanger (reference)
    - Sanger_Analyzed_Intervals: Intervals analyzed by Sanger
    - FIS_Sample_ID: Corresponding FIS NGS sample identifier
    - FIS_All_Variants: ALL variants detected by FIS method (not filtered)
    - FIS_Perfect_Match: Whether FIS exactly matches Sanger within analyzed intervals (Y/N)
      - only considering variants within analyzed intervals
    - FIS_Extra_Variants: Variants detected by FIS within analyzed intervals but not in Sanger
    - FIS_Missed_Variants: Sanger variants missed by FIS within analyzed intervals
    - FASTA_Sample_ID: Corresponding FASTA sample identifier
    - FASTA_All_Variants: ALL variants detected by FASTA method (not filtered)
    - FASTA_Analyzed_Intervals: Intervals analyzed by FASTA method
    - FASTA_Perfect_Match: Whether FASTA exactly matches Sanger within analyzed intervals (Y/N)
      - only considering variants within analyzed intervals
    - FASTA_Extra_Variants: Variants detected by FASTA within analyzed intervals but not in Sanger
    - FASTA_Missed_Variants: Sanger variants missed by FASTA within analyzed intervals
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd
from loguru import logger

from src.core.variants import is_position_in_intervals, normalize_position, parse_variant_position_allele, pos_base


def parse_intervals(intervals_str: str) -> list[tuple[int, int]]:
    """
    Parse analyzed intervals string into list of (start, end) tuples.

    Args:
        intervals_str: String like "16024-16365, 73-340, 438-576"

    Returns:
        List of (start, end) tuples
    """
    if pd.isna(intervals_str) or str(intervals_str).strip() in ["", "None", "N/A"]:
        return []

    intervals = []
    for raw_part in str(intervals_str).split(","):
        part = raw_part.strip()
        if "-" in part:
            try:
                start, end = map(int, part.split("-"))
                intervals.append((start, end))
            except ValueError:
                logger.warning(f"Could not parse interval: {part}")
                continue

    return intervals


def extract_position_from_variant(variant: str) -> int | str | None:
    """Extract a canonical position from simplified variant notation.

    Insertion coordinates retain their decimal suffix so interval filtering can
    distinguish an insertion after an interval endpoint from the endpoint base.
    """
    if not variant or variant in ["None", "", "N/A"]:
        return None

    position, _ = parse_variant_position_allele(variant.strip())
    try:
        normalized = normalize_position(position)
        pos_base(normalized)
    except (TypeError, ValueError):
        return None
    return normalized


def filter_variants_by_intervals(variants_set: set, intervals: list[tuple[int, int]]) -> set:
    """
    Filter variants to only include those within analyzed intervals.

    Args:
        variants_set: Set of variant strings
        intervals: List of (start, end) tuples for analyzed intervals

    Returns:
        Set of variants within intervals
    """
    if not intervals:
        return variants_set

    filtered_variants = set()
    for variant in variants_set:
        position = extract_position_from_variant(variant)
        if position is not None and is_position_in_intervals(position, intervals):
            filtered_variants.add(variant)

    return filtered_variants


def arg_parser() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(
        description="Merge FIS-Sanger and FIS-FASTA comparison results showing ALL variants detected",
    )
    parser.add_argument("-s", "--sanger_tsv", type=Path, required=True, help="Path to fis_sanger_comparison.tsv file")
    parser.add_argument("-f", "--fasta_tsv", type=Path, required=True, help="Path to FIS_FASTA_comparison.tsv file")
    parser.add_argument(
        "-o",
        "--output_file",
        type=Path,
        required=True,
        help="Path to output merged TSV file (showing all variants)",
    )
    return parser.parse_args()


def load_comparison_data(sanger_tsv: Path, fasta_tsv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load both comparison TSV files.

    Args:
        sanger_tsv: Path to Sanger comparison TSV
        fasta_tsv: Path to FASTA comparison TSV

    Returns:
        Tuple of (sanger_df, fasta_df)
    """
    logger.info(f"Loading Sanger comparison data from {sanger_tsv}")
    sanger_df = pd.read_csv(sanger_tsv, sep="\t")
    logger.info(f"Loaded {len(sanger_df)} Sanger comparison records")

    logger.info(f"Loading FASTA comparison data from {fasta_tsv}")
    fasta_df = pd.read_csv(fasta_tsv, sep="\t")
    logger.info(f"Loaded {len(fasta_df)} FASTA comparison records")

    return sanger_df, fasta_df


def merge_comparison_results(sanger_df: pd.DataFrame, fasta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge the Sanger and FASTA comparison results showing ALL variants detected by each method.

    Perfect Match Logic:
    - FIS/FASTA Perfect Match = 'Y' only if the method detects EXACTLY the same variants
      as Sanger within the analyzed intervals (no extra variants, no missing variants)
    - Only variants within Sanger analyzed intervals are considered for comparison
    - Variants outside analyzed intervals do not affect the perfect match status

    Args:
        sanger_df: DataFrame with Sanger comparison results (FIS vs Sanger)
        fasta_df: DataFrame with FASTA comparison results (FIS vs FASTA)

    Returns:
        DataFrame with comprehensive comparison results showing all variants
    """
    logger.info("Merging comparison results showing ALL variants from FIS and FASTA...")

    # Start with Sanger data as the base (reference)
    base_df = sanger_df.copy()

    # Create mapping from FIS_Sample to match with FASTA data
    fasta_lookup = fasta_df.set_index("FIS_Sample")

    # Initialize result DataFrame with updated column names

    merged_df = pd.DataFrame()

    # Process each Sanger sample as the reference
    for _, sanger_row in base_df.iterrows():
        fis_sample = sanger_row["FIS_Sample"]

        # Get Sanger variants and analyzed intervals
        sanger_variants_set = (
            set(str(sanger_row["Sanger_Variants"]).split()) if pd.notna(sanger_row["Sanger_Variants"]) else set()  # type: ignore[reportGeneralTypeIssues]  # type: ignore[reportGeneralTypeIssues]
        )
        sanger_intervals = parse_intervals(str(sanger_row["Analyzed_Intervals"]))

        # Get ALL FIS variants
        fis_variants_set = (
            set(str(sanger_row["FIS_Variants"]).split()) if pd.notna(sanger_row["FIS_Variants"]) else set()  # type: ignore[reportGeneralTypeIssues]
        )

        # Remove 'None' and empty strings
        sanger_variants_set.discard("None")
        sanger_variants_set.discard("")
        fis_variants_set.discard("None")
        fis_variants_set.discard("")

        # Filter FIS variants to only those within Sanger analyzed intervals for comparison
        fis_variants_in_intervals = filter_variants_by_intervals(fis_variants_set, sanger_intervals)

        # Calculate FIS concordance and differences (only within analyzed intervals)
        fis_extra_in_intervals = fis_variants_in_intervals - sanger_variants_set
        fis_missing = sanger_variants_set - fis_variants_in_intervals

        # Determine FIS perfect match status (only considers variants within analyzed intervals)
        # Perfect match means: no missing variants AND no extra variants within analyzed intervals
        fis_perfect_match = "Y" if len(fis_missing) == 0 and len(fis_extra_in_intervals) == 0 else "N"

        # Create base row with updated column names showing ALL FIS variants
        result_row = {
            "Sanger_Sample_ID": sanger_row["Sanger_Sample"],
            "Sanger_Batch": sanger_row["Sanger_Batch"],
            "Sanger_Variants": sanger_row["Sanger_Variants"],
            "Sanger_Analyzed_Intervals": sanger_row["Analyzed_Intervals"],
            "FIS_Sample_ID": fis_sample,
            "FIS_All_Variants": sanger_row["FIS_Variants"],  # Show ALL FIS variants
            "FIS_Perfect_Match": fis_perfect_match,
            "FIS_Extra_Variants": " ".join(sorted(fis_extra_in_intervals)) if fis_extra_in_intervals else "None",
            "FIS_Missed_Variants": " ".join(sorted(fis_missing)) if fis_missing else "None",
        }

        # Add FASTA data if available for this FIS sample
        if fis_sample in fasta_lookup.index:
            fasta_row = fasta_lookup.loc[fis_sample]

            # Get FASTA analyzed intervals (if available)
            parse_intervals(str(fasta_row["Analyzed_Intervals"]) if "Analyzed_Intervals" in fasta_row else "")

            # For FASTA comparison, use Sanger intervals as reference (same as FIS comparison)
            # This ensures consistent comparison criteria
            comparison_intervals = sanger_intervals

            # Calculate FASTA vs Sanger comparison
            sanger_variants_set = (
                set(str(sanger_row["Sanger_Variants"]).split()) if pd.notna(sanger_row["Sanger_Variants"]) else set()  # type: ignore[reportGeneralTypeIssues]  # type: ignore[reportGeneralTypeIssues]
            )
            # Get ALL FASTA variants
            fasta_variants_set = (
                set(str(fasta_row["fasta_Variants"]).split()) if pd.notna(fasta_row["fasta_Variants"]) else set()
            )

            # Remove 'None' and empty strings
            sanger_variants_set.discard("None")
            sanger_variants_set.discard("")
            fasta_variants_set.discard("None")
            fasta_variants_set.discard("")

            # Filter FASTA variants to only those within Sanger analyzed intervals for comparison
            fasta_variants_in_intervals = filter_variants_by_intervals(fasta_variants_set, comparison_intervals)

            # Calculate concordance and differences (only within analyzed intervals)
            fasta_extra_in_intervals = fasta_variants_in_intervals - sanger_variants_set
            fasta_missing = sanger_variants_set - fasta_variants_in_intervals

            # Determine FASTA perfect match status (only considers variants within analyzed intervals)
            # Perfect match means: no missing variants AND no extra variants within analyzed intervals
            fasta_perfect_match = "Y" if len(fasta_missing) == 0 and len(fasta_extra_in_intervals) == 0 else "N"

            result_row.update(
                {
                    "FASTA_Sample_ID": fasta_row["fasta_Sample"],
                    "FASTA_All_Variants": fasta_row["fasta_Variants"],  # Show ALL FASTA variants
                    "FASTA_Analyzed_Intervals": fasta_row.get("Analyzed_Intervals", "N/A"),
                    "FASTA_Perfect_Match": fasta_perfect_match,
                    "FASTA_Extra_Variants": " ".join(sorted(fasta_extra_in_intervals))
                    if fasta_extra_in_intervals
                    else "None",
                    "FASTA_Missed_Variants": " ".join(sorted(fasta_missing)) if fasta_missing else "None",
                },
            )
        else:
            # No FASTA data available for this sample
            result_row.update(
                {
                    "FASTA_Sample_ID": "N/A",
                    "FASTA_All_Variants": "N/A",
                    "FASTA_Analyzed_Intervals": "N/A",
                    "FASTA_Perfect_Match": "N/A",
                    "FASTA_Extra_Variants": "N/A",
                    "FASTA_Missed_Variants": "N/A",
                },
            )

        # Add row to result DataFrame
        merged_df = pd.concat([merged_df, pd.DataFrame([result_row])], ignore_index=True)

    # Reorder columns for better readability (showing all variants)
    column_order = [
        "Sanger_Sample_ID",
        "Sanger_Variants",
        "Sanger_Analyzed_Intervals",
        "FIS_Sample_ID",
        "FIS_All_Variants",
        "FIS_Perfect_Match",
        "FIS_Extra_Variants",
        "FIS_Missed_Variants",
        "FASTA_Sample_ID",
        "FASTA_All_Variants",
        "FASTA_Analyzed_Intervals",
        "FASTA_Perfect_Match",
        "FASTA_Extra_Variants",
        "FASTA_Missed_Variants",
    ]

    # Select only existing columns in the specified order
    existing_columns = [col for col in column_order if col in merged_df.columns]
    merged_df = merged_df[existing_columns]

    # Sort by Sanger_Sample_ID for consistent output
    merged_df = merged_df.sort_values("Sanger_Sample_ID")  # type: ignore[reportCallIssue]

    logger.info(f"Merged data contains {len(merged_df)} records (showing all variants)")

    return merged_df


def generate_summary_statistics(merged_df: pd.DataFrame) -> dict[str, int]:
    """
    Generate summary statistics for the merged comparison results showing all variants.

    Args:
        merged_df: Merged comparison DataFrame showing all variants

    Returns:
        Dictionary with summary statistics
    """
    stats = {}

    # Total samples with Sanger data (reference)
    stats["Total_Sanger_Samples"] = len(merged_df)

    # FIS vs Sanger comparison stats
    fis_available = merged_df["FIS_Sample_ID"].notna().sum()
    stats["FIS_vs_Sanger_Available"] = fis_available

    if fis_available > 0:
        fis_perfect_match = (merged_df["FIS_Perfect_Match"] == "Y").sum()
        stats["FIS_vs_Sanger_Perfect_Match"] = fis_perfect_match

    # FASTA vs Sanger comparison stats
    fasta_available = (merged_df["FASTA_Sample_ID"] != "N/A").sum()
    stats["FASTA_vs_Sanger_Available"] = fasta_available

    if fasta_available > 0:
        fasta_perfect_match = (merged_df["FASTA_Perfect_Match"] == "Y").sum()
        stats["FASTA_vs_Sanger_Perfect_Match"] = fasta_perfect_match

    # Both methods available for comparison
    both_available = ((merged_df["FIS_Sample_ID"].notna()) & (merged_df["FASTA_Sample_ID"] != "N/A")).sum()
    stats["Both_Methods_Available"] = both_available

    return stats


def save_merged_results(merged_df: pd.DataFrame, output_file: Path, stats: dict[str, int]) -> None:
    """
    Save the merged results showing all variants to TSV file and log summary statistics.

    Args:
        merged_df: Merged comparison DataFrame showing all variants
        output_file: Path to output TSV file
        stats: Summary statistics dictionary
    """
    logger.info(f"Saving merged results (all variants) to {output_file}")
    merged_df.to_csv(output_file, sep="\t", index=False)

    # Log summary statistics
    logger.info("=== MERGE SUMMARY STATISTICS (ALL VARIANTS) ===")
    logger.info(f"Total Sanger samples (reference): {stats['Total_Sanger_Samples']}")
    logger.info(f"Samples with FIS data: {stats['FIS_vs_Sanger_Available']}")
    logger.info(f"Samples with FASTA data: {stats['FASTA_vs_Sanger_Available']}")
    logger.info(f"Samples with both FIS and FASTA data: {stats['Both_Methods_Available']}")

    if stats["FIS_vs_Sanger_Available"] > 0:
        logger.info(f"FIS vs Sanger perfect match: {stats['FIS_vs_Sanger_Perfect_Match']}")

    if stats["FASTA_vs_Sanger_Available"] > 0:
        logger.info(f"FASTA vs Sanger perfect match: {stats['FASTA_vs_Sanger_Perfect_Match']}")

    # Calculate concordance rates
    if stats["FIS_vs_Sanger_Available"] > 0:
        fis_concordance_rate = (stats["FIS_vs_Sanger_Perfect_Match"] / stats["FIS_vs_Sanger_Available"]) * 100
        logger.info(f"FIS vs Sanger perfect match rate: {fis_concordance_rate:.1f}%")

    if stats["FASTA_vs_Sanger_Available"] > 0:
        fasta_concordance_rate = (stats["FASTA_vs_Sanger_Perfect_Match"] / stats["FASTA_vs_Sanger_Available"]) * 100
        logger.info(f"FASTA vs Sanger perfect match rate: {fasta_concordance_rate:.1f}%")


def merge_comparison_files(sanger_tsv: Path, fasta_tsv: Path, output_file: Path) -> None:
    """
    Main function to merge FIS-Sanger and FIS-FASTA comparison results showing ALL variants.

    Args:
        sanger_tsv: Path to Sanger comparison TSV
        fasta_tsv: Path to FASTA comparison TSV
        output_file: Path to output merged TSV file showing all variants
    """
    try:
        # Load comparison data
        sanger_df, fasta_df = load_comparison_data(sanger_tsv, fasta_tsv)

        # Merge the results showing all variants
        merged_df = merge_comparison_results(sanger_df, fasta_df)

        # Generate summary statistics
        stats = generate_summary_statistics(merged_df)

        # Save results
        save_merged_results(merged_df, output_file, stats)

        logger.success(f"Successfully merged comparison results showing all variants and saved to {output_file}")

    except Exception as e:
        logger.error(f"Error merging comparison results: {e}")
        raise


def main() -> None:
    """Main entry point."""
    args = arg_parser()

    logger.info("Starting merge of FIS-Sanger and FIS-FASTA comparison results showing ALL variants")

    merge_comparison_files(sanger_tsv=args.sanger_tsv, fasta_tsv=args.fasta_tsv, output_file=args.output_file)


if __name__ == "__main__":
    main()
