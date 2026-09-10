#!/usr/bin/env python3

"""
This script compares variant analysis results between FIS NGS pipeline
(merged_variants_transformed.tsv) and FASTA pipeline (statistic_fullbatch.json).
It generates a TSV file to facilitate quality control by identifying concordance
between the two analysis methods.

This is a simplified version that assumes direct sample name matching between
the FIS and FASTA datasets, eliminating the need for complex ID mapping chains.

The script includes ENHANCED CONCORDANCE ANALYSIS that considers:
1. IUPAC ambiguity codes (e.g., 199C matches 199Y because Y includes C and T)
2. Alignment-based concordance for specific cases:
   - 525.1A 525.2C matches 524.1A 524.2C (insertions shifted by 1)
   - 524T matches 522T (SNP shifted by 2)
   - 494DEL matches 498DEL (deletion shifted by 4)
   - 16189DEL matches 16183C 16184A (complex alignment)

The script shows ALL FIS variants detected by the NGS pipeline, not filtered to
FASTA analyzed intervals. This provides a comprehensive view of FIS detection
capabilities while comparing concordance within FASTA's analyzed regions.

The script extracts the actual analyzed intervals from each FASTA sample
(from the 'intervals' field in statistic_fullbatch.json) and uses these intervals
for concordance comparison, but shows all FIS variants in the output for
complete analysis coverage.

Example usage:
    python compare_fis_fasta_direct.py \\
        -f results/merged_variants_transformed.tsv \\
        -s temp/statistic_fullbatch.json \\
        -o temp/fis_fasta_comparison_direct.tsv

Output columns:
    - FIS_Sample: FIS NGS sample identifier
    - fasta_Sample: Corresponding FASTA sample ID (should match FIS_Sample)
    - fasta_Batch: Batch information from FASTA analysis
    - Analyzed_Intervals: Actual intervals analyzed by FASTA
    - FIS_Variants: ALL variants found by FIS NGS (not filtered by analyzed intervals)
    - fasta_Variants: Variants found by FASTA (within analyzed intervals)
    - Concordant: Whether the variant profiles match using exact matching (Y/N)
    - Concordant_Enhanced: Whether the variant profiles match using IUPAC codes and alignment rules (Y/N)
    - FIS_Unique: Variants found only by FIS NGS (original method)
    - fasta_Unique: Variants found only by FASTA (original method)
    - FIS_Unique_Enhanced: Variants found only by FIS NGS (enhanced method)
    - fasta_Unique_Enhanced: Variants found only by FASTA (enhanced method)
    - FIS_Count: Number of variants found by FIS NGS
    - fasta_Count: Number of variants found by FASTA
    - FIS_Flag: Whether FIS variants contain problematic patterns (Y/N)
    - FIS_Flag_Reasons: Specific reasons for flagging
"""

# Standard library
import json
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any

# Third-party
import pandas as pd
from loguru import logger

from src.core.variants import (
    format_variants_list,
    is_position_in_intervals,
    iupac_to_bases,
    normalize_position,
    parse_variant_position_allele,
    pos_base,
    pos_sort_key,
)

_INTERVAL_PAIR_LEN = 2
_HV_REGION_START = 16180
_HV_REGION_END = 16193


def arg_parser() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Compare variants between FIS NGS and FASTA pipelines (direct matching)")
    parser.add_argument(
        "-f",
        "--fis_tsv",
        type=Path,
        required=True,
        help="Path to merged_variants_transformed.tsv file (FIS NGS)",
    )
    parser.add_argument(
        "-s",
        "--fasta_json",
        type=Path,
        required=True,
        help="Path to statistic_fullbatch.json file (FASTA)",
    )
    parser.add_argument("-o", "--output_file", type=Path, required=True, help="Path to output TSV file")
    parser.add_argument(
        "--separator",
        type=str,
        default="space",
        choices=["comma", "space"],
        help="Separator for variant listings: comma or space. Default is space",
    )
    return parser.parse_args()


def parse_genotypes_from_tsv_string(genotype_string: object | None) -> list[str]:
    """
    Parse genotypes from TSV string format.

    Args:
        genotype_string: Space-separated genotypes like "152C 263G 315.1C 477C 16519C"

    Returns:
        List of genotype strings
    """
    # Handle NaN, None, or empty values
    if genotype_string is None or pd.isna(genotype_string):  # type: ignore[reportGeneralTypeIssues]
        return []

    # Convert to string if it's not already
    genotype_str = str(genotype_string).strip()

    if not genotype_str or genotype_str == "":
        return []

    # Split by space and clean up
    return [g.strip() for g in genotype_str.split() if g.strip()]


def get_fasta_intervals(sample_data: dict[str, Any]) -> list[tuple]:  # noqa: C901
    """
    Extract analyzed intervals from FASTA sample data.

    Args:
        sample_data: Sample data dictionary from FASTA JSON

    Returns:
        List of (start, end) tuples representing analyzed intervals
    """
    intervals = []

    if "intervals" in sample_data:
        intervals_data = sample_data["intervals"]

        # Handle different interval formats
        if isinstance(intervals_data, dict):
            for region_intervals in intervals_data.values():
                if isinstance(region_intervals, list):
                    for interval in region_intervals:
                        if isinstance(interval, list) and len(interval) == _INTERVAL_PAIR_LEN:
                            try:
                                start, end = int(interval[0]), int(interval[1])
                                intervals.append((start, end))
                            except (ValueError, TypeError):
                                logger.warning(f"Could not parse interval: {interval}")
        elif isinstance(intervals_data, list):
            # Format: ["73-340", "16024-16365"] (string format)
            for interval_str in intervals_data:
                try:
                    # Parse interval like "73-340" or "16024-16365"
                    start, end = map(int, interval_str.split("-"))
                    intervals.append((start, end))
                except (ValueError, AttributeError):
                    logger.warning(f"Could not parse interval: {interval_str}")

    return intervals


def get_all_variants_from_fasta_with_intervals(sample_data: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Get all variants (SNPs, insertions, deletions) from FASTA sample data.
    Only includes variants within the actual analyzed intervals for this sample.

    Args:
        sample_data: Sample data dictionary from FASTA JSON

    Returns:
        Combined list of all variants with normalized format, filtered by actual intervals
    """
    all_variants = []

    # Get the actual analyzed intervals for this sample
    intervals = get_fasta_intervals(sample_data)

    if "variants" in sample_data:
        variants_data = sample_data["variants"]

        # Process SNPs
        if "snps" in variants_data:
            all_variants.extend(
                {
                    "pos": variant["pos"],
                    "ref": variant.get("ref", ""),
                    "seq": variant.get("seq", ""),
                    "type": "snp",
                }
                for variant in variants_data["snps"]
                if is_position_in_intervals(variant["pos"], intervals)
            )

        # Process insertions
        if "insertions" in variants_data:
            all_variants.extend(
                {
                    "pos": variant["pos"],
                    "ref": variant.get("ref", ""),
                    "seq": variant.get("seq", ""),
                    "type": "insertion",
                }
                for variant in variants_data["insertions"]
                if is_position_in_intervals(variant["pos"], intervals)
            )

        # Process deletions
        if "deletions" in variants_data:
            all_variants.extend(
                {
                    "pos": variant["pos"],
                    "ref": variant.get("ref", ""),
                    "seq": variant.get("seq", ""),
                    "type": "deletion",
                }
                for variant in variants_data["deletions"]
                if is_position_in_intervals(variant["pos"], intervals)
            )

    return all_variants


def filter_fis_variants_by_intervals(genotypes: list[str], intervals: list[tuple]) -> list[str]:
    """
    Filter FIS genotypes to only include those within the specified intervals.

    Args:
        genotypes: List of genotype strings (e.g., ["152C", "263G", "315.1C"])
        intervals: List of (start, end) tuples representing analyzed intervals

    Returns:
        Filtered list of genotypes within the intervals
    """
    if not intervals:
        return genotypes  # If no intervals specified, return all genotypes

    filtered_genotypes = []

    for genotype in genotypes:
        # Extract position from genotype
        pos_str = ""
        for _, char in enumerate(genotype):
            if char.isdigit() or char == ".":
                pos_str += char
            else:
                break

        if pos_str and is_position_in_intervals(pos_str, intervals):
            filtered_genotypes.append(genotype)

    return filtered_genotypes


def convert_variants_to_genotypes(variants: list[dict[str, Any]]) -> list[str]:
    """
    Convert variant list to genotype format.

    Args:
        variants: List of variant dictionaries

    Returns:
        List of genotype strings
    """
    genotypes = []

    for variant in sorted(variants, key=lambda x: pos_sort_key(x["pos"])):
        pos = str(normalize_position(variant["pos"]))

        # Check if it's a deletion (seq is "-" or type is "deletion")
        if variant.get("type") == "deletion" or variant.get("seq") == "-":
            genotype = f"{pos}DEL"
        else:
            # SNP or insertion
            seq = variant.get("seq", variant.get("alt", ""))
            genotype = f"{pos}{seq}"

        genotypes.append(genotype)

    return genotypes


def find_unique_genotypes(genotypes_a: list[str], genotypes_b: list[str]) -> list[str]:
    """
    Find genotypes that are unique to genotypes_a (not present in genotypes_b).

    Args:
        genotypes_a: First list of genotypes
        genotypes_b: Second list of genotypes

    Returns:
        List of genotypes unique to genotypes_a
    """
    if not genotypes_b:
        return genotypes_a

    set_a = set(genotypes_a)
    set_b = set(genotypes_b)

    unique_genotypes = list(set_a - set_b)
    return sorted(unique_genotypes)


def is_concordant_genotypes(genotypes_a: list[str], genotypes_b: list[str]) -> str:
    """
    Check if two genotype lists are concordant (same genotypes).

    Args:
        genotypes_a: Genotypes from first source
        genotypes_b: Genotypes from second source

    Returns:
        "Y" if concordant, "N" if not concordant, "N/A" if one is empty
    """
    # If either has no genotypes, return N/A
    if not genotypes_a and not genotypes_b:
        return "Y"  # Both empty is concordant
    if not genotypes_a or not genotypes_b:
        return "N/A"

    # Convert to sets for comparison
    set_a = set(genotypes_a)
    set_b = set(genotypes_b)

    # Check if sets are equal
    if set_a == set_b:
        return "Y"
    return "N"


def is_iupac_concordant(allele1: str, allele2: str) -> bool:
    """
    Check if two alleles are concordant considering IUPAC codes.

    Args:
        allele1: First allele (e.g., "C")
        allele2: Second allele (e.g., "Y")

    Returns:
        True if alleles are concordant considering IUPAC codes
    """
    # Handle deletions
    if allele1 == "DEL" or allele2 == "DEL":
        return allele1 == allele2

    # Get possible bases for each allele
    bases1 = set(iupac_to_bases(allele1))
    bases2 = set(iupac_to_bases(allele2))

    # If either is empty, they're not concordant
    if not bases1 or not bases2:
        return False

    # Check if there's any overlap between the possible bases
    return bool(bases1 & bases2)


def is_alignment_concordant(variant1: str, variant2: str) -> bool:
    """
    Check if two variants are concordant based on alignment rules.

    Args:
        variant1: First variant (e.g., "525.1A")
        variant2: Second variant (e.g., "524.1A")

    Returns:
        True if variants are concordant based on alignment rules
    """
    pos1, allele1 = parse_variant_position_allele(variant1)
    pos2, allele2 = parse_variant_position_allele(variant2)

    # Define alignment rules
    alignment_rules = [
        (("525.1", "A"), ("524.1", "A")),
        (("525.2", "C"), ("524.2", "C")),
        (("525.3", "A"), ("524.3", "A")),
        (("525.4", "C"), ("524.4", "C")),
        (("524", "T"), ("522", "T")),
        (("494", "DEL"), ("498", "DEL")),
    ]

    # Check direct alignment rules
    for (rule_pos1, rule_allele1), (rule_pos2, rule_allele2) in alignment_rules:
        if (pos1 == rule_pos1 and allele1 == rule_allele1 and pos2 == rule_pos2 and allele2 == rule_allele2) or (
            pos1 == rule_pos2 and allele1 == rule_allele2 and pos2 == rule_pos1 and allele2 == rule_allele1
        ):
            return True

    return False


def check_16189_deletion_alignment(variants1: list[str], variants2: list[str]) -> bool:
    """
    Check if 16189DEL in one set matches 16183C 16184A in another set.

    Args:
        variants1: First set of variants
        variants2: Second set of variants

    Returns:
        True if the alignment rule is satisfied
    """
    has_16189_del_1 = "16189DEL" in variants1
    has_16189_del_2 = "16189DEL" in variants2

    has_16183c_16184a_1 = "16183C" in variants1 and "16184A" in variants1
    has_16183c_16184a_2 = "16183C" in variants2 and "16184A" in variants2

    # Check if 16189DEL in one matches 16183C+16184A in the other
    return (has_16189_del_1 and has_16183c_16184a_2) or (has_16189_del_2 and has_16183c_16184a_1)


def is_enhanced_concordant_genotypes(genotypes_a: list[str], genotypes_b: list[str]) -> str:  # noqa: C901,PLR0912
    """
    Check if two genotype lists are concordant using enhanced matching (IUPAC + alignment).

    Args:
        genotypes_a: Genotypes from first source
        genotypes_b: Genotypes from second source

    Returns:
        "Y" if concordant, "N" if not concordant, "N/A" if one is empty
    """
    # If either has no genotypes, return N/A
    if not genotypes_a and not genotypes_b:
        return "Y"  # Both empty is concordant
    if not genotypes_a or not genotypes_b:
        return "N/A"

    # Start with copies of the lists
    remaining_a = genotypes_a.copy()
    remaining_b = genotypes_b.copy()

    # Remove exact matches first
    for variant in genotypes_a:
        if variant in remaining_b:
            remaining_a.remove(variant)
            remaining_b.remove(variant)

    # Remove IUPAC matches
    for variant_a in remaining_a.copy():
        pos_a, allele_a = parse_variant_position_allele(variant_a)

        for variant_b in remaining_b.copy():
            pos_b, allele_b = parse_variant_position_allele(variant_b)

            if pos_a == pos_b and is_iupac_concordant(allele_a, allele_b):
                remaining_a.remove(variant_a)
                remaining_b.remove(variant_b)
                break

    # Remove alignment matches
    for variant_a in remaining_a.copy():
        for variant_b in remaining_b.copy():
            if is_alignment_concordant(variant_a, variant_b):
                remaining_a.remove(variant_a)
                remaining_b.remove(variant_b)
                break

    # Handle special case for 16189DEL vs 16183C+16184A
    if check_16189_deletion_alignment(genotypes_a, genotypes_b):
        if "16189DEL" in remaining_a:
            remaining_a.remove("16189DEL")
        if "16189DEL" in remaining_b:
            remaining_b.remove("16189DEL")
        # Remove the matching variants from the other set
        for variant in ["16183C", "16184A"]:
            if variant in remaining_a:
                remaining_a.remove(variant)
            if variant in remaining_b:
                remaining_b.remove(variant)

    # If both lists are empty after all matching, they're concordant
    if not remaining_a and not remaining_b:
        return "Y"
    return "N"


def find_unique_genotypes_enhanced(genotypes_a: list[str], genotypes_b: list[str]) -> list[str]:  # noqa: C901,PLR0912
    """
    Find genotypes unique to genotypes_a using enhanced matching (IUPAC + alignment).

    Args:
        genotypes_a: First list of genotypes
        genotypes_b: Second list of genotypes

    Returns:
        List of genotypes unique to genotypes_a after enhanced matching
    """
    if not genotypes_b:
        return genotypes_a

    # Start with copies of the lists
    remaining_a = genotypes_a.copy()
    remaining_b = genotypes_b.copy()

    # Remove exact matches first
    for variant in genotypes_a:
        if variant in remaining_b:
            remaining_a.remove(variant)
            remaining_b.remove(variant)

    # Remove IUPAC matches
    for variant_a in remaining_a.copy():
        pos_a, allele_a = parse_variant_position_allele(variant_a)

        for variant_b in remaining_b.copy():
            pos_b, allele_b = parse_variant_position_allele(variant_b)

            if pos_a == pos_b and is_iupac_concordant(allele_a, allele_b):
                remaining_a.remove(variant_a)
                remaining_b.remove(variant_b)
                break

    # Remove alignment matches
    for variant_a in remaining_a.copy():
        for variant_b in remaining_b.copy():
            if is_alignment_concordant(variant_a, variant_b):
                remaining_a.remove(variant_a)
                remaining_b.remove(variant_b)
                break

    # Handle special case for 16189DEL vs 16183C+16184A
    if check_16189_deletion_alignment(genotypes_a, genotypes_b):
        if "16189DEL" in remaining_a:
            remaining_a.remove("16189DEL")
        # Remove matching variants
        for variant in ["16183C", "16184A"]:
            if variant in remaining_a:
                remaining_a.remove(variant)

    return remaining_a


def should_flag_sample(genotypes: list[str]) -> tuple[bool, list[str]]:  # noqa: C901
    """
    Check if sample contains FIS variants that should trigger flagging.
    This function examines FIS-generated variants to help validate pipeline performance.

    Args:
        genotypes: List of FIS-generated genotype strings

    Returns:
        Tuple of (should_flag, flag_reasons)
    """
    if not genotypes:
        return False, []

    flag_reasons = []
    region_16180_16193_positions = []

    # Convert genotypes to variant format for analysis
    variants = []
    for genotype in genotypes:
        try:
            pos_str = ""
            for _, char in enumerate(genotype):
                if char.isdigit() or char == ".":
                    pos_str += char
                else:
                    break

            if pos_str:
                pos = normalize_position(pos_str)
                allele = genotype[len(pos_str) :]

                variants.append({"pos": pos, "allele": allele, "original": genotype})
        except (ValueError, IndexError):
            continue

    # Check for problematic variants with specific exceptions
    for v in variants:
        pos = v["pos"]
        allele = v["allele"]

        # Track positions in the 16180-16193 region
        if _HV_REGION_START <= pos_base(pos) <= _HV_REGION_END:
            region_16180_16193_positions.append(v["original"])

        # Flag insertions at position 16182 (but not others in this region)
        if normalize_position(pos) == "16182.1" and allele != "DEL":
            flag_reasons.append(f"Insertion at 16182: {v['original']}")

    if region_16180_16193_positions:
        flag_reasons.append(f"16180-16193 region variants: {', '.join(region_16180_16193_positions)}")

    return len(flag_reasons) > 0, flag_reasons


def filter_excluded_variants(genotypes: list[str]) -> list[str]:
    """
    Filter out variants that should be excluded from comparison.
    Currently excludes: 309.1, 309.2, 309.3 variants.

    Args:
        genotypes: List of genotype strings

    Returns:
        Filtered list of genotypes with excluded variants removed
    """
    excluded_positions = ["309.1", "309.2", "309.3"]
    filtered_genotypes = []

    for genotype in genotypes:
        # Extract position from genotype
        pos_str = ""
        for _, char in enumerate(genotype):
            if char.isdigit() or char == ".":
                pos_str += char
            else:
                break

        # Only include if position is not in excluded list
        if pos_str not in excluded_positions:
            filtered_genotypes.append(genotype)

    return filtered_genotypes


def compare_fis_fasta_direct(  # noqa: PLR0915
    fis_tsv: Path,
    fasta_json: Path,
    output_file: Path,
    variant_separator: str = "space",
) -> None:
    """
    Compare variants between FIS NGS and FASTA pipelines using direct sample matching.

    Args:
        fis_tsv: Path to merged_variants_transformed.tsv (FIS NGS)
        fasta_json: Path to statistic_fullbatch.json (FASTA)
        output_file: Path to output TSV file
        variant_separator: Separator for variant display
    """
    logger.info("Starting FIS NGS vs FASTA comparison (direct matching)")

    # Load FIS TSV file
    logger.info(f"Loading FIS NGS variants from {fis_tsv}")
    fis_data = pd.read_csv(fis_tsv, sep="\t")

    # Load FASTA JSON file
    logger.info(f"Loading FASTA variants from {fasta_json}")
    with Path(fasta_json).open() as f:
        fasta_data = json.load(f)

    # Get sample lists
    fis_samples = fis_data["Sample_ID"].tolist()
    fasta_samples = list(fasta_data.keys())

    logger.info(f"Found {len(fis_samples)} samples in FIS and {len(fasta_samples)} samples in FASTA")

    # Prepare comparison results
    comparison_results = []

    # Process each sample in FIS data
    for _, fis_row in fis_data.iterrows():
        fis_sample_id = fis_row["Sample_ID"]
        fis_genotypes_str = fis_row.get("Genotype", "")

        # Parse FIS genotypes
        fis_genotypes = parse_genotypes_from_tsv_string(fis_genotypes_str)

        # Look for matching FASTA sample (direct matching)
        fasta_sample_id = fis_sample_id if fis_sample_id in fasta_data else "Not found"

        if fasta_sample_id != "Not found":
            # Get the actual analyzed intervals for this FASTA sample
            fasta_intervals = get_fasta_intervals(fasta_data[fasta_sample_id])

            # Keep ALL FIS genotypes (do not filter by FASTA intervals)
            fis_genotypes_all = fis_genotypes  # Show ALL FIS variants

            # For concordance comparison, filter FIS variants to match FASTA analyzed intervals
            fis_genotypes_filtered = filter_fis_variants_by_intervals(fis_genotypes, fasta_intervals)

            # Get FASTA variants (already filtered by intervals)
            fasta_variants = get_all_variants_from_fasta_with_intervals(fasta_data[fasta_sample_id])
            fasta_genotypes = convert_variants_to_genotypes(fasta_variants)

            # Filter out excluded variants (309.1, 309.2, 309.3) from filtered FIS and FASTA
            fis_genotypes_filtered = filter_excluded_variants(fis_genotypes_filtered)
            fasta_genotypes = filter_excluded_variants(fasta_genotypes)

            # Also filter excluded variants from ALL FIS variants for display
            fis_genotypes_all = filter_excluded_variants(fis_genotypes_all)

            # Format intervals for display
            interval_ranges = []
            for start, end in fasta_intervals:
                interval_ranges.append(f"{start}-{end}")
            interval_str = ", ".join(interval_ranges) if interval_ranges else "None"

            # Extract batch information from FASTA data
            batch_info = fasta_data[fasta_sample_id].get("batch", "Unknown")
        else:
            # No FASTA match - show all FIS variants
            fis_genotypes_all = filter_excluded_variants(fis_genotypes)
            fis_genotypes_filtered = []
            fasta_genotypes = []
            interval_str = "N/A"
            batch_info = "N/A"

        # Compare variants - both original and enhanced methods
        concordant = is_concordant_genotypes(fis_genotypes_filtered, fasta_genotypes)
        concordant_enhanced = is_enhanced_concordant_genotypes(fis_genotypes_filtered, fasta_genotypes)

        # Find unique variants
        fis_unique = find_unique_genotypes(fis_genotypes_filtered, fasta_genotypes)
        fasta_unique = find_unique_genotypes(fasta_genotypes, fis_genotypes_filtered)

        fis_unique_enhanced = find_unique_genotypes_enhanced(fis_genotypes_filtered, fasta_genotypes)
        fasta_unique_enhanced = find_unique_genotypes_enhanced(fasta_genotypes, fis_genotypes_filtered)

        # Flag analysis for FIS variants (use ALL FIS variants for flagging)
        should_flag, flag_reasons = should_flag_sample(fis_genotypes_all)
        flag_status = "Y" if should_flag else "N"
        flag_reasons_str = "; ".join(flag_reasons) if flag_reasons else "None"

        # Format results
        result = {
            "FIS_Sample": fis_sample_id,
            "fasta_Sample": fasta_sample_id,
            "fasta_Batch": batch_info,
            "Analyzed_Intervals": interval_str,
            "FIS_Variants": format_variants_list(fis_genotypes_all, variant_separator),  # Show ALL FIS variants
            "fasta_Variants": format_variants_list(fasta_genotypes, variant_separator),
            "Concordant": concordant,
            "Concordant_Enhanced": concordant_enhanced,
            "FIS_Unique": format_variants_list(fis_unique, variant_separator),
            "fasta_Unique": format_variants_list(fasta_unique, variant_separator),
            "FIS_Unique_Enhanced": format_variants_list(fis_unique_enhanced, variant_separator),
            "fasta_Unique_Enhanced": format_variants_list(fasta_unique_enhanced, variant_separator),
            "FIS_Count": len(fis_genotypes_all),  # Count ALL FIS variants
            "fasta_Count": len(fasta_genotypes),
            "FIS_Flag": flag_status,
            "FIS_Flag_Reasons": flag_reasons_str,
        }

        comparison_results.append(result)

    # Create DataFrame
    df = pd.DataFrame(comparison_results)

    # Sort by sample name for better organization
    df = df.sort_values("FIS_Sample")

    # Output TSV file
    logger.info(f"Writing comparison results to {output_file}")
    df.to_csv(output_file, sep="\t", index=False)

    # Log summary statistics
    total_samples = len(df)
    concordant_samples = len(df[df["Concordant"] == "Y"])
    discordant_samples = len(df[df["Concordant"] == "N"])
    na_samples = len(df[df["Concordant"] == "N/A"])

    # Enhanced concordance statistics
    concordant_samples_enhanced = len(df[df["Concordant_Enhanced"] == "Y"])
    discordant_samples_enhanced = len(df[df["Concordant_Enhanced"] == "N"])
    na_samples_enhanced = len(df[df["Concordant_Enhanced"] == "N/A"])

    # Count flagging statistics
    flagged_samples = len(df[df["FIS_Flag"] == "Y"])
    unflagged_samples = len(df[df["FIS_Flag"] == "N"])

    # Count mapping statistics
    samples_with_fasta_match = len(df[df["fasta_Sample"] != "Not found"])
    samples_without_fasta_match = len(df[df["fasta_Sample"] == "Not found"])

    logger.info("FIS vs FASTA comparison completed:")
    logger.info(f"  Total FIS samples processed: {total_samples}")
    logger.info(
        f"  Samples with FASTA match: {samples_with_fasta_match} "
        f"({samples_with_fasta_match / total_samples * 100:.1f}%)",
    )
    logger.info(
        f"  Samples without FASTA match: {samples_without_fasta_match} "
        f"({samples_without_fasta_match / total_samples * 100:.1f}%)",
    )
    logger.info("  Concordance analysis (original method):")
    logger.info(f"    - Concordant: {concordant_samples} ({concordant_samples / total_samples * 100:.1f}%)")
    logger.info(f"    - Discordant: {discordant_samples} ({discordant_samples / total_samples * 100:.1f}%)")
    logger.info(f"    - N/A (missing data): {na_samples} ({na_samples / total_samples * 100:.1f}%)")
    logger.info("  Concordance analysis (enhanced method with IUPAC & alignment):")
    logger.info(
        f"    - Concordant: {concordant_samples_enhanced} ({concordant_samples_enhanced / total_samples * 100:.1f}%)",
    )
    logger.info(
        f"    - Discordant: {discordant_samples_enhanced} ({discordant_samples_enhanced / total_samples * 100:.1f}%)",
    )
    logger.info(f"    - N/A (missing data): {na_samples_enhanced} ({na_samples_enhanced / total_samples * 100:.1f}%)")
    logger.info("  FIS variant flagging analysis:")
    logger.info(f"    - Flagged samples: {flagged_samples} ({flagged_samples / total_samples * 100:.1f}%)")
    logger.info(f"    - Unflagged samples: {unflagged_samples} ({unflagged_samples / total_samples * 100:.1f}%)")


def main() -> None:
    """Main function."""
    args = arg_parser()

    try:
        compare_fis_fasta_direct(
            fis_tsv=args.fis_tsv,
            fasta_json=args.fasta_json,
            output_file=args.output_file,
            variant_separator=args.separator,
        )
        logger.success("FIS vs FASTA comparison completed successfully")
    except Exception as e:
        logger.error(f"Error during comparison: {e}")
        raise


if __name__ == "__main__":
    main()
