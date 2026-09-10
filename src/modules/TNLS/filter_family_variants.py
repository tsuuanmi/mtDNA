#!/usr/bin/env python
"""
Filter mtDNA family variant data based on specific criteria:
1. Filter out samples with "No ID" in variant
2. Keep only Family_ID with 2 or more people in family
3. Filter all families where people in family have the same variant profile
4. Output samples sorted by Family_ID where people in family have different variant profiles

This script uses loguru for structured logging.
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from loguru import logger

from src.core.variants import iupac_to_bases, normalize_position

_MIN_FAMILY_SIZE = 2


def parse_variant(variant_str: str) -> tuple:
    """
    Parse a variant string to extract position and base(s).

    Args:
        variant_str: Variant string like "513A", "16092Y", "249", "315.1C", etc.

    Returns:
        Tuple of (position, base) or (None, None) if parsing fails
    """
    variant_str = variant_str.strip()
    if not variant_str:
        return None, None

    # Handle different variant formats:
    # 1. Standard format: position + base (e.g., "513A", "16092Y")
    match = re.match(r"^(\d+)([A-Za-z]+)$", variant_str)
    if match:
        position = int(match.group(1))
        base = match.group(2).upper()
        return position, base

    # 2. Position with decimal + base (e.g., "315.1C")
    match = re.match(r"^(\d+\.\d+)([A-Za-z]+)$", variant_str)
    if match:
        position = normalize_position(match.group(1))
        base = match.group(2).upper()
        return position, base

    # 3. Position only (e.g., "249", "523", "524")
    match = re.match(r"^(\d+)$", variant_str)
    if match:
        position = int(match.group(1))
        base = None  # No specific base, indicates deletion or insertion
        return position, base

    # 4. Position with decimal only (e.g., "315.1")
    match = re.match(r"^(\d+\.\d+)$", variant_str)
    if match:
        position = normalize_position(match.group(1))
        base = None
        return position, base

    return None, None


def normalize_variant_with_iupac(variant_str: str) -> str:
    """
    Normalize a variant considering IUPAC codes.
    Converts hetero/homo variants at same position to a canonical form.

    Args:
        variant_str: Variant string like "513A", "513R", "249", etc.

    Returns:
        Normalized variant string
    """
    position, base = parse_variant(variant_str)
    if position is None:
        return variant_str  # Return original if can't parse

    # If no base (e.g., "249"), return as-is
    if base is None:
        return variant_str

    # Get all possible bases for this IUPAC code
    possible_bases = iupac_to_bases(base)

    # Create a canonical representation using sorted bases
    canonical_bases = "".join(sorted(set(possible_bases)))

    # Map common combinations to their IUPAC code for consistency
    base_to_iupac = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "AC": "M",
        "AG": "R",
        "AT": "W",
        "CG": "S",
        "CT": "Y",
        "GT": "K",
        "ACG": "V",
        "ACT": "H",
        "AGT": "D",
        "CGT": "B",
        "ACGT": "N",
    }

    canonical_iupac = base_to_iupac.get(canonical_bases, canonical_bases)
    return f"{position}{canonical_iupac}"


def variants_are_equivalent(variant1: str, variant2: str) -> bool:
    """
    Check if two variants are equivalent considering IUPAC codes.

    Args:
        variant1: First variant string
        variant2: Second variant string

    Returns:
        True if variants are at same position and have overlapping bases
    """
    pos1, base1 = parse_variant(variant1)
    pos2, base2 = parse_variant(variant2)

    if pos1 != pos2 or pos1 is None:
        return False

    # If both variants have no base (e.g., "249" vs "249"), they are equivalent
    if base1 is None and base2 is None:
        return True

    # If one has base and other doesn't, they are not equivalent
    if base1 is None or base2 is None:
        return False

    # Get possible bases for each variant and check for overlap
    bases1 = set(iupac_to_bases(base1))
    bases2 = set(iupac_to_bases(base2))

    # Variants are equivalent if they have any overlapping bases at the same position
    return bool(bases1.intersection(bases2))


def normalize_variant_list_with_iupac(variant_list: list[str]) -> set[str]:
    """
    Normalize a list of variants considering IUPAC equivalence.
    Converts each variant to its canonical IUPAC form.

    Args:
        variant_list: List of variant strings

    Returns:
        Set of normalized variant strings
    """
    if not variant_list:
        return set()

    normalized_variants = set()
    unparseable = []

    for variant in variant_list:
        pos, _base = parse_variant(variant)
        if pos is not None:
            # Normalize each variant individually to its canonical IUPAC form
            normalized = normalize_variant_with_iupac(variant)
            normalized_variants.add(normalized)
        else:
            unparseable.append(variant)

    # Add unparseable variants as-is
    normalized_variants.update(unparseable)

    return normalized_variants


def parse_family_ids(family_id_str: str) -> list[str]:
    """Parse Family_ID column which can contain multiple comma-separated family IDs"""
    if pd.isna(family_id_str) or family_id_str == "":
        return []
    return [fam_id.strip() for fam_id in str(family_id_str).split(",")]


def create_family_signature(family_variant_strings: dict[str, str]) -> str:
    """
    Create a family-wide signature that combines all variants from all family members.
    This allows families where members differ only in hetero/homo variants at the same
    positions to be considered equivalent.
    """
    if not family_variant_strings:
        return ""

    # Collect all variants across the family
    all_position_variants = defaultdict(set)
    all_position_no_base = set()  # For variants like "249", "523"
    all_unparseable = set()

    for variant_str in family_variant_strings:
        variant_list = parse_variants_to_list(variant_str)

        for variant in variant_list:
            pos, _base = parse_variant(variant)
            if pos is not None:
                if _base is not None:
                    # Add all possible bases for this IUPAC code to the position
                    possible_bases = iupac_to_bases(_base)
                    all_position_variants[pos].update(possible_bases)
                else:
                    # Position without base (e.g., "249", "523")
                    all_position_no_base.add(pos)
            else:
                all_unparseable.add(variant)

    # Create canonical representation for the entire family
    canonical_variants = []

    # Handle positions with bases
    for pos in sorted(all_position_variants.keys()):
        bases = all_position_variants[pos]
        canonical_bases = "".join(sorted(bases))

        # Map to IUPAC code for consistency
        base_to_iupac = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "AC": "M",
            "AG": "R",
            "AT": "W",
            "CG": "S",
            "CT": "Y",
            "GT": "K",
            "ACG": "V",
            "ACT": "H",
            "AGT": "D",
            "CGT": "B",
            "ACGT": "N",
        }
        canonical_iupac = base_to_iupac.get(canonical_bases, canonical_bases)
        canonical_variants.append(f"{pos}{canonical_iupac}")

    # Handle positions without bases
    canonical_variants.extend(str(pos) for pos in sorted(all_position_no_base))

    # Add unparseable variants
    canonical_variants.extend(sorted(all_unparseable))

    return " ".join(canonical_variants)


def normalize_variants(variant_str: str) -> str:
    """Normalize variant string by replacing commas with spaces, preserving original order"""
    if pd.isna(variant_str) or variant_str == "":
        return ""

    # Simply replace commas with spaces and clean up whitespace
    # This preserves the original order from the input file
    if "," in str(variant_str):
        # Replace commas with spaces and clean up extra whitespace
        normalized = " ".join([v.strip() for v in str(variant_str).split(",") if v.strip()])
    else:
        # Already space-separated, just clean up whitespace
        normalized = " ".join([v.strip() for v in str(variant_str).split() if v.strip()])

    return normalized


def normalize_variants_for_comparison(variant_str: str) -> str:
    """Normalize variant string for comparison using IUPAC-aware normalization"""
    if pd.isna(variant_str) or variant_str == "":
        return ""

    # Parse variants and normalize with IUPAC consideration
    variant_list = parse_variants_to_list(variant_str)
    normalized_set = normalize_variant_list_with_iupac(variant_list)

    # Return sorted, space-separated string for comparison purposes
    return " ".join(sorted(normalized_set))


def parse_variants_to_list(variant_str: str) -> list[str]:
    """Parse variant string into a list of individual variants"""
    if pd.isna(variant_str) or variant_str == "":
        return []
    # Handle both comma and space separated variants
    if "," in str(variant_str):
        return [v.strip() for v in str(variant_str).split(",") if v.strip()]
    return [v.strip() for v in str(variant_str).split() if v.strip()]


def find_common_variants_in_family(family_variants_list: list[str]) -> set[str]:  # noqa: C901
    """Find variants that are common to all members in a family using IUPAC-aware comparison"""
    if not family_variants_list:
        return set()

    # Convert each variant string to a normalized set of variants
    variant_sets = []
    for variant_str in family_variants_list:
        variant_list = parse_variants_to_list(variant_str)
        normalized_variants = normalize_variant_list_with_iupac(variant_list)
        variant_sets.append(normalized_variants)

    if not variant_sets:
        return set()

    # Find variants that are equivalent across all family members
    # A variant is "common" if there's an equivalent variant in all family members
    common_variants = set()

    # Start with variants from the first family member
    for variant1 in variant_sets[0]:
        pos1, base1 = parse_variant(variant1)
        if pos1 is None:
            continue

        # Check if this variant has an equivalent in all other family members
        is_common = True
        canonical_variant = variant1  # Will be updated to most inclusive IUPAC code

        for other_variant_set in variant_sets[1:]:
            # Look for equivalent variant at the same position
            found_equivalent = False
            for variant2 in other_variant_set:
                if variants_are_equivalent(variant1, variant2):
                    found_equivalent = True
                    # Update canonical variant to be more inclusive
                    pos2, base2 = parse_variant(variant2)
                    if pos1 == pos2:
                        # Combine bases from both variants
                        bases1 = set(iupac_to_bases(base1))
                        bases2 = set(iupac_to_bases(base2))
                        combined_bases = bases1.union(bases2)
                        canonical_bases = "".join(sorted(combined_bases))

                        # Map to IUPAC code
                        base_to_iupac = {
                            "A": "A",
                            "C": "C",
                            "G": "G",
                            "T": "T",
                            "AC": "M",
                            "AG": "R",
                            "AT": "W",
                            "CG": "S",
                            "CT": "Y",
                            "GT": "K",
                            "ACG": "V",
                            "ACT": "H",
                            "AGT": "D",
                            "CGT": "B",
                            "ACGT": "N",
                        }
                        canonical_iupac = base_to_iupac.get(canonical_bases, canonical_bases)
                        canonical_variant = f"{pos1}{canonical_iupac}"
                    break

            if not found_equivalent:
                is_common = False
                break

        if is_common:
            common_variants.add(canonical_variant)

    return common_variants


def find_truly_unique_variants_among_family(sample_variants: str, all_family_variants: list[str]) -> list[str]:
    """
    Find variants that are truly unique to a sample (not present in any other family member).

    Args:
        sample_variants: Variant string for the current sample
        all_family_variants: List of variant strings for all family members (including current sample)

    Returns:
        List of variants that are unique to this sample only
    """
    if pd.isna(sample_variants) or sample_variants == "":
        return []

    # Get normalized variants for the current sample
    sample_variant_list = parse_variants_to_list(sample_variants)
    normalize_variant_list_with_iupac(sample_variant_list)

    # Get all variants from other family members (excluding the current sample)
    other_family_variants = set()
    for other_variant_str in all_family_variants:
        if other_variant_str != sample_variants:  # Skip variants identical to current sample
            other_variant_list = parse_variants_to_list(other_variant_str)
            other_normalized_variants = normalize_variant_list_with_iupac(other_variant_list)
            other_family_variants.update(other_normalized_variants)

    # Find variants in current sample that are not equivalent to any variant in other family members
    # Preserve the original order from the sample
    truly_unique_variants = []
    for sample_variant in sample_variant_list:
        # Normalize the current variant for comparison
        normalized_sample_variant = normalize_variant_with_iupac(sample_variant)

        is_unique = True
        for other_variant in other_family_variants:
            if variants_are_equivalent(normalized_sample_variant, other_variant):
                is_unique = False
                break
        if is_unique:
            # Add the original variant (not normalized) to preserve original format
            truly_unique_variants.append(sample_variant)

    return truly_unique_variants


def find_truly_unique_variants_for_sample(
    sample_variants: str,
    all_family_variants: list[str],
    sample_index: int,
) -> list[str]:
    """
    Find variants that are truly unique to a sample (not present in any other family member).

    Args:
        sample_variants: Variant string for the current sample
        all_family_variants: List of variant strings for all family members
        sample_index: Index of the current sample in the family list

    Returns:
        List of variants that are unique to this sample only
    """
    if pd.isna(sample_variants) or sample_variants == "":
        return []

    # Get normalized variants for the current sample
    sample_variant_list = parse_variants_to_list(sample_variants)
    sample_normalized_variants = normalize_variant_list_with_iupac(sample_variant_list)

    # Get normalized variants for all other family members
    other_family_variants = set()
    for i, other_variant_str in enumerate(all_family_variants):
        if i != sample_index:  # Skip the current sample
            other_variant_list = parse_variants_to_list(other_variant_str)
            other_normalized_variants = normalize_variant_list_with_iupac(other_variant_list)

            # Check for IUPAC equivalence with each other family member's variants
            for sample_variant in sample_normalized_variants.copy():
                for other_variant in other_normalized_variants:
                    if variants_are_equivalent(sample_variant, other_variant):
                        other_family_variants.add(sample_variant)
                        break

    # Find variants in sample that are not equivalent to any variant in other family members
    truly_unique_variants = sample_normalized_variants - other_family_variants

    return sorted(truly_unique_variants)


def find_unique_variants_for_sample(sample_variants: str, common_variants: list[str]) -> list[str]:
    """Find variants unique to a sample (not in common variants) using IUPAC-aware comparison"""
    if pd.isna(sample_variants) or sample_variants == "":
        return []

    sample_variant_list = parse_variants_to_list(sample_variants)

    # Find variants in sample that are not equivalent to any common variant
    # Preserve the original order from the sample
    unique_variants = []

    for sample_variant in sample_variant_list:
        # Normalize the current variant for comparison
        normalized_sample_variant = normalize_variant_with_iupac(sample_variant)

        is_unique = True
        # Check if this sample variant is equivalent to any common variant
        for common_variant in common_variants:
            if variants_are_equivalent(normalized_sample_variant, common_variant):
                is_unique = False
                break

        if is_unique:
            # Add the original variant (not normalized) to preserve original format
            unique_variants.append(sample_variant)

    return unique_variants


def filter_family_variants(input_file: str, output_file: str) -> None:  # noqa: C901,PLR0912,PLR0915
    """Main filtering function"""
    logger.info("Starting mtDNA family variant filtering")
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output file: {output_file}")

    # Read the TSV file, skipping comment lines that start with //
    try:
        df = pd.read_csv(input_file, sep="\t", comment="#", low_memory=False)
        # If the file starts with a comment line, we need to handle it differently
        with Path(input_file).open() as f:
            first_line = f.readline().strip()
            if first_line.startswith("//"):
                # Skip the comment line and read again
                df = pd.read_csv(input_file, sep="\t", skiprows=1, low_memory=False)
    except Exception as e:
        logger.error(f"Error reading file {input_file}: {e}")
        raise

    logger.info(f"Successfully loaded data: {len(df)} rows, {len(df.columns)} columns")

    # Step 1: Filter out samples with "No ID" in variants
    initial_count = len(df)
    df_filtered = df[df["Variants"] != "No ID"].copy()
    no_id_filtered = len(df_filtered)
    logger.info(
        f"Step 1: Filtered out 'No ID' variants: {initial_count - no_id_filtered} rows removed, "
        f"{no_id_filtered} rows remaining",
    )

    # Step 2: Expand rows for multiple Family_IDs and count family sizes
    family_member_count = defaultdict(set)
    family_variants = defaultdict(set)
    family_variant_strings = defaultdict(list)  # Store actual variant strings for each family

    # Track which rows belong to which families
    row_families = {}

    logger.info("Step 2: Processing family relationships and variant profiles")

    for idx, row in df_filtered.iterrows():
        family_ids = parse_family_ids(row["Family_ID"])  # type: ignore[reportArgumentType]
        row_families[idx] = family_ids

        for family_id in family_ids:
            # Use a combination of identifiers to uniquely identify each person
            person_id = f"{row['Name']}_{row['mtDNA_barcode']}_{idx}"
            family_member_count[family_id].add(person_id)
            # We'll compute family signatures later, after collecting all family members
            family_variant_strings[family_id].append(row["Variants"])  # Store original variant string

    logger.info(f"Found {len(family_member_count)} unique families")

    # Create family signatures using IUPAC-aware combination
    logger.info("Step 2.5: Creating family signatures with IUPAC normalization")
    family_signatures = {}
    for family_id, variant_strings in family_variant_strings.items():
        signature = create_family_signature(variant_strings)  # type: ignore[reportArgumentType]
        family_signatures[family_id] = signature

    # Group families by their signatures
    signature_to_families = defaultdict(list)
    for family_id, signature in family_signatures.items():
        signature_to_families[signature].append(family_id)

    # Step 3: Filter families with 2+ members
    families_with_multiple_members = {
        fam_id for fam_id, members in family_member_count.items() if len(members) >= _MIN_FAMILY_SIZE
    }

    single_member_families = len(family_member_count) - len(families_with_multiple_members)
    logger.info(
        f"Step 3: Families with 2+ members: {len(families_with_multiple_members)} "
        f"(excluded {single_member_families} single-member families)",
    )

    # Step 4: For families with 2+ members, check if they have truly different variant profiles
    # after IUPAC normalization. A family has "different profiles" if family members
    # have variants at different positions (not just hetero vs homo at same position)
    families_with_different_variants = set()
    families_with_same_variants = set()

    for fam_id in families_with_multiple_members:
        family_variant_list = family_variant_strings[fam_id]

        # Check if all family members have variants at the same positions
        all_positions = set()
        member_positions = []

        for variant_str in family_variant_list:
            variant_list = parse_variants_to_list(variant_str)
            positions = set()
            for variant in variant_list:
                pos, _base = parse_variant(variant)
                if pos is not None:
                    positions.add(pos)
            member_positions.append(positions)
            all_positions.update(positions)

        # Check if any family member has variants at different positions
        has_different_positions = False
        if len(member_positions) > 1:
            first_positions = member_positions[0]
            for other_positions in member_positions[1:]:
                if first_positions != other_positions:
                    has_different_positions = True
                    break

        if has_different_positions:
            families_with_different_variants.add(fam_id)
        else:
            families_with_same_variants.add(fam_id)

    same_variant_families = len(families_with_same_variants)
    logger.info(
        f"Step 4: Families with different variant profiles: {len(families_with_different_variants)} "
        f"(excluded {same_variant_families} families with identical variant positions after IUPAC normalization)",
    )

    # Step 5: Filter rows that belong to families with different variants
    logger.info("Step 5: Filtering samples from families with variant differences")
    filtered_rows = []

    for idx, _row in df_filtered.iterrows():
        family_ids = row_families[idx]
        # Keep row if it belongs to at least one family with different variants
        if any(fam_id in families_with_different_variants for fam_id in family_ids):
            filtered_rows.append(idx)

    df_final = df_filtered.loc[filtered_rows].copy()

    excluded_samples = len(df_filtered) - len(df_final)
    logger.info(f"Filtered data: {len(df_final)} samples remaining ({excluded_samples} samples excluded)")

    # Step 5.5: Calculate truly unique variants for each sample
    logger.info("Step 5.5: Calculating truly unique variants for each sample")

    # Create a mapping of sample to family variants for better tracking
    sample_to_family_variants = {}
    for idx, row in df_final.iterrows():
        family_ids = row_families[idx]
        if family_ids:
            primary_family_id = family_ids[0]
            if primary_family_id in family_variant_strings:
                sample_to_family_variants[idx] = {
                    "family_id": primary_family_id,
                    "family_variants": family_variant_strings[primary_family_id],
                    "sample_variants": row["Variants"],
                }

    # Add unique variants column
    unique_variants_list = []
    for idx, _row in df_final.iterrows():
        if idx in sample_to_family_variants:
            family_data = sample_to_family_variants[idx]
            sample_variants = family_data["sample_variants"]
            family_variants = family_data["family_variants"]

            # Find truly unique variants by comparing with all other family members
            unique_variants = find_truly_unique_variants_among_family(sample_variants, family_variants)
            unique_variants_str = " ".join(unique_variants) if unique_variants else "None"
        else:
            unique_variants_str = "None"

        unique_variants_list.append(unique_variants_str)

    # Add the unique variants column to the dataframe
    df_final["Unique_Variants"] = unique_variants_list

    # Step 5.6: Normalize the Variants column to use space separation consistently
    df_final["Variants"] = df_final["Variants"].apply(
        lambda x: normalize_variants(x) if pd.notna(x) and x != "No ID" else x,
    )

    # Step 6: Sort by Family_ID first, then by batch
    logger.info("Step 6: Sorting results by Family_ID, then by batch")

    def get_first_family_id(family_id_str: str) -> str:
        family_ids = parse_family_ids(family_id_str)
        return family_ids[0] if family_ids else "ZZZ"  # Put empty family IDs at the end

    df_final["_sort_family"] = df_final["Family_ID"].apply(get_first_family_id)
    df_final = df_final.sort_values(["_sort_family", "mtDNA_batch", "Name"])
    df_final = df_final.drop("_sort_family", axis=1)

    # Save the filtered data
    logger.info(f"Saving filtered data to {output_file}")
    try:
        df_final.to_csv(output_file, sep="\t", index=False)
        logger.success(f"Successfully saved {len(df_final)} rows to {output_file}")
    except Exception as e:
        logger.error(f"Error saving file {output_file}: {e}")
        raise

    # Print summary statistics
    logger.info("=" * 50)
    logger.info("FILTERING SUMMARY")
    logger.info("=" * 50)
    logger.info(f"Input rows: {len(df)}")
    logger.info(f"After removing 'No ID': {len(df_filtered)}")
    logger.info(f"Final output rows: {len(df_final)}")
    logger.info(
        f"Reduction: {len(df) - len(df_final)} rows ({((len(df) - len(df_final)) / len(df) * 100):.1f}%)"
    )
    logger.info("=" * 50)


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter mtDNA family variant data")
    parser.add_argument("input_file", help="Input TSV file path")
    parser.add_argument("output_file", help="Output TSV file path")

    args = parser.parse_args()

    try:
        filter_family_variants(args.input_file, args.output_file)
        logger.success("Filtering completed successfully!")
        logger.info(f"Output saved to: {args.output_file}")
    except (FileNotFoundError, pd.errors.ParserError, OSError, ValueError, KeyError) as e:
        logger.error(f"Error during filtering: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
