#!/usr/bin/env python
"""
This script compares variant analysis results between FIS NGS pipeline
(transformed TSV with Sample_ID column) and Sanger pipeline (merged_statistics.json).
It generates TSV and Excel files to facilitate quality control by identifying
concordance between the two analysis methods.

This version reads sample IDs directly from the transformed TSV and uses a compact
MT-to-STR/HID mapping file to resolve run-specific FIS sample names to Sanger IDs.

Variants are compared directly by position and allele. Exact alleles match, and
IUPAC ambiguity codes match compatible alleles at the same position (for example,
199C matches 199Y because Y represents C or T). Uninformative N alleles and the
established always-concordant positions are excluded.

The script shows ALL FIS variants detected by the NGS pipeline, not filtered to
Sanger analyzed intervals. This provides a comprehensive view of FIS detection
capabilities while comparing concordance within Sanger's analyzed regions.

The script extracts the actual analyzed intervals from each Sanger sample
(from the 'intervals' field in merged_statistics.json) and uses these intervals
for concordance comparison, but shows all FIS variants in the output for
complete analysis coverage.

Example usage:
    python compare_fis_sanger_direct.py \
        -f results/FIS_transformed.tsv \
        -s results/merged_statistics.json \
        -t results/modules/NGS/mapping.tsv \
        -o results/FIS_sanger_comparison_direct.tsv

Output columns:
    - FIS_Sample: FIS NGS sample identifier (Sample_ID from TSV)
    - Sample_ID: Base sample ID; omitted when identical to FIS_Sample
    - Sanger_Sample: Corresponding Sanger sample ID (matched by base ID)
    - Sanger_Batch: Batch information from Sanger analysis (e.g., "20241223_mtDNA_12_L1")
    - Analyzed_Intervals: Actual intervals analyzed by Sanger (e.g., "73-340, 16024-16365")
    - FIS_Variants: Nomenclature-corrected FIS variants (not filtered by analyzed intervals)
    - Sanger_Variants: Variants found by Sanger (within analyzed intervals)
    - Concordant: Whether the profiles match using exact or same-position IUPAC-compatible alleles (Y/N)
    - FIS_Unique: Informative variants found only by FIS NGS
    - Sanger_Unique: Informative variants found only by Sanger
    - FIS_QC: Low depth QC flags from FIS analysis
    - FIS_Flag: Combined QC flag (Y if any QC issues from R2 or flagging analysis)
    - FIS_Flag_Level: Flag severity level (0=None, 1=Low, 2=Medium, 3=High, 4=Critical)
    - FIS_Flag_Info: QC flag details from R2 analysis
    - FIS_Flag_Reasons: Specific reasons for flagging (e.g., "Insertion at 16182, 16180-16193 region")
    - Sanger_Flag: Whether Sanger variants contain problematic patterns (Y/N)
    - Sanger_Flag_Reasons: Specific reasons for flagging
"""

import json
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any, NamedTuple

import pandas as pd
from loguru import logger

from src.core.flagging import SampleFlagger
from src.core.variants import (
    format_variants_list,
    is_position_in_intervals,
    iupac_to_bases,
    normalize_position,
    parse_variant_position_allele,
    pos_base,
)
from src.tools.fis.etl import load_reference_bases, profile_to_variants
from src.tools.fis.utils import clean_sample_id, filter_genotypes_by_intervals, parse_genotypes

# Configuration constants
IUPAC_AMBIGUITY_CODES = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
SPECIAL_POSITIONS = ["455", "463", "573", "309"]
ALWAYS_CONCORDANT_POSITIONS = {"309.1", "309.2", "573.1", "16193.1"}
FLAGGING_REGION_START = 16180
FLAGGING_REGION_END = 16193
MAPPING_MT_COLUMN = "mtDNA"
MAPPING_STR_COLUMN = "STR"
FIS_NOMENCLATURE_GENOTYPE_COLUMN = "FIS_Nomenclature correction genotype"


def _variants_for_flagging(genotypes: list[str]) -> list:
    """Convert legacy genotype lists through the canonical FIS adapter."""
    return profile_to_variants(" ".join(genotypes), load_reference_bases(Path("ref/rCRS.fasta")))


def arg_parser() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Compare variants between FIS NGS and Sanger pipelines (direct Sample_ID)")
    parser.add_argument(
        "-f",
        "--fis_tsv",
        type=Path,
        required=True,
        help="Path to transformed TSV file with Sample_ID column (FIS NGS)",
    )
    parser.add_argument(
        "-s",
        "--sanger_json",
        type=Path,
        required=True,
        help="Path to merged_statistics.json file (Sanger)",
    )
    parser.add_argument(
        "-t",
        "--mapping-file",
        "--tsv_mapping",
        dest="mapping_file",
        type=Path,
        help="Path to the MT-to-STR/HID mapping TSV file",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=Path,
        required=True,
        help="Path to output TSV file; an Excel file is written alongside it",
    )
    parser.add_argument(
        "--separator",
        type=str,
        default="space",
        choices=["comma", "space"],
        help="Separator for variant listings: comma or space. Default is space",
    )
    return parser.parse_args()


def parse_sample_mapping(mapping_path: Path | None) -> dict[str, str]:
    """Parse MT and STR/HID aliases into a canonical Sanger mtDNA ID mapping."""
    if mapping_path is None:
        return {}

    mapping_data = pd.read_csv(mapping_path, sep="\t", dtype=str)
    required_columns = {MAPPING_MT_COLUMN, MAPPING_STR_COLUMN}
    missing_columns = required_columns - set(mapping_data.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        msg = f"Mapping file {mapping_path} is missing required columns: {missing}"
        raise ValueError(msg)

    mapping: dict[str, str] = {}
    for mt_value, str_value in mapping_data[[MAPPING_MT_COLUMN, MAPPING_STR_COLUMN]].itertuples(index=False):
        mt_id = clean_sample_id(mt_value)
        str_id = clean_sample_id(str_value)
        if not mt_id or not str_id:
            continue

        for alias in {mt_id, str_id}:
            existing_mt_id = mapping.get(alias)
            if existing_mt_id is not None and existing_mt_id != mt_id:
                msg = f"Mapping alias {alias!r} points to both {existing_mt_id!r} and {mt_id!r}"
                raise ValueError(msg)
            mapping[alias] = mt_id

    return mapping


def get_sanger_intervals(sample_data: dict[str, Any]) -> list[tuple]:
    """
    Extract analyzed intervals from Sanger sample data.

    Args:
        sample_data: Sample data dictionary from Sanger JSON

    Returns:
        List of (start, end) tuples representing analyzed intervals
    """
    intervals = []

    if "intervals" in sample_data:
        interval_data = sample_data["intervals"]

        # Process each region (HV1, HV2, HV3, etc.)
        for ranges in interval_data.values():
            if ranges:  # Check if there are intervals for this region
                intervals.extend((int(start), int(end)) for start, end in ranges)

    return intervals


def get_all_variants_from_sanger_with_intervals(sample_data: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Get all variants (SNPs, insertions, deletions) from Sanger sample data.
    Only includes variants within the actual analyzed intervals for this sample.

    Args:
        sample_data: Sample data dictionary from Sanger JSON

    Returns:
        Combined list of all variants with normalized format, filtered by actual intervals
    """
    all_variants = []

    # Get the actual analyzed intervals for this sample
    intervals = get_sanger_intervals(sample_data)

    if "variants" in sample_data:
        variants = sample_data["variants"]

        # Process SNPs, insertions, and deletions (filtered by analyzed intervals)
        for variant_type, key in [("snp", "snps"), ("ins", "insertions"), ("del", "deletions")]:
            all_variants.extend(
                {
                    "pos": normalize_position(variant["pos"]),
                    "type": variant_type,
                    "ref": variant.get("ref", ""),
                    "seq": variant.get("seq", ""),
                }
                for variant in variants.get(key, [])
                if "pos" in variant and is_position_in_intervals(normalize_position(variant["pos"]), intervals)
            )

    return all_variants


def convert_variants_to_genotypes(variants: list[dict[str, Any]]) -> list[str]:
    """
    Convert variant list to genotype format.

    Args:
        variants: List of variant dictionaries

    Returns:
        List of genotype strings
    """
    genotypes = []

    for variant in sorted(variants, key=lambda x: pos_base(x["pos"])):
        pos = variant["pos"]
        seq = variant.get("seq", "")
        variant_type = variant.get("type", "")

        # Create genotype string: position + sequence
        # For deletions, convert "-" to "DEL" format
        genotype = (f"{pos}DEL" if seq == "-" or variant_type == "del" else f"{pos}{seq}") if seq else pos

        genotypes.append(genotype)

    return genotypes


def _filter_informative_genotypes(genotypes: list[str]) -> list[str]:
    """Remove duplicates, uninformative N alleles, and always-concordant positions."""
    informative_genotypes = []
    for genotype in genotypes:
        position, allele = parse_variant_position_allele(genotype)
        if position not in ALWAYS_CONCORDANT_POSITIONS and not allele.upper().endswith("N"):
            informative_genotypes.append(genotype)
    return list(dict.fromkeys(informative_genotypes))


def is_iupac_concordant(allele_a: str, allele_b: str) -> bool:
    """Check whether two alleles represent at least one shared IUPAC base."""
    if allele_a == allele_b:
        return True

    bases_a = set(iupac_to_bases(allele_a.upper()))
    bases_b = set(iupac_to_bases(allele_b.upper()))
    return bool(bases_a & bases_b)


def _find_unmatched_genotypes(genotypes_a: list[str], genotypes_b: list[str]) -> tuple[list[str], list[str]]:
    """Return variants unmatched by exact or same-position IUPAC-compatible comparison."""
    remaining_a = _filter_informative_genotypes(genotypes_a)
    remaining_b = _filter_informative_genotypes(genotypes_b)

    for variant_a in remaining_a.copy():
        pos_a, allele_a = parse_variant_position_allele(variant_a)
        for variant_b in remaining_b.copy():
            pos_b, allele_b = parse_variant_position_allele(variant_b)
            if pos_a == pos_b and is_iupac_concordant(allele_a, allele_b):
                remaining_a.remove(variant_a)
                remaining_b.remove(variant_b)
                break

    return remaining_a, remaining_b


def find_unique_genotypes(genotypes_a: list[str], genotypes_b: list[str]) -> list[str]:
    """Find informative variants in A without an exact or same-position IUPAC match in B."""
    unmatched_a, _ = _find_unmatched_genotypes(genotypes_a, genotypes_b)
    return unmatched_a


def is_concordant_genotypes(genotypes_a: list[str], genotypes_b: list[str]) -> str:
    """Compare informative variants using exact or same-position IUPAC-compatible alleles."""
    informative_a = _filter_informative_genotypes(genotypes_a)
    informative_b = _filter_informative_genotypes(genotypes_b)
    if not informative_a or not informative_b:
        return "N/A"

    unmatched_a, unmatched_b = _find_unmatched_genotypes(informative_a, informative_b)
    return "Y" if not unmatched_a and not unmatched_b else "N"


def find_sanger_sample_match(
    sample_id: str,
    sanger_samples: set[str],
    sample_mapping: dict[str, str] | None = None,
) -> str | None:
    """Find a Sanger sample by exact ID or the longest mapped FIS ID prefix."""
    clean_id = sample_id.replace("_L02", "").replace("_L01", "").strip().strip("_")
    if clean_id in sanger_samples:
        return clean_id

    if not sample_mapping:
        return None

    matching_aliases = [alias for alias in sample_mapping if clean_id == alias or clean_id.startswith(f"{alias}_")]
    if not matching_aliases:
        return None

    mapped_id = sample_mapping[max(matching_aliases, key=len)]
    return mapped_id if mapped_id in sanger_samples else None


class SangerMatch(NamedTuple):
    """Resolved Sanger match data for a FIS sample."""

    sanger_sample_id: str
    sanger_intervals: list[tuple]
    fis_genotypes_all: list[str]
    fis_genotypes_filtered: list[str]
    sanger_genotypes: list[str]
    batch_info: str


def _resolve_sanger_match(
    sanger_sample_id: str | None,
    sanger_data: dict[str, Any],
    fis_genotypes: list[str],
) -> SangerMatch:
    """Resolve Sanger match data, or return a no-match result with all FIS variants shown."""
    if sanger_sample_id and sanger_sample_id in sanger_data:
        sanger_intervals = get_sanger_intervals(sanger_data[sanger_sample_id])
        fis_genotypes_filtered = filter_genotypes_by_intervals(fis_genotypes, sanger_intervals)
        sanger_variants = get_all_variants_from_sanger_with_intervals(sanger_data[sanger_sample_id])
        sanger_genotypes = convert_variants_to_genotypes(sanger_variants)
        batch_info = sanger_data[sanger_sample_id].get("batch", "Unknown")
        return SangerMatch(
            sanger_sample_id, sanger_intervals, fis_genotypes, fis_genotypes_filtered, sanger_genotypes, batch_info
        )

    # No Sanger match - show all FIS variants
    return SangerMatch("Not found", [], fis_genotypes, [], [], "N/A")


def _compute_fis_flags(fis_row: pd.Series, fis_genotypes_all: list[str]) -> tuple[str, str]:
    """Compute the combined FIS flag and reasons string from precomputed columns or fallback flagging."""
    fis_flag_precomputed = fis_row.get("FIS_Flag", "")
    fis_flag_reasons_precomputed = fis_row.get("FIS_Flag_Reasons", "")
    fis_flag_info = fis_row.get("FIS_Flag_Info", "")

    if str(fis_flag_precomputed) and fis_flag_reasons_precomputed is not None:
        combined_fis_flag = str(fis_flag_precomputed)
        fis_flag_reasons_str = str(fis_flag_reasons_precomputed) if pd.notna(fis_flag_reasons_precomputed) else ""
        return combined_fis_flag, fis_flag_reasons_str

    # Fallback: compute flagging (for backward compatibility)
    fis_variants_for_flagging = _variants_for_flagging(fis_genotypes_all)
    _, flag_reasons_fis, _ = SampleFlagger(fis_variants_for_flagging).analyze()
    fis_flag_reasons_str = ", ".join(flag_reasons_fis) if flag_reasons_fis else ""

    # Combine FIS_Flag: Y if any QC issues from R2 or flagging analysis
    fis_flag_info_str = str(fis_flag_info).strip() if pd.notna(fis_flag_info) else ""  # pyright: ignore[reportGeneralTypeIssues]
    has_r2_flag = bool(fis_flag_info_str)
    has_flag_reasons = bool(flag_reasons_fis)
    combined_fis_flag = "Y" if (has_r2_flag or has_flag_reasons) else "N"
    return combined_fis_flag, fis_flag_reasons_str


def _format_interval_string(sanger_sample_id: str, sanger_data: dict[str, Any], sanger_intervals: list[tuple]) -> str:
    """Format analyzed intervals for display."""
    if sanger_sample_id != "Not found" and sanger_sample_id in sanger_data:
        interval_ranges = [f"{start}-{end}" for start, end in sanger_intervals]
        return ", ".join(interval_ranges) if interval_ranges else "None"
    return "N/A"


def _build_comparison_result(
    fis_sample_id: object,
    base_sample_id: object,
    match: SangerMatch,
    interval_str: str,
    concordant: str,
    fis_unique: list[str],
    sanger_unique: list[str],
    fis_qc_flag: object,
    combined_fis_flag: str,
    fis_flag_level: object,
    fis_flag_info: object,
    fis_flag_reasons_str: str,
    sanger_flag_status: str,
    sanger_flag_reasons_str: str,
    variant_separator: str,
) -> dict[str, Any]:
    """Build the comparison result row dictionary."""
    return {
        "FIS_Sample": fis_sample_id,
        "Sample_ID": base_sample_id,
        "Sanger_Sample": match.sanger_sample_id,
        "Sanger_Batch": match.batch_info,
        "Analyzed_Intervals": interval_str,
        "FIS_Variants": format_variants_list(match.fis_genotypes_all, variant_separator),
        "Sanger_Variants": format_variants_list(match.sanger_genotypes, variant_separator),
        "Concordant": concordant,
        "FIS_Unique": format_variants_list(fis_unique, variant_separator),
        "Sanger_Unique": format_variants_list(sanger_unique, variant_separator),
        "FIS_QC": fis_qc_flag,
        "FIS_Flag": combined_fis_flag,
        "FIS_Flag_Level": fis_flag_level,
        "FIS_Flag_Info": fis_flag_info,
        "FIS_Flag_Reasons": fis_flag_reasons_str,
        "Sanger_Flag": sanger_flag_status,
        "Sanger_Flag_Reasons": sanger_flag_reasons_str,
    }


def _get_sort_key(row: pd.Series) -> str:
    """Sort matched Sanger samples first and unmatched samples last."""
    sanger_sample = str(row["Sanger_Sample"])
    return sanger_sample if sanger_sample and sanger_sample != "Not found" else "zzz_not_found"


def _drop_redundant_sample_id_column(df_results: pd.DataFrame) -> pd.DataFrame:
    """Omit Sample_ID when it contains the same identifiers as FIS_Sample."""
    required_columns = {"FIS_Sample", "Sample_ID"}
    if not required_columns.issubset(df_results.columns):
        return df_results

    fis_sample_ids = df_results["FIS_Sample"].astype("string")
    sample_ids = df_results["Sample_ID"].astype("string")
    if not fis_sample_ids.equals(sample_ids):
        return df_results

    logger.info("FIS_Sample and Sample_ID are identical; omitting Sample_ID from outputs")
    return df_results.drop(columns="Sample_ID")


def _log_comparison_summary(df_results: pd.DataFrame, output_file: Path) -> None:
    """Log concise summary statistics for the comparison."""
    total_samples = len(df_results)
    matched_samples = len(df_results[df_results["Sanger_Sample"] != "Not found"])
    concordant_samples = len(df_results[df_results["Concordant"] == "Y"])
    fis_flagged_samples = len(df_results[df_results["FIS_Flag"] == "Y"])
    sanger_flagged_samples = len(df_results[df_results["Sanger_Flag"] == "Y"])

    matched_pct = matched_samples / total_samples * 100
    concordant_pct = concordant_samples / matched_samples * 100
    logger.success(
        f"Saved {output_file.name}: {total_samples} FIS, {matched_samples} matched ({matched_pct:.0f}%)",
    )
    logger.info(f"Concordance (exact/IUPAC): {concordant_samples} ({concordant_pct:.0f}%)")
    logger.info(f"Flagged: FIS={fis_flagged_samples}, Sanger={sanger_flagged_samples}")


def compare_fis_sanger_direct(
    fis_tsv: Path,
    sanger_json: Path,
    output_file: Path,
    mapping_file: Path | None = None,
    variant_separator: str = "space",
) -> None:
    """
    Compare variants between FIS NGS and Sanger pipelines using direct Sample_ID matching.

    Args:
        fis_tsv: Path to transformed TSV file with Sample_ID column (FIS NGS)
        sanger_json: Path to merged_statistics.json (Sanger)
        output_file: Path to output TSV file; Excel uses the same stem with an .xlsx suffix
        mapping_file: Optional path to the MT-to-STR/HID mapping TSV file
        variant_separator: Separator for variant display
    """
    logger.info(f"Comparing FIS ({fis_tsv.name}) vs Sanger ({sanger_json.name})")

    sample_mapping = parse_sample_mapping(mapping_file)

    # Load FIS TSV file
    fis_data = pd.read_csv(fis_tsv, sep="\t")

    # Load Sanger JSON file
    with sanger_json.open() as f:
        sanger_data = json.load(f)

    # Get sample lists - support both old (Sample_ID) and new (FIS_Sample) column names
    sample_col = "FIS_Sample" if "FIS_Sample" in fis_data.columns else "Sample_ID"
    sanger_samples = set(sanger_data)

    # Prepare comparison results
    comparison_results = []

    # Process each sample in FIS data
    for _, fis_row in fis_data.iterrows():
        fis_sample_id = fis_row[sample_col]
        # Compare the nomenclature-corrected profile; retain old inputs as compatibility fallbacks.
        fis_genotypes_str = fis_row.get(
            FIS_NOMENCLATURE_GENOTYPE_COLUMN,
            fis_row.get("FIS_Variants", fis_row.get("Genotype", "")),
        )
        fis_qc_flag = fis_row.get("FIS_QC", "")  # Get FIS_QC column for depth flags
        fis_flag_info = fis_row.get("FIS_Flag_Info", "")  # QC flag details
        base_sample_id = fis_row.get("FIS_sample_id_base", fis_row.get("sample_id_base", fis_row.get("Sample_ID", "")))
        fis_flag_level = fis_row.get("FIS_Flag_Level", "")  # Flag level (0-4)

        # Fallback to computing if columns not present (backward compatibility)
        if not base_sample_id:
            base_sample_id = str(fis_sample_id)

        # Parse FIS genotypes
        fis_genotypes = parse_genotypes(fis_genotypes_str)

        # Find matching Sanger sample
        sanger_sample_id = find_sanger_sample_match(str(base_sample_id), sanger_samples, sample_mapping)
        match = _resolve_sanger_match(sanger_sample_id, sanger_data, fis_genotypes)

        concordant = is_concordant_genotypes(match.fis_genotypes_filtered, match.sanger_genotypes)
        fis_unique = find_unique_genotypes(match.fis_genotypes_filtered, match.sanger_genotypes)
        sanger_unique = find_unique_genotypes(match.sanger_genotypes, match.fis_genotypes_filtered)

        # Flag analysis for FIS variants
        combined_fis_flag, fis_flag_reasons_str = _compute_fis_flags(fis_row, match.fis_genotypes_all)

        # Flag analysis for Sanger variants
        sanger_variants_for_flagging = _variants_for_flagging(match.sanger_genotypes)
        should_flag_sanger, flag_reasons_sanger, _ = SampleFlagger(sanger_variants_for_flagging).analyze()
        sanger_flag_status = "Y" if should_flag_sanger else "N"
        sanger_flag_reasons_str = ", ".join(flag_reasons_sanger) if flag_reasons_sanger else ""

        # Format intervals for display
        interval_str = _format_interval_string(match.sanger_sample_id, sanger_data, match.sanger_intervals)

        # Create comparison result
        result = _build_comparison_result(
            fis_sample_id,
            base_sample_id,
            match,
            interval_str,
            concordant,
            fis_unique,
            sanger_unique,
            fis_qc_flag,
            combined_fis_flag,
            fis_flag_level,
            fis_flag_info,
            fis_flag_reasons_str,
            sanger_flag_status,
            sanger_flag_reasons_str,
            variant_separator,
        )

        comparison_results.append(result)

    # Create DataFrame and save results
    df_results = pd.DataFrame(comparison_results)
    df_results["sort_key"] = df_results.apply(_get_sort_key, axis=1)
    df_results = df_results.sort_values("sort_key").drop("sort_key", axis=1).reset_index(drop=True)
    df_results = _drop_redundant_sample_id_column(df_results)

    # Save TSV and Excel files with the same stem.
    excel_file = output_file.with_suffix(".xlsx")
    df_results.to_csv(output_file, sep="\t", index=False)
    df_results.to_excel(excel_file, index=False)
    logger.success(f"Saved comparison results to {output_file} and {excel_file}")

    # Log summary statistics
    _log_comparison_summary(df_results, output_file)


def main() -> None:
    """Main function."""
    args = arg_parser()

    try:
        compare_fis_sanger_direct(
            fis_tsv=args.fis_tsv,
            sanger_json=args.sanger_json,
            output_file=args.output_file,
            mapping_file=args.mapping_file,
            variant_separator=args.separator,
        )
    except Exception as e:
        logger.error(f"Error during comparison: {e}")
        raise


if __name__ == "__main__":
    main()
