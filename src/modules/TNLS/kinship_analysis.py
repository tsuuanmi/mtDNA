#!/usr/bin/env python
"""
Kinship Analysis Module.

This module provides functions for analyzing kinship relationships between two individuals
based on their STR (Short Tandem Repeat) profiles. It calculates the likelihood of different
types of relationships and determines the most probable relationship.
"""

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

# Constants
KINSHIP_MATRIX = {
    "Identical": [0, 0, 1],
    "Parent-child": [0, 1, 0],
    "Full siblings": [1 / 4, 1 / 2, 1 / 4],
    "Half siblings,Uncle-nephew,Grandparent-grandchild": [1 / 2, 1 / 2, 0],
    "First cousins": [3 / 4, 1 / 4, 0],
    "Unrelated": [1, 0, 0],
    "Second cousin": [7 / 8, 1 / 8, 0],  # Added based on mapping_dict in get_key_with_highest_value
}

DEFAULT_ALLELE_FREQUENCY = 0.1
DEFAULT_THRESHOLD = 10**6
DEFAULT_THRESHOLD_2 = 5 * 10**5


def get_match(p1: list[float], p2: list[float]) -> tuple[list, list]:
    """
    Find matching and non-matching elements between two lists.

    Args:
        p1: List of alleles from first person
        p2: List of alleles from second person

    Returns:
        Tuple containing:
            - List of matching elements
            - List of non-matching elements
    """
    matches = p1 if p1 == p2 else list(set(p1) & set(p2))

    set1 = set(p1)
    set2 = set(p2)
    unmatched_elements = set1.symmetric_difference(set2)
    unmatches = list(unmatched_elements)

    return matches, unmatches


def get_probability(df: pd.DataFrame, marker: str, allele: float) -> float:
    """
    Get the probability of a specific allele at a marker location.

    Args:
        df: DataFrame containing allele frequencies
        marker: STR marker name
        allele: Allele value to look up

    Returns:
        Frequency of the allele or default value if not found
    """
    try:
        if marker in df.columns:
            df_af = df[df["Allele"] == float(allele)]
            af = df_af.iloc[0][marker] if not df_af.empty else DEFAULT_ALLELE_FREQUENCY
        else:
            af = DEFAULT_ALLELE_FREQUENCY
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error getting probability for {marker}:{allele}: {e}")
        return DEFAULT_ALLELE_FREQUENCY
    return 1 if af == 0 else af


def classifier_case(list1: list[float], list2: list[float]) -> tuple[int, dict[str, float]]:  # noqa: C901,PLR0911,RET503  # type: ignore[reportReturnType]
    """
    Classify the type of allele matching pattern between two individuals.

    Args:
        list1: List of alleles from first person
        list2: List of alleles from second person

    Returns:
        Tuple containing:
            - Case number (0-8)
            - Dictionary of relevant alleles for likelihood calculations
    """
    match_list, unmatch_list = get_match(list1, list2)
    match = len(match_list)

    # Case ii,ii (homozygous identical)
    if len(list1) > len(set(list1)) and len(list2) > len(set(list2)):
        if match == 2:  # noqa: PLR2004
            return 0, {"i": list1[0]}  # 'case ii ii'
        if match == 0:
            return 2, {"i": list1[0], "j": list2[0]}  # 'case ii jj'

    # Case ii,ij (first homozygous, second heterozygous)
    elif len(list1) > len(set(list1)) and list2[0] != list2[1]:
        if match == 1:
            return 1, {"i": match_list[0], "j": unmatch_list[0]}  # 'case ii ij'
        if match == 0:
            return 3, {"i": list1[0], "j": list2[0], "k": list2[1]}  # 'case ii jk'

    # Case ij,* (first heterozygous)
    elif list1[0] != list1[1]:
        if match == 2:  # noqa: PLR2004
            return 4, {"i": list1[0], "j": list1[1]}  # 'case ij ij'
        if match == 1:
            if list2[0] == list2[1]:
                return 7, {"i": match_list[0], "j": unmatch_list[0]}  # 'case ij ii or ii ij'
            return 5, {"i": match_list[0], "j": unmatch_list[0], "k": unmatch_list[0]}  # 'case ij ik
        if match == 0:
            if list2[0] == list2[1]:
                return 8, {"i": list2[0], "j": list1[0], "k": list1[1]}  # 'case ij kk'
            return 6, {"i": list1[0], "j": list1[1], "k": list2[0], "l": list2[1]}  # 'case ij kl'


def get_key_with_highest_value(  # noqa: RET503
    dictionary: dict[int, float],
    threshold: float = DEFAULT_THRESHOLD,
    threshold_2: float = DEFAULT_THRESHOLD_2,
) -> tuple[str, float]:  # type: ignore[reportReturnType]
    """
    Determine the highest value in a dictionary and return its key.

    Used to find the most likely relationship based on likelihood ratios.

    Args:
        dictionary: Dictionary of relationship indices and their likelihood values
        threshold: Primary threshold for determining relationship
        threshold_2: Secondary threshold for full siblings

    Returns:
        Tuple containing:
            - String describing the relationship
            - Likelihood ratio value
    """
    if not dictionary:
        return None, None  # type: ignore[reportReturnType]

    logger.debug(f"Input dictionary: {dictionary}")

    mapping_dict = {
        0: "Identical",
        1: "Parent-child",
        2: "Full siblings",
        3: "Half siblings,Uncle-nephew,Grandparent-grandchild",
        4: "Second cousin",
    }

    max_value = max(dictionary.values())
    result_key = None

    for key, value in dictionary.items():
        if value == max_value:
            result_key = key
            break

    if result_key is not None:
        if max_value >= threshold:
            return mapping_dict[result_key], max_value
        if max_value <= threshold and result_key == 2 and max_value >= threshold_2:  # noqa: PLR2004
            return mapping_dict[result_key], max_value
        return "Unrelated", max_value


def calc_likelihood(case_no: int, test_type: int, param_dict: dict[str, float]) -> float:  # noqa: PLR0911,RET503  # type: ignore[reportReturnType]
    """
    Calculate likelihood based on the case number, test type, and parameters.

    Args:
        case_no: Case number (0-8) from classifier_case
        test_type: Index of the relationship type being tested
        param_dict: Dictionary of allele frequencies

    Returns:
        Calculated likelihood value
    """
    index_rel = list(KINSHIP_MATRIX.keys())[test_type]
    k_rel = KINSHIP_MATRIX[index_rel]
    k0, k1, k2 = k_rel[0], 2 * k_rel[1], k_rel[2]

    if case_no == 0:  # case ii ii
        pi = param_dict["pi"]
        return k0 + 2 * k1 / pi + k2 / pi**2
    if case_no == 1:  # case ii ij
        pi, pj = param_dict["pi"], param_dict["pj"]
        return k0 + k1 / pi
    if case_no == 2:  # case ii jj  # noqa: PLR2004
        pi, pj, _pk = param_dict["pi"], param_dict["pj"], param_dict.get("pk", 0)
        return k0
    if case_no == 3:  # case ii jk  # noqa: PLR2004
        pi, pj, _pk = param_dict["pi"], param_dict["pj"], param_dict["pk"]
        return k0
    if case_no == 4:  # case ij ij  # noqa: PLR2004
        pi, pj, _pk = param_dict["pi"], param_dict["pj"], param_dict.get("pk", 0)
        return k0 + k1 * (pi + pj) / (2 * pi * pj) + k2 / (2 * pi * pj)
    if case_no == 5:  # case ij ik  # noqa: PLR2004
        pi, pj, _pk = param_dict["pi"], param_dict["pj"], param_dict["pk"]
        return k0 + k1 / (2 * pi)
    if case_no == 6:  # case ij kl  # noqa: PLR2004
        pi, pj, _pk, _pl = param_dict["pi"], param_dict["pj"], param_dict["pk"], param_dict["pl"]
        return k0
    if case_no == 7:  # case ij ii or ii ij  # noqa: PLR2004
        pi, pj = param_dict["pj"], param_dict["pi"]
        return k0 + k1 / pi
    if case_no == 8:  # case ij kk  # noqa: PLR2004
        _pk, pj, pi = param_dict["pi"], param_dict["pj"], param_dict["pk"]
        return k0


def calc_probability_rela(lr: float) -> float:
    """
    Calculate relationship probability from likelihood ratio.

    Args:
        lr: Likelihood ratio

    Returns:
        Probability as a percentage
    """
    pr = 0.5
    return (lr * pr * 100) / (lr * pr + (1 - pr))


def check_relation(input1: dict[str, dict], input2: dict[str, dict], str_af_file_path: str) -> tuple[str, float, float]:
    """
    Determine the most likely relationship between two individuals.

    Args:
        input1: STR profile for first individual
        input2: STR profile for second individual
        str_af_file_path: Path to the STR allele frequency file

    Returns:
        Tuple containing:
            - String describing the closest relationship
            - Log10 of the likelihood ratio
            - Probability of the relationship
    """
    try:
        str_af_vn = pd.read_excel(str_af_file_path)
        logger.info(f"Successfully loaded allele frequency data from {str_af_file_path}")
    except Exception as e:
        logger.error(f"Failed to load allele frequency data from {str_af_file_path}: {e}")
        raise

    shared_keys = set(input1.keys()) & set(input2.keys())
    logger.info(f"Found {len(shared_keys)} shared STR markers")

    input1 = {key: input1[key] for key in shared_keys}
    input2 = {key: input2[key] for key in shared_keys}

    likelihood_list = {}

    for test_rel in range(5):
        likelihood = 1
        for (key1, _), (key2, _) in zip(input1.items(), input2.items(), strict=False):
            p_dict = {"pi": 0.0, "pj": 0.0, "pk": 0.0, "pl": 0.0}

            a1_a2_dict = list(input1[key1].values())
            p1 = [a1_a2_dict[0], a1_a2_dict[1]]

            a1_a2_dict = list(input2[key2].values())
            p2 = [a1_a2_dict[0], a1_a2_dict[1]]

            case, allele_dict = classifier_case(p1, p2)

            for index_a, i in enumerate(allele_dict.keys()):
                allele_p = get_probability(str_af_vn, key1, allele_dict[i])
                p_dict[list(p_dict.keys())[index_a]] = allele_p

            likelihood *= calc_likelihood(case, test_rel, p_dict)

        likelihood_list[test_rel] = likelihood

    closest_relationship, lr = get_key_with_highest_value(likelihood_list, DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_2)
    logger.info(f"Closest relationship: {closest_relationship} (likelihood: {lr:.2e})")

    if closest_relationship == "Unrelated":
        probability_rela = -1
    else:
        probability_rela = calc_probability_rela(lr)
        logger.info(f"Most likely relationship: {closest_relationship} (probability: {probability_rela:.2f}%)")

    return closest_relationship, math.log10(lr), probability_rela


def load_samples_from_json(json_path: str, lid1: str, lid2: str) -> tuple[dict, dict]:
    """
    Load sample data from merged JSON file by LID.

    Args:
        json_path: Path to the merged JSON file
        lid1: LID of the first sample to compare
        lid2: LID of the second sample to compare

    Returns:
        Tuple containing:
            - First sample data in the required format
            - Second sample data in the required format
    """
    try:
        with Path(json_path).open() as f:
            data = json.load(f)

        logger.info(f"Successfully loaded merged JSON data with {len(data)} samples")

        sample1 = None
        sample2 = None

        for sample in data:
            if sample.get("LID") == lid1:
                sample1 = sample
            elif sample.get("LID") == lid2:
                sample2 = sample

            if sample1 and sample2:
                break

        if not sample1:
            msg = f"Sample with LID {lid1} not found in the dataset"
            raise ValueError(msg)  # noqa: TRY301
        if not sample2:
            msg = f"Sample with LID {lid2} not found in the dataset"
            raise ValueError(msg)  # noqa: TRY301

        logger.info(
            f"Found samples: {lid1} (SampleCode: {sample1['SampleCode']}) and "
            f"{lid2} (SampleCode: {sample2['SampleCode']})",
        )

        formatted_sample1 = format_sample_for_analysis(sample1)
        formatted_sample2 = format_sample_for_analysis(sample2)

    except Exception as e:
        logger.error(f"Error loading samples from JSON: {e}")
        raise

    return formatted_sample1, formatted_sample2


def format_sample_for_analysis(sample: dict) -> dict:
    """
    Convert sample from JSON format to the format required by the check_relation function.

    Args:
        sample: A sample dictionary from the merged JSON

    Returns:
        Formatted dictionary suitable for kinship analysis
    """
    formatted_data = {}

    skip_fields = {
        "SampleCode",
        "LID",
        "batch",
        "barcode",
        "Penta_D",
        "SE33",
        "DYS391",
        "DYS576",
        "DYS570",
        "Yindel",
        "AMEL",
    }

    for key, value in sample.items():
        if key in skip_fields or value is None:
            continue

        if isinstance(value, str) and "," in value:
            alleles = value.split(",")
            if len(alleles) == 2:  # noqa: PLR2004
                try:
                    a1 = float(alleles[0])
                    a2 = float(alleles[1])
                    formatted_data[key] = {"a1": a1, "a2": a2}
                except ValueError:
                    logger.warning(f"Could not convert allele values to float for {key}: {value}")

    return formatted_data


def compare_one_with_all(
    lid1: str,
    json_path: str,
    str_af_file_path: str,
    min_likelihood: float = 0,
) -> list[dict[str, Any]]:
    """
    Compare one sample with all other samples in the dataset.

    Args:
        lid1: LID of the sample to compare with all others
        json_path: Path to the merged JSON file
        str_af_file_path: Path to the STR allele frequency file
        min_likelihood: Minimum log10 likelihood ratio to include in results

    Returns:
        List of dictionaries containing comparison results, sorted by likelihood ratio
    """
    try:
        with Path(json_path).open() as f:
            data = json.load(f)

        logger.info(f"Successfully loaded merged JSON data with {len(data)} samples")

        sample1 = None
        all_samples = []

        for sample in data:
            if sample.get("LID") == lid1:
                sample1 = sample
            else:
                all_samples.append(sample)

        if not sample1:
            msg = f"Sample with LID {lid1} not found in the dataset"
            raise ValueError(msg)  # noqa: TRY301

        logger.info(f"Found target sample: {lid1} (SampleCode: {sample1['SampleCode']})")
        logger.info(f"Will compare with {len(all_samples)} other samples")

        formatted_sample1 = format_sample_for_analysis(sample1)

        results = []
        for i, sample2 in enumerate(all_samples):
            if i % 100 == 0 and i > 0:
                logger.info(f"Processed {i}/{len(all_samples)} comparisons")

            try:
                formatted_sample2 = format_sample_for_analysis(sample2)
                if not formatted_sample2:
                    continue

                relationship, log10_likelihood_ratio, probability_rela = check_relation(
                    formatted_sample1,
                    formatted_sample2,
                    str_af_file_path,
                )

                if log10_likelihood_ratio >= min_likelihood and relationship != "Unrelated":
                    results.append(
                        {
                            "lid": sample2.get("LID"),
                            "sample_code": sample2.get("SampleCode"),
                            "relationship": relationship,
                            "log10_lr": log10_likelihood_ratio,
                            "probability": probability_rela,
                        },
                    )
            except Exception as e:  # noqa: BLE001
                logger.warning(f"Error comparing with sample {sample2.get('LID')}: {e}")
                continue

        results.sort(key=lambda x: x["log10_lr"], reverse=True)

    except Exception as e:
        logger.error(f"Error in compare_one_with_all: {e}")
        raise

    return results


def clean_id(id_value: object) -> str:
    """
    Clean ID values by removing .0 suffix from numeric IDs.

    Args:
        id_value: The ID value to clean

    Returns:
        Cleaned ID as string with .0 suffix removed if present
    """
    id_str = str(id_value)
    if id_str.endswith(".0"):
        return id_str[:-2]
    return id_str


def parse_tsv_mapping(tsv_path: str) -> dict[str, str]:
    """
    Parse the TSV file to create a mapping from mtDNA ID to STR ID.

    Args:
        tsv_path: Path to the TSV file containing the mapping

    Returns:
        A dictionary mapping mtDNA IDs to STR IDs
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t")

        # Create a mapping dictionary from mtDNA_ID to STR_ID
        mtdna_to_str_mapping = {}
        cleaned_count = 0

        for _, row in df.iterrows():
            if pd.notna(row.get("mtDNA_ID")) and pd.notna(row.get("STR_ID")):  # type: ignore[reportGeneralTypeIssues]
                # Clean IDs by removing .0 suffix if present
                mtdna_id = clean_id(row["mtDNA_ID"])
                str_id = clean_id(row["STR_ID"])

                # Track if we cleaned any IDs for logging
                if str(row["mtDNA_ID"]).endswith(".0") or str(row["STR_ID"]).endswith(".0"):
                    cleaned_count += 1

                mtdna_to_str_mapping[mtdna_id] = str_id

        logger.info(
            f"Loaded mapping for {len(mtdna_to_str_mapping)} mtDNA IDs (cleaned {cleaned_count} IDs with .0 suffix)",
        )
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error parsing TSV file: {e}")
        return {}

    return mtdna_to_str_mapping


def get_str_ids_from_mtdna_ids(mtdna_id1: str, mtdna_id2: str, mapping: dict[str, str]) -> tuple[str, str] | None:
    """
    Get STR IDs corresponding to the given mtDNA IDs.

    Args:
        mtdna_id1: First mtDNA ID
        mtdna_id2: Second mtDNA ID
        mapping: Dictionary mapping mtDNA IDs to STR IDs

    Returns:
        Tuple of (str_id1, str_id2) if both mappings exist, otherwise None
    """
    str_id1 = mapping.get(mtdna_id1)
    str_id2 = mapping.get(mtdna_id2)

    if str_id1 is None:
        logger.warning(f"No STR data found for mtDNA ID: {mtdna_id1}")

    if str_id2 is None:
        logger.warning(f"No STR data found for mtDNA ID: {mtdna_id2}")

    if str_id1 is None or str_id2 is None:
        return None

    return (str_id1, str_id2)


def main() -> int:  # noqa: C901,PLR0911,PLR0912,PLR0915
    """
    Main function to parse arguments and run kinship analysis.
    """
    parser = argparse.ArgumentParser(description="Kinship analysis between two samples")
    parser.add_argument("lid1", help="LID of the first sample")
    parser.add_argument("lid2", help='LID of the second sample or "all" to compare with all samples')
    parser.add_argument("--json", required=True, help="Path to merged JSON file")
    parser.add_argument("--str-af", help="Path to STR allele frequency file)")
    parser.add_argument(
        "--min-likelihood",
        type=float,
        default=3.0,
        help="Minimum log10 likelihood ratio to include in results when comparing with all (default: 3.0)",
    )
    parser.add_argument("--top", type=int, default=0, help="Limit results to top N matches (0 = no limit, default: 0)")
    parser.add_argument(
        "--mode",
        choices=["str", "mtdna"],
        default="str",
        help="Mode of operation: 'str' for STR IDs, 'mtdna' for mtDNA IDs",
    )
    parser.add_argument(
        "--mapping-file",
        help="Path to the TSV file with STR to mtDNA ID mapping",
    )

    args = parser.parse_args()

    # Store original IDs for logging purposes
    original_id1 = args.lid1
    original_id2 = args.lid2

    # If mode is mtDNA, convert mtDNA IDs to STR IDs
    if args.mode == "mtdna":
        logger.info("Operating in mtDNA mode, converting mtDNA IDs to STR IDs")
        mapping = parse_tsv_mapping(args.mapping_file)

        if not mapping:
            logger.error("Failed to load mtDNA to STR ID mapping")
            return 1

        # Check if compare with all samples
        if args.lid2.lower() == "all":
            str_id = mapping.get(args.lid1)
            if str_id is None:
                logger.error(f"No STR data found for mtDNA ID: {args.lid1}")
                return 1

            logger.info(f"Mapped mtDNA ID {args.lid1} to STR ID {str_id}")
            args.lid1 = str_id
            # args.lid2 remains "all"
        else:
            str_ids = get_str_ids_from_mtdna_ids(args.lid1, args.lid2, mapping)

            if str_ids is None:
                logger.error("Cannot proceed with kinship analysis: STR data not available for one or both mtDNA IDs")
                return 1

            logger.info(f"Mapped mtDNA IDs {args.lid1} and {args.lid2} to STR IDs {str_ids[0]} and {str_ids[1]}")
            args.lid1, args.lid2 = str_ids
    else:
        logger.info("Operating in STR mode, using provided STR IDs")

    # Continue with existing kinship analysis using the updated args.lid1 and args.lid2
    logger.info(f"Performing kinship analysis for{'':1}{'mtDNA' if args.mode == 'mtdna' else 'STR'} IDs:")
    if args.mode == "mtdna":
        logger.info(f"  Original mtDNA ID 1: {original_id1} → STR ID: {args.lid1}")
        if original_id2.lower() != "all":
            logger.info(f"  Original mtDNA ID 2: {original_id2} → STR ID: {args.lid2}")
    else:
        logger.info(f"  STR ID 1: {args.lid1}")
        if args.lid2.lower() != "all":
            logger.info(f"  STR ID 2: {args.lid2}")

    try:
        if args.lid2.lower() == "all":
            logger.info(f"Comparing sample {args.lid1} with all other samples")
            results = compare_one_with_all(args.lid1, args.json, args.str_af, args.min_likelihood)

            if not results:
                logger.info(f"No relationships found for {args.lid1} above the minimum threshold")
                return 0

            if args.top > 0 and len(results) > args.top:
                results = results[: args.top]

            logger.info(f"Found {len(results)} potential relationships for {args.lid1}:")
            logger.info(f"{'LID':<10} {'SampleCode':<15} {'Relationship':<40} {'Log10 LR':<10} {'Probability':<10}")
            logger.info("-" * 90)

            for result in results:
                logger.info(
                    f"{result['lid']:<10} {result['sample_code']:<15} {result['relationship']:<40} "
                    f"{result['log10_lr']:.2f}{'':>3} {result['probability']:.2f}%",
                )

            return 0
        input1, input2 = load_samples_from_json(args.json, args.lid1, args.lid2)
        if not input1 or not input2:
            logger.error("Could not load valid sample data")
            return 1

        relationship, log10_likelihood_ratio, probability_rela = check_relation(input1, input2, args.str_af)

        logger.info("Kinship Analysis Results:")
        if args.mode == "mtdna":
            logger.info(f"Sample 1 mtDNA ID: {original_id1} (STR ID: {args.lid1})")
            logger.info(f"Sample 2 mtDNA ID: {original_id2} (STR ID: {args.lid2})")
        else:
            logger.info(f"Sample 1 STR ID: {args.lid1}")
            logger.info(f"Sample 2 STR ID: {args.lid2}")
        logger.info(f"Relationship: {relationship}")
        logger.info(f"Log10 Likelihood Ratio: {log10_likelihood_ratio:.2f}")
        if probability_rela > 0:
            logger.info(f"Probability: {probability_rela:.2f}%")
        else:
            logger.info("Probability: Not applicable (unrelated)")

    except Exception as e:  # noqa: BLE001
        logger.exception(f"Error in kinship analysis: {e}")
        return 1

    return 0


if __name__ == "__main__":
    logger.info("Starting kinship analysis")
    sys.exit(main())
