#!/usr/bin/env python3
"""Calculate frequency of a specific variant from merged statistics JSON.

# Run with defaults (uses bundled default JSON and variant 16037 A G)
python src/modules/NGS/calculate_variant_freq.py

# Run specifying merged JSON and variant (POSITION REF ALT)
python src/modules/NGS/calculate_variant_freq.py merged_statistics.json 225 G A
"""

import json
import sys
from pathlib import Path
from typing import Any

from loguru import logger

from src.core.variants import pos_base

_DEFAULT_JSON = "/mnt/bca/mtDNA/science/results/merged/merged_statistics.json"
_DEFAULT_POSITION = "16037"
_DEFAULT_REF = "A"
_DEFAULT_ALT = "G"
_MIN_ARGS = 2
_ARG_POSITION = 2
_ARG_REF = 3
_ARG_ALT = 4


def calculate_variant_frequency(json_file: str, position: str, ref: str, alt: str) -> dict[str, Any]:
    """
    Calculate the frequency of a specific variant across all samples.

    Args:
        json_file: Path to the merged statistics JSON file
        position: Position of the variant (e.g., "16037")
        ref: Reference allele (e.g., "A")
        alt: Alternate allele (e.g., "G")

    Returns:
        Dictionary with frequency statistics
    """
    variant_key = f"{position} {ref} {alt}"
    pos = pos_base(position)

    with Path(json_file).open() as f:
        data = json.load(f)

    total_samples = len(data)
    samples_with_variant = 0
    sample_list: list[str] = []

    for sample_id, sample_data in data.items():
        # Check if the sample has variant information
        if "variants" in sample_data:
            variants = sample_data["variants"]

            # Check SNPs for the specific variant
            if "snps" in variants:
                for snp in variants["snps"]:
                    if snp.get("pos") == pos and snp.get("ref") == ref and snp.get("seq") == alt:
                        samples_with_variant += 1
                        sample_list.append(sample_id)
                        break  # Found the variant in this sample, move to next sample

    frequency = (samples_with_variant / total_samples * 100) if total_samples > 0 else 0

    return {
        "variant": variant_key,
        "total_samples": total_samples,
        "samples_with_variant": samples_with_variant,
        "frequency_percent": round(frequency, 2),
        "sample_ids": sample_list,
    }


def main() -> None:
    """Main function to calculate variant frequency."""
    json_file = _DEFAULT_JSON if len(sys.argv) < _MIN_ARGS else sys.argv[1]

    # Default variant: 16037 A G
    position = sys.argv[_ARG_POSITION] if len(sys.argv) > _ARG_POSITION else _DEFAULT_POSITION
    ref = sys.argv[_ARG_REF] if len(sys.argv) > _ARG_REF else _DEFAULT_REF
    alt = sys.argv[_ARG_ALT] if len(sys.argv) > _ARG_ALT else _DEFAULT_ALT

    logger.info(f"Analyzing variant: {position} {ref} {alt}")
    logger.info(f"Reading file: {json_file}")
    logger.info("-" * 60)

    results = calculate_variant_frequency(json_file, position, ref, alt)

    logger.info("Results:")
    logger.info(f"  Variant: {results['variant']}")
    logger.info(f"  Total samples: {results['total_samples']}")
    logger.info(f"  Samples with variant: {results['samples_with_variant']}")
    logger.info(f"  Frequency: {results['frequency_percent']}%")

    if results["samples_with_variant"] > 0:
        logger.info("Samples with this variant:")
        for i, sample_id in enumerate(results["sample_ids"], 1):
            logger.info(f"  {i}. {sample_id}")


if __name__ == "__main__":
    main()
