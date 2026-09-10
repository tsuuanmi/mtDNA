#!/usr/bin/env python3
"""
Blind Copy Files Script

This script copies AB1 and FASTA files for samples listed in a blind mapping TSV,
renames them using blind IDs, and organizes them into folders.
It also copies Sequencher TXT variant tables (filter by batch first, like AB1/FASTA)
into a 'sequencher' subfolder, and creates a filtered JSON statistics file from a
predefined merged statistics file.

Usage:
    python -m src.modules.statistics.blind_copy \
        --data-dir DATA_DIR --results-dir RESULTS_DIR \
        [--mapping MAPPING_FILE] [--output OUTPUT_DIR]

Examples:
    python -m src.modules.statistics.blind_copy --data-dir /data --results-dir /results
    python -m src.modules.statistics.blind_copy \
        --data-dir /data --results-dir /results \
        --mapping custom_mapping.tsv --output /path/to/output
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

from loguru import logger


def load_metadata_file(metadata_file: Path) -> dict[str, str]:
    """
    Load the metadata TSV file and return a dictionary of sample_id -> batch_id.

    Args:
        metadata_file: Path to the metadata TSV file

    Returns:
        Dictionary mapping sample IDs to batch IDs
    """
    metadata = {}

    if not metadata_file.exists():
        return metadata

    try:
        with metadata_file.open() as f:
            # Skip header line
            f.readline().strip()

            for line_num, raw_line in enumerate(f, 2):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 2:  # noqa: PLR2004
                    logger.warning(f"Invalid line {line_num} in metadata file: {line}")
                    continue

                batch_id, sample_id = parts[0], parts[1]
                metadata[sample_id] = batch_id

        # Only log if metadata found
        if metadata:
            logger.info(f"Loaded metadata for {len(metadata)} samples")

    except Exception as e:  # noqa: BLE001
        logger.error(f"Error reading metadata file {metadata_file}: {e}")
        logger.exception("Full traceback:")

    return metadata


def read_mapping_file(mapping_file: Path) -> dict[str, str]:
    """
    Read the mapping file and return a dictionary of blind_id -> sample_id.

    Args:
        mapping_file: Path to the TSV mapping file

    Returns:
        Dictionary mapping blind IDs to sample IDs
    """
    mapping = {}

    try:
        with mapping_file.open() as f:
            for line_num, raw_line in enumerate(f, 1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) != 2:  # noqa: PLR2004
                    logger.warning(f"Invalid line {line_num} in mapping file: {line}")
                    continue

                blind_id, sample_id = parts
                mapping[blind_id] = sample_id

    except Exception as e:  # noqa: BLE001
        logger.error(f"Error reading mapping file {mapping_file}: {e}")
        sys.exit(1)

    return mapping


def find_files(directory: Path, pattern: str, *, recursive: bool = True) -> list[Path]:
    """Find files matching a pattern in a directory."""
    if not directory.exists():
        return []

    try:
        files = list(directory.rglob(pattern)) if recursive else list(directory.glob(pattern))
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error searching for files with pattern '{pattern}' in {directory}: {e}")
        return []
    return files


def find_sample_files(  # noqa: C901
    data_dir: Path,
    results_dir: Path,
    sample_id: str,
    sample_batch_id: str | None = None,
    fallback_batch_id: str | None = None,
) -> tuple[list[Path], list[Path], list[Path]]:
    """
    Find AB1, FASTA, and Sequencher TXT files for a given sample ID.

    Sequencher TXT variant tables live under
    ``results_dir/manual_pipeline/manual/<batch>/`` and are matched by batch
    first (like AB1/FASTA), then by an anchored LID prefix so numeric LIDs
    that are substrings of each other (e.g. 2428 vs 24285) do not cross-match.

    Args:
        data_dir: Base data directory
        results_dir: Base results directory
        sample_id: Sample ID to search for
        sample_batch_id: Specific batch ID for this sample from metadata
        fallback_batch_id: Fallback batch ID from command line/env

    Returns:
        Tuple of (ab1_files, fasta_files, txt_files)
    """
    ab1_files = []
    fasta_files = []
    txt_files = []

    # Use sample-specific batch ID first, then fallback batch ID
    batch_id = sample_batch_id or fallback_batch_id

    try:
        # Search for AB1 files in data directory
        if batch_id:
            # Search in specific batch directory first
            batch_raw_dir = data_dir / "raw" / batch_id
            ab1_files.extend(find_files(batch_raw_dir, f"*{sample_id}*.ab1"))

        if not ab1_files and not sample_batch_id:
            # Only search all directories if we don't have sample-specific batch ID
            raw_dir = data_dir / "raw"
            ab1_files.extend(find_files(raw_dir, f"*{sample_id}*.ab1"))

        # Search for FASTA files in results directory
        if batch_id:
            # Search in specific batch directory first
            batch_fasta_dir = results_dir / "fasta" / batch_id
            fasta_files.extend(find_files(batch_fasta_dir, f"*{sample_id}*.fasta"))

        if not fasta_files and not sample_batch_id:
            # Only search all directories if we don't have sample-specific batch ID
            fasta_dir = results_dir / "fasta"
            fasta_files.extend(find_files(fasta_dir, f"*{sample_id}*.fasta"))

        # Fallback: search archived results if not found in active results
        if not fasta_files:
            if batch_id:
                batch_archive_dir = results_dir / "archive" / batch_id
                fasta_files.extend(find_files(batch_archive_dir, f"*{sample_id}*.fasta"))

            if not fasta_files and not sample_batch_id:
                archive_dir = results_dir / "archive"
                fasta_files.extend(find_files(archive_dir, f"*{sample_id}*.fasta"))

        # Search for Sequencher TXT variant tables (filter by batch first).
        # Names are LID-prefixed: '<LID>.TXT' or '<LID>-<ranges>.TXT', so anchor
        # the LID at the start to avoid substring cross-matches.
        if batch_id:
            batch_txt_dir = results_dir / "manual_pipeline" / "manual" / batch_id
            txt_files.extend(find_files(batch_txt_dir, f"{sample_id}.TXT"))
            txt_files.extend(find_files(batch_txt_dir, f"{sample_id}-*.TXT"))
        if not txt_files and not sample_batch_id:
            txt_root = results_dir / "manual_pipeline" / "manual"
            txt_files.extend(find_files(txt_root, f"{sample_id}.TXT"))
            txt_files.extend(find_files(txt_root, f"{sample_id}-*.TXT"))

    except Exception as e:  # noqa: BLE001
        logger.error(f"Error searching for files for sample {sample_id}: {e}")
        logger.exception("Full traceback:")

    return ab1_files, fasta_files, txt_files


def parse_fasta(fasta_file: Path) -> dict[str, str]:
    """
    Parse a FASTA file and extract sequences by HV region.
    Includes HV1, HV2, and HV3.

    Args:
        fasta_file: Path to the FASTA file

    Returns:
        Dictionary mapping HV region names to sequences (HV1, HV2, and HV3)
    """
    sequences = {}
    current_region = None
    current_seq = []

    try:
        with fasta_file.open() as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Save previous sequence if exists and is HV1, HV2, or HV3
                    if current_region and current_seq and current_region in ["HV1", "HV2", "HV3"]:
                        sequences[current_region] = "".join(current_seq)
                        current_seq = []

                    # Extract region name (HV1, HV2, HV3)
                    header = line[1:]  # Remove the '>' character
                    region_match = re.search(r"(HV\d)", header)
                    if region_match:
                        region = region_match.group(1)
                        # Process HV1, HV2, and HV3
                        if region in ["HV1", "HV2", "HV3"]:
                            current_region = region
                            current_seq = []
                        else:
                            current_region = None
                            current_seq = []
                    else:
                        current_region = header.split("|")[0] if "|" in header else header
                        current_seq = []
                elif current_region and current_region in ["HV1", "HV2", "HV3"]:
                    current_seq.append(line)

        # Save the last sequence if it's HV1, HV2, or HV3
        if current_region and current_seq and current_region in ["HV1", "HV2", "HV3"]:
            sequences[current_region] = "".join(current_seq)

        # Reduced logging

    except Exception as e:  # noqa: BLE001
        logger.error(f"Error parsing FASTA file {fasta_file}: {e}")

    return sequences


def split_fasta_by_region(
    fasta_file: Path,
    output_dir: Path,
    sample_id: str,
    *,
    remove_original: bool = True,
) -> tuple[list[Path], list[str]]:
    """
    Split a FASTA file into separate files for each HV region.
    Processes HV1, HV2, and HV3.

    Args:
        fasta_file: Path to the FASTA file
        output_dir: Directory to write the split files
        sample_id: Sample ID to use in the output filenames
        remove_original: Whether to remove the original FASTA file after splitting

    Returns:
        Tuple of (created files list, HV regions list)
    """
    created_files = []
    hv_regions = []

    try:
        # Parse the FASTA file
        sequences = parse_fasta(fasta_file)

        if not sequences:
            logger.warning(f"No HV1/HV2/HV3 sequences found in {fasta_file}")
            return created_files, hv_regions, hv_regions  # type: ignore[reportReturnType]

        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write each region to a separate file
        for region, sequence in sequences.items():
            output_file = output_dir / f"{sample_id}_{region}.fasta"

            with output_file.open("w") as f:
                f.write(f">{sample_id} {region}\n")
                f.write(sequence)

            created_files.append(output_file)
            hv_regions.append(region)
            # Remove verbose logging - only essential info

        # Remove original FASTA file if requested and split was successful
        if remove_original and created_files:
            try:
                fasta_file.unlink()
                # Remove verbose logging
            except Exception as e:  # noqa: BLE001
                logger.warning(f"Could not remove original FASTA file {fasta_file}: {e}")

    except Exception as e:  # noqa: BLE001
        logger.error(f"Error splitting FASTA file {fasta_file}: {e}")

    return created_files, hv_regions


def copy_and_rename_file(src_file: Path, dest_dir: Path, old_pattern: str, new_pattern: str) -> bool:
    """
    Copy a file and rename it by replacing the old pattern with the new pattern.

    Args:
        src_file: Source file path
        dest_dir: Destination directory
        old_pattern: Pattern to replace in filename
        new_pattern: New pattern to use in filename

    Returns:
        True if successful, False otherwise
    """
    try:
        # Create destination directory if it doesn't exist
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Generate new filename
        old_filename = src_file.name
        new_filename = old_filename.replace(old_pattern, new_pattern)

        dest_file = dest_dir / new_filename

        # Copy the file
        shutil.copy2(src_file, dest_file)
        # Remove verbose logging - only log errors
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error copying {src_file}: {e}")
        return False
    return True


def clean_sample_variants(sample_data: dict) -> dict:
    """
    Remove 'file' field from all variants in sample data to maintain complete blinding.

    Args:
        sample_data: Dictionary containing sample data

    Returns:
        Cleaned sample data with 'file' fields removed from variants
    """
    if not isinstance(sample_data, dict):
        return sample_data

    # Create a deep copy to avoid modifying the original
    cleaned_data = sample_data.copy()

    if "variants" in cleaned_data and isinstance(cleaned_data["variants"], dict):
        cleaned_variants = {}

        for variant_type, variant_list in cleaned_data["variants"].items():
            if isinstance(variant_list, list):
                cleaned_variant_list = []
                for variant in variant_list:
                    if isinstance(variant, dict):
                        # Create a copy without the 'file' field
                        cleaned_variant = {k: v for k, v in variant.items() if k != "file"}
                        cleaned_variant_list.append(cleaned_variant)
                    else:
                        cleaned_variant_list.append(variant)
                cleaned_variants[variant_type] = cleaned_variant_list
            else:
                cleaned_variants[variant_type] = variant_list

        cleaned_data["variants"] = cleaned_variants

    return cleaned_data


def process_samples(  # noqa: C901,PLR0912,PLR0915
    mapping: dict[str, str],
    data_dir: Path,
    results_dir: Path,
    output_dir: Path,
    metadata: dict[str, str] | None = None,
    fallback_batch_id: str | None = None,
    *,
    merged_stats_file: Path | None = None,
) -> None:
    """
    Process all samples in the mapping file.

    Args:
        mapping: Dictionary of blind_id -> sample_id
        data_dir: Base data directory
        results_dir: Base results directory
        output_dir: Output directory for blinded files
        metadata: Dictionary of sample_id -> batch_id from metadata file
        fallback_batch_id: Fallback batch ID from command line/env
        merged_stats_file: Path to merged statistics JSON file
    """
    total_samples = len(mapping)
    processed_samples = 0
    total_files_copied = 0
    total_files_split = 0
    total_txt_copied = 0

    # Data completeness tracking
    samples_with_complete_data = 0  # 4 AB1 + 3 FASTA (HV1, HV2, HV3)
    samples_with_sufficient_ab1 = 0  # >= 4 AB1 files
    samples_with_all_hv_regions = 0  # All 3 HV regions
    samples_missing_data = []

    # Detailed tracking
    ab1_counts = []
    fasta_counts = []
    hv_region_stats = {"HV1": 0, "HV2": 0, "HV3": 0}

    # Create centralized JSON directory
    json_output_dir = output_dir / "JSON"
    json_output_dir.mkdir(parents=True, exist_ok=True)

    # Create centralized Sequencher TXT directory (blinded variant tables)
    sequencher_output_dir = output_dir / "sequencher"
    sequencher_output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Processing {total_samples} samples...")

    for i, (blind_id, sample_id) in enumerate(mapping.items(), 1):
        # Concise progress logging
        progress_prefix = f"[{i}/{total_samples}]"

        try:
            # Find files for this sample
            sample_batch_id = metadata.get(sample_id) if metadata else None

            ab1_files, fasta_files, txt_files = find_sample_files(
                data_dir,
                results_dir,
                sample_id,
                sample_batch_id,
                fallback_batch_id,
            )

            # Count files
            num_ab1 = len(ab1_files)
            num_fasta = len(fasta_files)
            num_txt = len(txt_files)
            ab1_counts.append(num_ab1)
            fasta_counts.append(num_fasta)

            # Check AB1 files
            if num_ab1 >= 4:  # noqa: PLR2004
                samples_with_sufficient_ab1 += 1

            if not ab1_files and not fasta_files and not txt_files:
                logger.warning(f"{progress_prefix} {sample_id}: No files found")
                samples_missing_data.append(
                    {
                        "sample_id": sample_id,
                        "blind_id": blind_id,
                        "ab1_count": 0,
                        "fasta_count": 0,
                        "txt_count": 0,
                        "hv_regions": [],
                    },
                )
                continue
        except Exception as e:  # noqa: BLE001
            logger.error(f"Error searching for files for sample {sample_id}: {e}")
            continue

        # Create output directory for this blind ID
        sample_output_dir = output_dir / blind_id

        files_copied = 0
        files_split = 0
        sample_hv_regions = set()

        try:
            # Copy AB1 files
            for ab1_file in ab1_files:
                if copy_and_rename_file(ab1_file, sample_output_dir, sample_id, blind_id):
                    files_copied += 1

            # Split each FASTA file by HV region (splitting is always enabled)
            for fasta_file in fasta_files:
                split_files, hv_regions = split_fasta_by_region(
                    fasta_file,
                    sample_output_dir,
                    blind_id,
                    remove_original=False,
                )
                files_split += len(split_files)
                sample_hv_regions.update(hv_regions)
                if split_files:
                    files_copied += len(split_files)  # Count split files as "copied"

            # Copy Sequencher TXT variant tables, renamed LID -> blind ID.
            for txt_file in txt_files:
                if copy_and_rename_file(txt_file, sequencher_output_dir, sample_id, blind_id):
                    files_copied += 1
                    total_txt_copied += 1

            # Check HV regions and completeness
            if sample_hv_regions:
                hv_region_stats["HV1"] += 1 if "HV1" in sample_hv_regions else 0
                hv_region_stats["HV2"] += 1 if "HV2" in sample_hv_regions else 0
                hv_region_stats["HV3"] += 1 if "HV3" in sample_hv_regions else 0

                if len(sample_hv_regions) == 3:  # noqa: PLR2004
                    samples_with_all_hv_regions += 1

            # Check for complete data and log status
            has_complete_data = num_ab1 >= 4 and len(sample_hv_regions) == 3  # noqa: PLR2004
            if has_complete_data:
                samples_with_complete_data += 1
                logger.success(
                    f"{progress_prefix} {sample_id}: ✓ Complete (AB1:{num_ab1}, HV:{len(sample_hv_regions)})",
                )
            else:
                # Log warning for incomplete data
                status_parts = []
                if num_ab1 < 4:  # noqa: PLR2004
                    status_parts.append(f"AB1:{num_ab1}/4")
                else:
                    status_parts.append(f"AB1:{num_ab1}")

                if len(sample_hv_regions) < 3:  # noqa: PLR2004
                    missing_hv = {"HV1", "HV2", "HV3"} - sample_hv_regions
                    status_parts.append(f"Missing:{','.join(sorted(missing_hv))}")
                else:
                    status_parts.append(f"HV:{len(sample_hv_regions)}")

                logger.warning(f"{progress_prefix} {sample_id}: ⚠ Incomplete ({', '.join(status_parts)})")

                missing_info = {
                    "sample_id": sample_id,
                    "blind_id": blind_id,
                    "ab1_count": num_ab1,
                    "fasta_count": num_fasta,
                    "txt_count": num_txt,
                    "hv_regions": sorted(sample_hv_regions),
                }
                samples_missing_data.append(missing_info)

        except Exception as e:  # noqa: BLE001
            logger.error(f"Error processing files for sample {sample_id}: {e}")
            logger.exception("Full traceback:")
            continue

        if files_copied > 0:
            processed_samples += 1
            total_files_copied += files_copied
            total_files_split += files_split
        else:
            logger.warning(f"{progress_prefix} {sample_id}: No files copied")

    # Concise Summary
    logger.info(f"{'=' * 80}")
    logger.success("PROCESSING COMPLETE")
    logger.info(f"{'=' * 80}")
    logger.info(f"Processed: {processed_samples}/{total_samples} samples | Files: {total_files_copied}")
    logger.info(f"Sequencher TXT copied: {total_txt_copied}")
    logger.info(
        f"Complete data (4 AB1 + 3 HV): {samples_with_complete_data}/{total_samples} "
        f"({(samples_with_complete_data / total_samples) * 100:.1f}%)",
    )
    logger.info(
        f"Sufficient AB1 (≥4): {samples_with_sufficient_ab1}/{total_samples} "
        f"({(samples_with_sufficient_ab1 / total_samples) * 100:.1f}%)",
    )
    logger.info(
        f"All HV regions: {samples_with_all_hv_regions}/{total_samples} "
        f"({(samples_with_all_hv_regions / total_samples) * 100:.1f}%)",
    )

    # HV regions coverage
    logger.info(f"HV Coverage - HV1:{hv_region_stats['HV1']} HV2:{hv_region_stats['HV2']} HV3:{hv_region_stats['HV3']}")

    # Show incomplete samples only if there are any
    if samples_missing_data:
        logger.warning(f"⚠ Incomplete data: {len(samples_missing_data)} samples")
        for _idx, sample_info in enumerate(samples_missing_data[:10], 1):
            parts = []
            if sample_info["ab1_count"] < 4:  # noqa: PLR2004
                parts.append(f"AB1:{sample_info['ab1_count']}/4")
            hv_regions = sample_info.get("hv_regions", [])
            if len(hv_regions) < 3:  # noqa: PLR2004
                missing = ",".join(sorted({"HV1", "HV2", "HV3"} - set(hv_regions)))
                parts.append(f"Missing:{missing}")
            logger.warning(f"  {sample_info['sample_id']}: {' | '.join(parts)}")
        if len(samples_missing_data) > 10:  # noqa: PLR2004
            logger.warning(f"  ... and {len(samples_missing_data) - 10} more")

    logger.info(f"{'=' * 80}")

    # Create individual JSON files and merged JSON file
    if merged_stats_file:
        create_json_files(merged_stats_file, mapping, json_output_dir)


def clean_variants_data(variants_dict: dict) -> dict:
    """
    Remove 'file' field from all variants to maintain complete blinding.

    Args:
        variants_dict: Dictionary containing variants data

    Returns:
        Cleaned variants dictionary with 'file' fields removed
    """
    cleaned = {}

    for variant_type, variant_list in variants_dict.items():
        if isinstance(variant_list, list):
            cleaned_variants = []
            for variant in variant_list:
                if isinstance(variant, dict):
                    # Create a copy without the 'file' field
                    cleaned_variant = {k: v for k, v in variant.items() if k != "file"}
                    cleaned_variants.append(cleaned_variant)
                else:
                    cleaned_variants.append(variant)
            cleaned[variant_type] = cleaned_variants
        else:
            cleaned[variant_type] = variant_list

    return cleaned


def create_json_files(merged_stats_file: Path, mapping: dict[str, str], json_output_dir: Path) -> bool:
    """
    Create individual JSON files for each sample and a merged JSON file in the JSON directory.
    All files use blind IDs and have 'file' fields removed from variants.

    Args:
        merged_stats_file: Path to the merged statistics JSON file
        mapping: Dictionary of blind_id -> sample_id
        json_output_dir: Directory to save the JSON files

    Returns:
        True if successful, False otherwise
    """

    if merged_stats_file and merged_stats_file.exists():
        try:
            with merged_stats_file.open() as f:
                all_stats = json.load(f)
        except (json.JSONDecodeError, Exception) as e:  # noqa: BLE001
            logger.error(f"Error loading merged statistics: {e}")
            all_stats = {}
    else:
        all_stats = {}

    # Create individual JSON files and collect data for merged file
    merged_data = {}
    created_files = 0
    missing_samples = []

    try:
        for blind_id, sample_id in mapping.items():
            if sample_id in all_stats:
                # Get the sample data and clean it
                sample_data = all_stats[sample_id].copy()

                # Clean variants to remove 'file' fields
                if "variants" in sample_data:
                    sample_data["variants"] = clean_variants_data(sample_data["variants"])

                # Create individual JSON file with blind_id
                individual_file = json_output_dir / f"{blind_id}.json"
                with individual_file.open("w") as f:
                    json.dump(sample_data, f, indent=2)

                # Add to merged data using blind_id as key
                merged_data[blind_id] = sample_data
                created_files += 1
            else:
                missing_samples.append(sample_id)

        # Create merged JSON file
        if merged_data:
            merged_file = json_output_dir / "merged.json"
            with merged_file.open("w") as f:
                json.dump(merged_data, f, indent=2)

        logger.info(f"{'=' * 80}")
        logger.success(f"JSON: {created_files} individual + 1 merged ({len(merged_data)}/{len(mapping)} samples)")

        if missing_samples:
            logger.warning(f"⚠ {len(missing_samples)} samples not in merged statistics")

        logger.info(f"{'=' * 80}")
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error creating JSON files: {e}")
        logger.exception("Full traceback:")
        return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Copy and blind AB1 and FASTA files based on mapping file, and create filtered JSON statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m src.modules.statistics.blind_copy --data-dir /data --results-dir /results
  python -m src.modules.statistics.blind_copy \
    --data-dir /data --results-dir /results \
    --mapping custom_mapping.tsv --output /path/to/output
        """,
    )

    parser.add_argument(
        "--mapping",
        type=str,
        help="Path to mapping TSV file (default: data/modules/statistics/sample_mapping.tsv)",
    )
    parser.add_argument("--output", type=str, help="Output directory (default: blinded_output)")
    parser.add_argument(
        "--metadata",
        type=str,
        help="Path to metadata TSV file (default: /results/merged/merged_metadata.tsv)",
    )
    parser.add_argument(
        "--merged-stats",
        type=str,
        help="Path to merged statistics JSON file (default: /results/merged/merged_statistics.json)",
    )
    parser.add_argument("--data-dir", type=str, required=True, help="Data directory path")
    parser.add_argument("--results-dir", type=str, required=True, help="Results directory path")

    args = parser.parse_args()

    # Get script and repository directories
    script_dir = Path(__file__).parent
    project_root = script_dir.parents[2]

    # Set paths
    mapping_file = (
        Path(args.mapping)
        if args.mapping
        else project_root / "data" / "modules" / "statistics" / "sample_mapping.tsv"
    )
    output_dir = Path(args.output) if args.output else script_dir / "blinded_output"
    data_dir = Path(args.data_dir)
    results_dir = Path(args.results_dir)
    merged_dir = results_dir / "merged"
    metadata_file = Path(args.metadata) if args.metadata else merged_dir / "merged_metadata.tsv"
    merged_stats_file = Path(args.merged_stats) if args.merged_stats else merged_dir / "merged_statistics.json"

    # Check if mapping file exists
    if not mapping_file.exists():
        logger.error(f"Mapping file not found: {mapping_file}")
        sys.exit(1)

    # Read mapping file
    mapping = read_mapping_file(mapping_file)
    if not mapping:
        logger.error("No valid mappings found in mapping file")
        sys.exit(1)

    logger.info(f"Processing {len(mapping)} samples from {mapping_file.name}")

    # Load metadata file
    metadata = load_metadata_file(metadata_file)

    logger.info(f"Processing {len(mapping)} samples from {mapping_file.name}")

    # Process samples
    try:
        process_samples(
            mapping,
            data_dir,
            results_dir,
            output_dir,
            metadata,
            None,
            merged_stats_file=merged_stats_file,
        )
        logger.success(f"Success! Blinded files are available in: {output_dir}")
    except KeyboardInterrupt:
        logger.warning("Operation cancelled by user")
        sys.exit(1)
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error: {e}")
        logger.exception("Full traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()
