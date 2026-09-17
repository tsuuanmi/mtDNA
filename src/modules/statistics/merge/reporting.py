"""Merged TSV/report writers and TSV-vs-JSON discrepancy analysis."""

import csv
from collections import Counter
from pathlib import Path

from loguru import logger

from src.modules.statistics.merge.config import MIN_TSV_PARTS
from src.modules.statistics.merge.models import BatchSummary, DirectoryScanResult, MergePaths


def create_merged_tsv(
    successful_samples: dict[str, list[str]],
    metadata: dict[str, dict[str, str]],
    paths: MergePaths,
) -> tuple[dict[str, int], int]:
    """
    Create a merged TSV file with batch, sample_id, and barcode for all successful samples.

    Args:
        successful_samples: Dictionary of successful sample IDs per batch
        metadata: Metadata dictionary by batch
        paths: Merge input/output path contract

    Returns:
        Per-batch successful sample counts and the total number of samples written
    """
    samples_per_batch: dict[str, int] = {}
    total_samples = 0

    paths.merged_metadata_tsv.parent.mkdir(parents=True, exist_ok=True)
    with paths.merged_metadata_tsv.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["batch", "sample_id", "barcode"])

        for batch_name, samples in sorted(successful_samples.items()):
            samples_per_batch[batch_name] = 0

            if batch_name in metadata:
                batch_metadata = metadata[batch_name]

                for sample_id in sorted(samples):
                    if sample_id in batch_metadata:
                        barcode = batch_metadata[sample_id]
                        writer.writerow([batch_name, sample_id, barcode])
                        samples_per_batch[batch_name] += 1
                        total_samples += 1
    return samples_per_batch, total_samples


def batch_statistics(
    samples_per_batch: dict[str, int],
    successful_samples: dict[str, list[str]],
    metadata: dict[str, dict[str, str]],
    scan_result: DirectoryScanResult,
    paths: MergePaths,
) -> None:
    """
    Print statistics about the number of samples in each batch with separate columns for regenerate and FASTA.
    Also saves the statistics to a TSV file.

    Args:
        samples_per_batch: Count of successful samples per batch
        successful_samples: Dictionary of successful samples by batch
        metadata: Metadata dictionary by batch
        scan_result: Directory scan results
        paths: Merge input/output path contract
    """
    logger.info("Batch Statistics:")
    logger.info("-" * 150)
    header = f"{'Batch Name':<25} {'Successful':<20} {'Metadata':<20} {'Regenerate':<20} {'Fasta':<20} {'Status':<20}"
    logger.info(header)
    logger.info("-" * 150)

    # Get all unique batch names
    all_batches = set(successful_samples.keys()) | set(metadata.keys()) | set(scan_result.regenerate_samples.keys())
    sorted_batches = sorted(all_batches)

    # Create batch summaries
    batch_summaries = []
    total_summary = BatchSummary("Total")

    for batch_name in sorted_batches:
        summary = BatchSummary(
            batch_name=batch_name,
            successful_count=samples_per_batch.get(batch_name, 0),
            metadata_count=len(metadata.get(batch_name, {})),
            regenerate_count=len(scan_result.regenerate_samples.get(batch_name, set())),
            fasta_count=len(scan_result.fasta_samples.get(batch_name, set())),
        )

        batch_summaries.append(summary)

        # Update totals
        total_summary.successful_count += summary.successful_count
        total_summary.metadata_count += summary.metadata_count
        total_summary.regenerate_count += summary.regenerate_count
        total_summary.fasta_count += summary.fasta_count
        # Log batch statistics
        logger.info(
            f"{summary.batch_name:<25} {summary.successful_count:<20} {summary.metadata_count:<20} "
            f"{summary.regenerate_count:<20} {summary.fasta_count:<20} {summary.status:<20}",
        )

    logger.info("-" * 150)
    logger.info(
        f"{total_summary.batch_name:<25} {total_summary.successful_count:<20} {total_summary.metadata_count:<20} "
        f"{total_summary.regenerate_count:<20} {total_summary.fasta_count:<20}",
    )
    logger.info("-" * 150)

    # Save batch statistics to TSV file
    output_batch_stats = paths.batch_statistics_tsv
    output_batch_stats.parent.mkdir(parents=True, exist_ok=True)

    with output_batch_stats.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Batch", "Successful", "Metadata", "Regenerate", "Fasta", "Status"])

        for summary in batch_summaries:
            writer.writerow(
                [
                    summary.batch_name,
                    summary.successful_count,
                    summary.metadata_count,
                    summary.regenerate_count,
                    summary.fasta_count,
                    summary.status,
                ],
            )

        # Write total row
        writer.writerow(
            [
                total_summary.batch_name,
                total_summary.successful_count,
                total_summary.metadata_count,
                total_summary.regenerate_count,
                total_summary.fasta_count,
                "",
            ],
        )

    logger.success(f"Batch statistics saved to: {output_batch_stats}")


def _discrepancy_flags(tsv_rows: int, *, in_tsv: bool, in_json: bool) -> list[str]:
    """Return every discrepancy flag that applies to one sample."""
    flags: list[str] = []
    if tsv_rows > 1:
        flags.append("duplicate_rows")
    if in_tsv and not in_json:
        flags.append("in_tsv_not_json")
    if in_json and not in_tsv:
        flags.append("in_json_not_tsv")
    return flags or ["matched"]


def _write_sample_discrepancy_tsv(
    output_file: Path,
    merged_statistics: dict[str, dict],
    tsv_samples: set[str],
    sample_counts: Counter[str],
    sample_batches: dict[str, set[str]],
) -> None:
    """
    Save the per-sample discrepancy table: one row per unique sample ID across
    the merged TSV and JSON, with every applicable flag so duplicates and
    one-sided sample IDs are trackable as data instead of console logs.
    """
    json_samples = set(merged_statistics)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "discrepancies", "tsv_rows", "tsv_batches", "json_batch"])
        for sample_id in sorted(tsv_samples | json_samples):
            flags = _discrepancy_flags(
                sample_counts[sample_id],
                in_tsv=sample_id in tsv_samples,
                in_json=sample_id in json_samples,
            )
            json_batch = str(merged_statistics.get(sample_id, {}).get("batch", ""))
            writer.writerow(
                [
                    sample_id,
                    ";".join(flags),
                    sample_counts[sample_id],
                    ",".join(sorted(sample_batches.get(sample_id, set()))),
                    json_batch,
                ],
            )


def analyze_sample_discrepancies(merged_statistics: dict[str, dict], total_samples: int, paths: MergePaths) -> None:
    """
    Analyze and report discrepancies between TSV and JSON data.

    TSV counts include rows and may contain the same sample ID in multiple batches,
    while merged JSON uses unique sample IDs as keys. Compare unique IDs in both
    directions instead of treating the row-count difference as a missing-sample
    count. Also saves a per-sample discrepancy TSV covering every unique sample
    ID across both outputs.

    Args:
        merged_statistics: Merged JSON statistics data
        total_samples: Number of rows in TSV file
        paths: Merge input/output path contract
    """
    json_sample_count = len(merged_statistics)
    tsv_sample_count = total_samples

    merged_metadata_tsv = paths.merged_metadata_tsv
    logger.info(f"Merged JSON file created: {paths.merged_statistics_json}.")
    logger.info(f"Total unique samples in merged JSON: {json_sample_count}")
    logger.info(f"Merged TSV file created: {merged_metadata_tsv}. Total rows in merged TSV: {tsv_sample_count}")

    # Read the TSV once so row counts, duplicate rows, and unique-ID differences
    # are all derived from the same data.
    tsv_samples: set[str] = set()
    sample_batches: dict[str, set[str]] = {}
    sample_counts: Counter[str] = Counter()
    try:
        with merged_metadata_tsv.open(newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)  # Skip header
            for parts in reader:
                if len(parts) >= MIN_TSV_PARTS and parts[1]:
                    batch_id, sample_id = parts[0], parts[1]
                    tsv_samples.add(sample_id)
                    sample_counts[sample_id] += 1
                    sample_batches.setdefault(sample_id, set()).add(batch_id)
    except FileNotFoundError:
        logger.error(f"TSV file {merged_metadata_tsv} not found")
        return

    duplicate_rows = sum(count - 1 for count in sample_counts.values() if count > 1)
    duplicate_ids = sum(count > 1 for count in sample_counts.values())
    if duplicate_rows:
        logger.info(
            f"TSV contains {duplicate_rows} duplicate rows across {duplicate_ids} sample IDs; "
            "merged JSON keeps one key per sample ID.",
        )

    json_samples = set(merged_statistics)
    missing_from_json = tsv_samples - json_samples
    missing_from_tsv = json_samples - tsv_samples

    _write_sample_discrepancy_tsv(
        paths.sample_discrepancies_tsv,
        merged_statistics,
        tsv_samples,
        sample_counts,
        sample_batches,
    )

    if missing_from_json:
        logger.info(f"Found {len(missing_from_json)} unique TSV sample IDs missing from JSON:")
        for sample_id in sorted(missing_from_json):
            batch_ids = ", ".join(sorted(sample_batches[sample_id]))
            logger.info(f"  - {sample_id} (batch: {batch_ids})")

    if missing_from_tsv:
        logger.info(f"Found {len(missing_from_tsv)} JSON sample IDs not present in TSV.")

    if not missing_from_json and not missing_from_tsv:
        logger.info("TSV and JSON contain the same unique sample IDs.")

    logger.success(f"Sample discrepancies saved to: {paths.sample_discrepancies_tsv}")
