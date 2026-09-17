"""CLI for merging mtDNA batch metadata and results into consolidated exports."""

import argparse
import sys
from pathlib import Path

from loguru import logger

from src.modules.statistics.merge.config import (
    DEFAULT_DATA_DIR,
    DEFAULT_EXCLUDED_BATCHES,
    DEFAULT_RESULTS_DIR,
)
from src.modules.statistics.merge.metadata import process_excel_files, read_metadata_files
from src.modules.statistics.merge.models import MergePaths
from src.modules.statistics.merge.reporting import analyze_sample_discrepancies, batch_statistics, create_merged_tsv
from src.modules.statistics.merge.scanner import DirectoryScanner, get_successful_samples
from src.modules.statistics.merge.statistics_json import merge_statistics_json


def main() -> None:
    """
    Main function to orchestrate the data merging process.
    """
    parser = argparse.ArgumentParser(description="Merge mtDNA data from multiple sources")
    parser.add_argument("excluded_batches", nargs="*", default=DEFAULT_EXCLUDED_BATCHES, help="Batch names to exclude")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Root data directory containing metadata/ (default: %(default)s)",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Results root containing manual_pipeline/regenerate/, fasta/ and archive/ (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for merged outputs and the scan cache (default: <results-dir>/merged)",
    )
    parser.add_argument("--force-refresh", action="store_true", help="Force refresh of cached directory scans")

    args = parser.parse_args()
    excluded_batches = args.excluded_batches
    force_refresh = args.force_refresh
    output_dir = args.output_dir or args.results_dir / "merged"
    paths = MergePaths(data_dir=args.data_dir, results_dir=args.results_dir, output_dir=output_dir)

    logger.info(f"Merging mtDNA data: data {paths.data_dir}, results {paths.results_dir}, outputs {paths.output_dir}")

    try:
        if excluded_batches != DEFAULT_EXCLUDED_BATCHES:
            logger.info(f"Excluding the following batches: {', '.join(excluded_batches)}")

        if force_refresh:
            logger.info("Force refresh enabled - ignoring cache")

        # Process Excel files to TSV (only batches with missing or unreadable TSVs)
        process_excel_files(excluded_batches, paths)

        # Single comprehensive directory scan reused by every downstream step
        logger.info("Performing comprehensive directory scan")
        scanner = DirectoryScanner(paths, excluded_batches)
        scan_result = scanner.scan_all_directories(force_refresh=force_refresh)

        # Load successful samples and metadata using cached scan results
        logger.info("Processing successful samples from cached scan results")
        successful_samples = get_successful_samples(scan_result)

        logger.info("Reading metadata files")
        metadata = read_metadata_files(excluded_batches, paths, force_refresh=force_refresh)

        # Create merged outputs
        logger.info("Creating merged TSV file")
        samples_per_batch, total_samples = create_merged_tsv(successful_samples, metadata, paths)

        # Print statistics using cached scan results
        batch_statistics(samples_per_batch, successful_samples, metadata, scan_result, paths)

        # Merge JSON statistics
        merged_statistics, _batch_count = merge_statistics_json(excluded_batches, paths, force_refresh=force_refresh)

        # Analyze discrepancies
        analyze_sample_discrepancies(merged_statistics, total_samples, paths)

        logger.info("-" * 150)

    except OSError as e:
        logger.error(f"Unexpected error: {e}")
        logger.exception("Traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()
