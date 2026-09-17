"""Centralized, cached batch/sample directory scanning for the merger."""

import json
from pathlib import Path
from typing import Any

from loguru import logger

from src.modules.statistics.merge.cache import load_cache, needs_rescan, save_cache, source_entry
from src.modules.statistics.merge.models import DirectoryScanResult, MergePaths


class DirectoryScanner:
    """
    Scan regenerate, FASTA and archive directories for per-batch sample IDs,
    cached per source by file metadata.

    Regenerate sample IDs are read from each batch's statistic_fullbatch.json keys;
    batches without a statistics JSON fall back to listing their sample directories.
    """

    def __init__(self, paths: MergePaths, excluded_batches: list[str] | None = None) -> None:
        self.paths = paths
        self.excluded_batches = set(excluded_batches or [])
        self._scan_cache: DirectoryScanResult | None = None
        self._cache_data: dict[str, dict[str, dict[str, Any]]] = {}
        self._force_refresh = False

    def _is_excluded(self, batch_name: str) -> bool:
        """Check if batch should be excluded."""
        return batch_name in self.excluded_batches

    def _scan_regenerate_batches(self) -> dict[str, set[str]]:
        """
        Read per-batch regenerate sample IDs from statistic_fullbatch.json keys,
        reusing cached IDs for unchanged sources.

        Returns:
            Dictionary mapping batch names to sets of sample IDs
        """
        result: dict[str, set[str]] = {}
        json_batches = 0
        dir_batches = 0
        cached_batches = 0
        regenerate_cache = self._cache_data.setdefault("regenerate", {})

        regenerate_dir = self.paths.regenerate_dir
        if not regenerate_dir.exists():
            logger.warning(f"Regenerate directory {regenerate_dir} does not exist")
            return result

        for batch_folder in regenerate_dir.iterdir():
            if not batch_folder.is_dir() or self._is_excluded(batch_folder.name):
                continue

            batch_name = batch_folder.name
            json_file = batch_folder / "statistic_fullbatch.json"
            source_path = json_file if json_file.exists() else batch_folder
            cached_batch = regenerate_cache.get(batch_name, {})

            if not needs_rescan(source_path, cached_batch, payload_key="samples", force_refresh=self._force_refresh):
                result[batch_name] = set(cached_batch.get("samples", []))
                cached_batches += 1
                continue

            if json_file.exists():
                try:
                    with json_file.open() as f:
                        batch_data = json.load(f)
                except (OSError, json.JSONDecodeError) as e:
                    logger.warning(f"Error reading JSON file {json_file}: {e}, falling back to directory scan")
                    sample_ids = {path.name for path in batch_folder.iterdir() if path.is_dir()}
                    result[batch_name] = sample_ids
                    regenerate_cache[batch_name] = source_entry(batch_folder, sorted(sample_ids), payload_key="samples")
                    dir_batches += 1
                    continue

                # Use JSON keys as sample IDs (much faster than directory scanning)
                sample_ids = set(batch_data.keys())
                result[batch_name] = sample_ids
                regenerate_cache[batch_name] = source_entry(json_file, sorted(sample_ids), payload_key="samples")
                json_batches += 1
            else:
                # No statistics JSON yet; count sample directories instead
                sample_ids = {path.name for path in batch_folder.iterdir() if path.is_dir()}
                result[batch_name] = sample_ids
                regenerate_cache[batch_name] = source_entry(batch_folder, sorted(sample_ids), payload_key="samples")
                dir_batches += 1

        total_batches = json_batches + dir_batches + cached_batches
        if total_batches > 0:
            logger.info(
                f"Regenerate scan: {cached_batches} cached, {json_batches} from JSON, {dir_batches} from directories",
            )

        return result

    def _scan_file_batches(self, cache_key: str, directory: Path, file_pattern: str) -> dict[str, set[str]]:
        """
        Scan per-batch file outputs (for example FASTA), reusing cached sample
        IDs for unchanged batch directories.

        Args:
            cache_key: Scan-cache category for this directory.
            directory: Directory containing one subdirectory per batch.
            file_pattern: Pattern sample files match (for example "*.fasta").

        Returns:
            Dictionary mapping batch names to sets of sample IDs
        """
        result: dict[str, set[str]] = {}
        category_cache = self._cache_data.setdefault(cache_key, {})
        cached_batches = 0
        scanned_batches = 0

        if not directory.exists():
            logger.warning(f"Directory {directory} does not exist")
            return result

        for batch_folder in directory.iterdir():
            if not batch_folder.is_dir() or self._is_excluded(batch_folder.name):
                continue

            cached_batch = category_cache.get(batch_folder.name, {})
            if not needs_rescan(batch_folder, cached_batch, payload_key="samples", force_refresh=self._force_refresh):
                cached_samples = set(cached_batch.get("samples", []))
                if cached_samples:
                    result[batch_folder.name] = cached_samples
                cached_batches += 1
                continue

            sample_ids = {path.stem for path in batch_folder.glob(file_pattern) if path.is_file()}
            if sample_ids:
                result[batch_folder.name] = sample_ids
            category_cache[batch_folder.name] = source_entry(batch_folder, sorted(sample_ids), payload_key="samples")
            scanned_batches += 1

        logger.info(f"{cache_key} scan: {cached_batches} cached, {scanned_batches} scanned")
        return result

    def scan_all_directories(self, *, force_refresh: bool = False) -> DirectoryScanResult:
        """
        Scan all relevant directories and cache results.

        Args:
            force_refresh: Force refresh of cached results

        Returns:
            DirectoryScanResult containing all scan data
        """
        self._force_refresh = force_refresh

        if self._scan_cache is not None and not force_refresh:
            return self._scan_cache

        # Load cache from file
        self._cache_data = load_cache(self.paths.scan_cache)

        logger.info("Scanning directories for samples")

        # Scan regenerate directory using JSON files for fast sample counting (with caching)
        regenerate_samples = self._scan_regenerate_batches()

        # Scan FASTA directories with per-batch caching; archive results add to the same view
        fasta_samples = self._scan_file_batches("fasta", self.paths.fasta_dir, "*.fasta")
        archive_fasta_samples = self._scan_file_batches("archive_fasta", self.paths.archive_dir, "*.fasta")
        for batch_name, sample_ids in archive_fasta_samples.items():
            fasta_samples.setdefault(batch_name, set()).update(sample_ids)

        save_cache(self.paths.scan_cache, self._cache_data)

        self._scan_cache = DirectoryScanResult(
            regenerate_samples=regenerate_samples,
            fasta_samples=fasta_samples,
        )

        logger.info(
            f"Directory scan complete. Found {len(regenerate_samples)} batches with samples (from JSON), "
            f"{len(fasta_samples)} with FASTA files across active and archive results",
        )

        return self._scan_cache


def get_successful_samples(scan_result: DirectoryScanResult) -> dict[str, list[str]]:
    """
    Combine regenerate and FASTA scan results into per-batch sorted sample lists.

    Args:
        scan_result: Directory scan results from DirectoryScanner.scan_all_directories()

    Returns:
        Dictionary mapping batch names to lists of successful sample IDs
    """
    successful_samples: dict[str, list[str]] = {}

    # Combine regenerate and fasta samples
    all_batches = set(scan_result.regenerate_samples.keys()) | set(scan_result.fasta_samples.keys())

    for batch_name in sorted(all_batches):
        samples: set[str] = set()

        # Add samples from regenerate directory
        if batch_name in scan_result.regenerate_samples:
            samples.update(scan_result.regenerate_samples[batch_name])

        # Add samples from fasta directory
        if batch_name in scan_result.fasta_samples:
            samples.update(scan_result.fasta_samples[batch_name])

        if samples:
            successful_samples[batch_name] = sorted(samples)

    return successful_samples
