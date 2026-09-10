#!/usr/bin/env python3

import csv
import json
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from os import environ
from pathlib import Path
from typing import Any, TypedDict

import dotenv
import pandas as pd
from loguru import logger

# Load environment variables from .env file
dotenv.load_dotenv()

_MIN_PARTS = 2
_SCAN_CACHE_VERSION = 2


# Dataclasses for structured data
@dataclass
class BatchSummary:
    """Summary of samples in a batch across different processing stages."""

    batch_name: str
    successful_count: int = 0
    metadata_count: int = 0
    regenerate_count: int = 0
    fasta_count: int = 0

    @property
    def status(self) -> str:
        """Determine processing status based on counts."""
        if self.metadata_count == 0:
            return "No Meta"
        if self.fasta_count == 0:
            return "No FASTA"
        if self.successful_count == 0:
            return "No Pass"
        if self.metadata_count == self.successful_count:
            return "OK"
        return "Mismatch"


@dataclass
class DirectoryScanResult:
    """Result of scanning directories for samples."""

    regenerate_samples: dict[str, set[str]]
    fasta_samples: dict[str, set[str]]


class SourceSignatureEntry(TypedDict):
    """File metadata used to validate merged-statistics cache freshness."""

    path: str
    mtime_ns: int
    ctime_ns: int
    size: int


class ConfigurationError(Exception):
    """Raised when environment configuration is invalid."""


def validate_environment() -> None:
    """
    Validate required environment variables are set.

    Raises:
        ConfigurationError: If required environment variables are missing or invalid
    """
    required_vars = ["DATA_DIR", "RESULTS_DIR"]
    missing_vars = []

    for var in required_vars:
        value = environ.get(var)
        if not value:
            missing_vars.append(var)
        elif not Path(value).exists():
            logger.warning(f"Path {value} for {var} does not exist")

    if missing_vars:
        msg = f"Missing required environment variables: {', '.join(missing_vars)}"
        raise ConfigurationError(msg)


# Validate environment and get paths
validate_environment()
DATA_DIR = Path(environ.get("DATA_DIR"))  # type: ignore[reportArgumentType]
RESULTS_DIR = Path(environ.get("RESULTS_DIR"))  # type: ignore[reportArgumentType]

METADATA_DIR = DATA_DIR / "metadata"
REGENERATE_DIR = RESULTS_DIR / "manual_pipeline/regenerate"
FASTA_DIR = RESULTS_DIR / "fasta"
ARCHIVE_DIR = RESULTS_DIR / "archive"
OUTPUT_FILE = RESULTS_DIR / "merged/merged_metadata.tsv"
OUTPUT_JSON = RESULTS_DIR / "merged/merged_statistics.json"
OUTPUT_JSON_MANIFEST = RESULTS_DIR / "merged/merged_statistics.manifest.json"
OUTPUT_BATCH_STATS = RESULTS_DIR / "merged/batch_statistics.tsv"
CACHE_FILE = RESULTS_DIR / "merged/scan_cache.json"

# Default list of batches to exclude
DEFAULT_EXCLUDED_BATCHES = [
    "20250429_mtDNA_100",
]


class DirectoryScanner:
    """
    Centralized directory scanning with caching and JSON-based optimization for sample counting.
    """

    def __init__(self, excluded_batches: list[str] | None = None) -> None:
        self.excluded_batches = set(excluded_batches or [])
        self._scan_cache: DirectoryScanResult | None = None
        self._cache_data: dict[str, dict[str, dict[str, Any]]] = {}
        self._force_refresh = False

    def _is_excluded(self, batch_name: str) -> bool:
        """Check if batch should be excluded."""
        return batch_name in self.excluded_batches

    def _load_cache(self) -> dict[str, dict[str, dict[str, Any]]]:
        """
        Load scan cache from JSON file.

        Returns:
            Dictionary mapping scan categories to batch cache data.
        """
        cache: dict[str, dict[str, dict[str, Any]]] = {}
        if CACHE_FILE.exists():
            try:
                with CACHE_FILE.open() as f:
                    data = json.load(f)

                if data.get("version") != _SCAN_CACHE_VERSION:
                    logger.info("Ignoring scan cache with incompatible schema version")
                    return cache

                raw_cache = data.get("scans", {})
                if not isinstance(raw_cache, dict):
                    logger.info("Ignoring malformed scan cache")
                    return cache

                cache = raw_cache
                cached_batches = sum(len(category_cache) for category_cache in cache.values())
                logger.info(f"Loaded cache with {cached_batches} batch entries")
            except (OSError, json.JSONDecodeError) as e:
                logger.warning(f"Error loading cache file: {e}")
        return cache

    def _save_cache(self, cache_data: dict[str, dict[str, dict[str, Any]]]) -> None:
        """
        Save scan cache to JSON file.

        Args:
            cache_data: Dictionary mapping scan categories to batch cache data.
        """
        CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
        data = {"version": _SCAN_CACHE_VERSION, "last_scan": datetime.now().isoformat(), "scans": cache_data}  # noqa: DTZ005
        try:
            with CACHE_FILE.open("w") as f:
                json.dump(data, f, indent=2)
            cached_batches = sum(len(category_cache) for category_cache in cache_data.values())
            logger.debug(f"Saved cache with {cached_batches} batch entries")
        except OSError as e:
            logger.warning(f"Error saving cache file: {e}")

    def _source_cache_entry(self, source_path: Path, sample_ids: set[str]) -> dict[str, Any]:
        """Build a cache entry for sample IDs derived from a source file or directory."""
        source_stat = source_path.stat()
        return {
            "source_path": str(source_path),
            "source_type": "file" if source_path.is_file() else "directory",
            "mtime_ns": source_stat.st_mtime_ns,
            "ctime_ns": source_stat.st_ctime_ns,
            "size": source_stat.st_size,
            "samples": sorted(sample_ids),
        }

    def _needs_rescan(self, source_path: Path, cached_data: dict[str, Any]) -> bool:
        """
        Check if a cached source needs to be re-scanned based on source metadata.

        Args:
            source_path: Directory or file whose entries/data were cached
            cached_data: Cached data for this source

        Returns:
            True if source needs re-scanning
        """
        needs_rescan = self._force_refresh or not cached_data
        if needs_rescan:
            return True

        try:
            source_stat = source_path.stat()
        except OSError:
            return True

        source_type = "file" if source_path.is_file() else "directory"
        required_keys = {"source_path", "source_type", "mtime_ns", "ctime_ns", "size", "samples"}
        if not required_keys.issubset(cached_data):
            return True

        return (
            cached_data.get("source_path") != str(source_path)
            or cached_data.get("source_type") != source_type
            or cached_data.get("mtime_ns") != source_stat.st_mtime_ns
            or cached_data.get("ctime_ns") != source_stat.st_ctime_ns
            or cached_data.get("size") != source_stat.st_size
        )

    def _scan_regenerate_fast(self) -> dict[str, set[str]]:
        """
        Fast regenerate directory scanning using JSON files to count samples.
        Uses caching to skip unchanged batches based on modification time.

        Returns:
            Dictionary mapping batch names to sets of sample IDs
        """
        result = {}
        json_batches = 0
        dir_batches = 0
        cached_batches = 0
        regenerate_cache = self._cache_data.setdefault("regenerate", {})

        if not REGENERATE_DIR.exists():
            logger.warning(f"Regenerate directory {REGENERATE_DIR} does not exist")
            return result

        try:
            for batch_folder in REGENERATE_DIR.iterdir():
                if not batch_folder.is_dir() or self._is_excluded(batch_folder.name):
                    continue

                batch_name = batch_folder.name
                json_file = batch_folder / "statistic_fullbatch.json"

                source_path = json_file if json_file.exists() else batch_folder
                cached_batch = regenerate_cache.get(batch_name, {})

                if not self._needs_rescan(source_path, cached_batch):
                    cached_samples = cached_batch.get("samples", [])
                    result[batch_name] = set(cached_samples)
                    cached_batches += 1
                    continue

                if json_file.exists():
                    try:
                        with json_file.open() as f:
                            batch_data = json.load(f)

                        # Use JSON keys as sample IDs (much faster than directory scanning)
                        sample_ids = set(batch_data.keys())
                        result[batch_name] = sample_ids
                        json_batches += 1

                        regenerate_cache[batch_name] = self._source_cache_entry(json_file, sample_ids)

                    except (OSError, json.JSONDecodeError) as e:
                        logger.warning(f"Error reading JSON file {json_file}: {e}, falling back to directory scan")
                        # Fallback to directory scanning for this batch
                        sample_dirs = [d for d in batch_folder.iterdir() if d.is_dir()]
                        sample_ids = {d.name for d in sample_dirs}
                        result[batch_name] = sample_ids
                        dir_batches += 1

                        regenerate_cache[batch_name] = self._source_cache_entry(batch_folder, sample_ids)
                else:
                    # No JSON file, scan directories
                    sample_dirs = [d for d in batch_folder.iterdir() if d.is_dir()]
                    sample_ids = {d.name for d in sample_dirs}
                    result[batch_name] = sample_ids
                    dir_batches += 1

                    regenerate_cache[batch_name] = self._source_cache_entry(batch_folder, sample_ids)

        except (OSError, ValueError) as e:
            logger.warning(f"Error scanning regenerate directory: {e}, falling back to subprocess")
            return self._scan_directory_fast(REGENERATE_DIR, "*")

        # Save updated cache
        self._save_cache(self._cache_data)

        total_batches = json_batches + dir_batches + cached_batches
        if total_batches > 0:
            logger.info(
                f"Regenerate scan: {cached_batches} cached, {json_batches} from JSON, {dir_batches} from directories",
            )

        return result

    def _scan_directory_python(self, directory: Path, file_pattern: str = "*") -> dict[str, set[str]]:
        """Scan two-level batch/sample directories with Python."""
        result: dict[str, set[str]] = {}
        for batch_folder in directory.iterdir():
            if not batch_folder.is_dir() or self._is_excluded(batch_folder.name):
                continue

            if file_pattern == "*":
                sample_ids = {path.name for path in batch_folder.iterdir() if path.is_dir()}
            else:
                sample_ids = {path.stem for path in batch_folder.glob(file_pattern) if path.is_file()}

            if sample_ids:
                result[batch_folder.name] = sample_ids

        return result

    def _scan_directory_fast(self, directory: Path, file_pattern: str = "*") -> dict[str, set[str]]:
        """
        Fast directory scanning using subprocess find command.

        Args:
            directory: Directory to scan
            file_pattern: Pattern to match (e.g., "*.fasta")

        Returns:
            Dictionary mapping batch names to sets of sample IDs
        """
        result: dict[str, set[str]] = {}

        if not directory.exists():
            logger.warning(f"Directory {directory} does not exist")
            return result

        try:
            # Use find command for faster directory traversal
            if file_pattern == "*":
                # For directories (regenerate samples)
                cmd = ["find", str(directory), "-mindepth", "2", "-maxdepth", "2", "-type", "d"]
            else:
                # For files (for example, FASTA)
                cmd = ["find", str(directory), "-mindepth", "2", "-maxdepth", "2", "-name", file_pattern, "-type", "f"]

            process = subprocess.run(cmd, capture_output=True, text=True, check=True)

            for line in process.stdout.splitlines():
                path = Path(line)
                batch_name = path.parent.name

                if self._is_excluded(batch_name):
                    continue

                if batch_name not in result:
                    result[batch_name] = set()

                # Extract sample ID (filename without extension or directory name)
                sample_id = path.stem if path.is_file() else path.name
                result[batch_name].add(sample_id)

        except subprocess.CalledProcessError as e:
            logger.warning(f"Fast directory scan failed for {directory}, falling back to Python: {e}")
            return self._scan_directory_python(directory, file_pattern)
        except FileNotFoundError:
            logger.warning("'find' command not available, using Python fallback")
            return self._scan_directory_python(directory, file_pattern)

        return result

    def _scan_file_batches_cached(self, cache_key: str, directory: Path, file_pattern: str) -> dict[str, set[str]]:
        """Scan batch file outputs, reusing cached sample IDs for unchanged batch directories."""
        result: dict[str, set[str]] = {}
        category_cache = self._cache_data.setdefault(cache_key, {})
        cached_batches = 0
        scanned_batches = 0

        if not directory.exists():
            logger.warning(f"Directory {directory} does not exist")
            return result

        try:
            if self._force_refresh or not category_cache:
                result = self._scan_directory_fast(directory, file_pattern)
                batch_folders = [
                    batch_folder
                    for batch_folder in directory.iterdir()
                    if batch_folder.is_dir() and not self._is_excluded(batch_folder.name)
                ]
                for batch_folder in batch_folders:
                    sample_ids = result.get(batch_folder.name, set())
                    category_cache[batch_folder.name] = self._source_cache_entry(batch_folder, sample_ids)
                logger.info(f"{cache_key} scan: 0 cached, {len(batch_folders)} scanned with find")
                self._save_cache(self._cache_data)
                return result

            for batch_folder in directory.iterdir():
                if not batch_folder.is_dir() or self._is_excluded(batch_folder.name):
                    continue

                cached_batch = category_cache.get(batch_folder.name, {})
                if not self._needs_rescan(batch_folder, cached_batch):
                    cached_samples = set(cached_batch.get("samples", []))
                    if cached_samples:
                        result[batch_folder.name] = cached_samples
                    cached_batches += 1
                    continue

                sample_ids = {path.stem for path in batch_folder.glob(file_pattern) if path.is_file()}
                if sample_ids:
                    result[batch_folder.name] = sample_ids
                category_cache[batch_folder.name] = self._source_cache_entry(batch_folder, sample_ids)
                scanned_batches += 1
        except OSError as e:
            logger.warning(f"Cached scan failed for {directory}, falling back to find: {e}")
            return self._scan_directory_fast(directory, file_pattern)

        logger.info(f"{cache_key} scan: {cached_batches} cached, {scanned_batches} scanned")
        self._save_cache(self._cache_data)
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
        self._cache_data = self._load_cache()

        logger.info("Scanning directories for samples")

        # Scan regenerate directory using JSON files for fast sample counting (with caching)
        regenerate_samples = self._scan_regenerate_fast()

        # Scan FASTA directories with per-batch caching.
        fasta_samples = self._scan_file_batches_cached("fasta", FASTA_DIR, "*.fasta")
        archive_fasta_samples = self._scan_file_batches_cached("archive_fasta", ARCHIVE_DIR, "*.fasta")
        for batch_name, sample_ids in archive_fasta_samples.items():
            fasta_samples.setdefault(batch_name, set()).update(sample_ids)

        self._scan_cache = DirectoryScanResult(
            regenerate_samples=regenerate_samples,
            fasta_samples=fasta_samples,
        )

        logger.info(
            f"Directory scan complete. Found {len(regenerate_samples)} batches with samples (from JSON), "
            f"{len(fasta_samples)} with FASTA files across active and archive results",
        )

        return self._scan_cache


def get_successful_samples(
    excluded_batches: list[str] | None = None,
    scan_result: DirectoryScanResult | None = None,
) -> dict[str, list[str]]:
    """
    Get a dictionary of batch names and their successful samples
    by scanning both the regenerate and fasta folders.

    Args:
        excluded_batches: List of batch names to exclude
        scan_result: Pre-computed directory scan results (optional, will scan if not provided)

    Returns:
        Dictionary mapping batch names to lists of successful sample IDs
    """
    # Use provided scan result or create new scanner (for backward compatibility)
    if scan_result is None:
        scanner = DirectoryScanner(excluded_batches)
        scan_result = scanner.scan_all_directories()

    successful_samples = {}

    # Combine regenerate and fasta samples
    all_batches = set(scan_result.regenerate_samples.keys()) | set(scan_result.fasta_samples.keys())

    for batch_name in sorted(all_batches):
        samples = set()

        # Add samples from regenerate directory
        if batch_name in scan_result.regenerate_samples:
            samples.update(scan_result.regenerate_samples[batch_name])

        # Add samples from fasta directory
        if batch_name in scan_result.fasta_samples:
            samples.update(scan_result.fasta_samples[batch_name])

        if samples:
            successful_samples[batch_name] = sorted(samples)

    return successful_samples


def read_metadata_files(excluded_batches: list[str] | None = None) -> dict[str, dict[str, str]]:
    """
    Read all TSV files from metadata directory and return a dictionary
    mapping batch name to a dict of sample_id -> barcode

    Args:
        excluded_batches: List of batch names to exclude
    """
    metadata = {}
    excluded_batches = excluded_batches or []

    metadata_path = Path(METADATA_DIR)
    if not metadata_path.exists():
        logger.error(f"Metadata directory {METADATA_DIR} not found")
        sys.exit(1)

    tsv_files = list(metadata_path.glob("*.tsv"))

    for tsv_file in tsv_files:
        batch_name = tsv_file.stem
        if batch_name in excluded_batches:
            continue

        metadata[batch_name] = {}

        with Path(tsv_file).open() as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) >= _MIN_PARTS:
                        sample_id = parts[0]
                        barcode = parts[1]
                        metadata[batch_name][sample_id] = barcode

    return metadata


def extract_batch_id(batch_name: str) -> str:
    """
    Extract a clean batch ID from a batch name.

    Args:
        batch_name: The full batch name (e.g., '20250509_mtDNA_44')

    Returns:
        A clean batch ID
    """
    # In this implementation, we return the batch name as is
    # If specific cleaning is needed, it can be implemented here
    return batch_name


def create_merged_tsv(
    successful_samples: dict[str, list[str]],
    metadata: dict[str, dict[str, str]],
) -> tuple[dict[str, int], int]:
    """
    Create a merged TSV file with batch, sample_id, and barcode for all successful samples
    """
    samples_per_batch = {}
    total_samples = 0

    with OUTPUT_FILE.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["batch", "sample_id", "barcode"])

        for batch_name, samples in sorted(successful_samples.items()):
            samples_per_batch[batch_name] = 0

            if batch_name in metadata:
                batch_metadata = metadata[batch_name]

                for sample_id in sorted(samples):
                    if sample_id in batch_metadata:
                        barcode = batch_metadata[sample_id]
                        clean_batch_name = extract_batch_id(batch_name)
                        writer.writerow([clean_batch_name, sample_id, barcode])
                        samples_per_batch[batch_name] += 1
                        total_samples += 1
    return samples_per_batch, total_samples


def batch_statistics(
    samples_per_batch: dict[str, int],
    successful_samples: dict[str, list[str]],
    metadata: dict[str, dict[str, str]],
    scan_result: DirectoryScanResult | None = None,
) -> None:
    """
    Print statistics about the number of samples in each batch with separate columns for regenerate and FASTA.
    Also saves the statistics to a TSV file.

    Args:
        samples_per_batch: Count of successful samples per batch
        successful_samples: Dictionary of successful samples by batch
        metadata: Metadata dictionary by batch
        scan_result: Pre-computed directory scan results (optional, will scan if not provided)
    """
    # Use provided scan result or create new scanner
    if scan_result is None:
        scanner = DirectoryScanner()
        scan_result = scanner.scan_all_directories()

    logger.info("Batch Statistics:")
    logger.info("-" * 150)
    header = (
        f"{'Batch Name':<25} {'Successful':<20} {'Metadata':<20} "
        f"{'Regenerate':<20} {'Fasta':<20} {'Status':<20}"
    )
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
    OUTPUT_BATCH_STATS.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_BATCH_STATS.open("w", newline="") as f:
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

    logger.success(f"Batch statistics saved to: {OUTPUT_BATCH_STATS}")


def _statistics_json_files(excluded_batches: set[str]) -> list[Path]:
    """Return batch statistics JSON files that should be included in the merge."""
    regenerate_path = Path(REGENERATE_DIR)
    if not regenerate_path.exists():
        return []

    return sorted(
        [
            batch_folder / "statistic_fullbatch.json"
            for batch_folder in regenerate_path.iterdir()
            if batch_folder.is_dir()
            and batch_folder.name not in excluded_batches
            and (batch_folder / "statistic_fullbatch.json").exists()
        ],
        key=lambda path: path.parent.name,
    )


def _statistics_source_signature(json_files: list[Path]) -> dict[str, SourceSignatureEntry]:
    """Build a source-file signature that detects changed, added, and removed batch JSON files."""
    source_signature: dict[str, SourceSignatureEntry] = {}
    for json_file in json_files:
        file_stat = json_file.stat()
        source_signature[json_file.parent.name] = {
            "path": str(json_file),
            "mtime_ns": file_stat.st_mtime_ns,
            "ctime_ns": file_stat.st_ctime_ns,
            "size": file_stat.st_size,
        }
    return source_signature


def _load_merged_statistics_if_fresh(
    source_signature: dict[str, SourceSignatureEntry],
    excluded_batches: set[str],
    *,
    force_refresh: bool = False,
) -> tuple[dict[str, dict], int] | None:
    """Load the existing merged statistics JSON when all source inputs are unchanged."""
    if force_refresh or not OUTPUT_JSON.exists() or not OUTPUT_JSON_MANIFEST.exists():
        return None

    try:
        with OUTPUT_JSON_MANIFEST.open() as f:
            manifest = json.load(f)

        if manifest.get("excluded_batches") != sorted(excluded_batches):
            return None
        if manifest.get("source_files") != source_signature:
            return None

        output_mtime_ns = OUTPUT_JSON.stat().st_mtime_ns
        if any(source["mtime_ns"] > output_mtime_ns for source in source_signature.values()):
            return None

        with OUTPUT_JSON.open() as f:
            merged_data = json.load(f)
    except (OSError, json.JSONDecodeError, ValueError) as e:
        logger.warning(f"Could not reuse merged statistics cache: {e}")
        return None
    else:
        batch_count = manifest.get("batch_count", len(source_signature))
        logger.success(
            f"Reused merged statistics file with {len(merged_data)} total samples from {batch_count} batches",
        )
        return merged_data, batch_count


def merge_statistics_json(
    excluded_batches: list[str] | None = None,
    *,
    force_refresh: bool = False,
) -> tuple[dict[str, dict], int]:
    """
    Find all statistic_fullbatch.json files in the regenerate folder,
    merge them into a single JSON file, and save to OUTPUT_JSON.

    The output format uses sample_id as the key without combining with batch name.

    Args:
        excluded_batches: List of batch names to exclude

    Returns:
        The merged statistics data and number of batches processed
    """
    merged_data = {}
    excluded_batch_set = set(excluded_batches or [])
    batch_count = 0

    regenerate_path = Path(REGENERATE_DIR)
    if not regenerate_path.exists():
        logger.error(f"Regenerate directory {REGENERATE_DIR} not found")
        return merged_data, 0

    json_files = _statistics_json_files(excluded_batch_set)
    try:
        source_signature = _statistics_source_signature(json_files)
    except OSError as e:
        logger.error(f"Error reading statistics source metadata: {e!s}")
        return merged_data, 0

    cached_statistics = _load_merged_statistics_if_fresh(
        source_signature,
        excluded_batch_set,
        force_refresh=force_refresh,
    )
    if cached_statistics is not None:
        return cached_statistics

    for json_file in json_files:
        batch_name = json_file.parent.name
        try:
            with json_file.open() as f:
                batch_data = json.load(f)

            clean_batch_name = extract_batch_id(batch_name)

            for sample_id, sample_data in batch_data.items():
                sample_data["batch"] = clean_batch_name
                merged_data[sample_id] = sample_data

            batch_count += 1

        except (OSError, json.JSONDecodeError, ValueError) as e:
            logger.error(f"Error processing statistics for batch {batch_name}: {e!s}")

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_JSON.open("w") as f:
        json.dump(merged_data, f, indent=2)

    with OUTPUT_JSON_MANIFEST.open("w") as f:
        json.dump(
            {
                "created_at": datetime.now().isoformat(),  # noqa: DTZ005
                "excluded_batches": sorted(excluded_batch_set),
                "source_files": source_signature,
                "batch_count": batch_count,
            },
            f,
            indent=2,
        )

    logger.success(f"Created merged statistics file with {len(merged_data)} total samples from {batch_count} batches")
    return merged_data, batch_count


def export_tsv_data(  # noqa: C901,PLR0911,PLR0912
    input_file: str,
    output_file: str,
    excluded_batches: list[str] | None = None,
    barcode_len: int | None = None,
) -> int:
    """
    Export Excel or TSV data to TSV format, ensuring numeric IDs and barcodes
    preserve leading zeros. Uses the same approach as extract_sample_info.

    Args:
        input_file: Path to the input Excel or TSV file
        output_file: Path to the output TSV file
        excluded_batches: List of batch names to exclude
        barcode_len: Target length for barcodes (if None, will be determined from data)

    Returns:
        Number of samples exported

    Raises:
        FileNotFoundError: If the input file doesn't exist
        pd.errors.EmptyDataError: If the input file is empty
        ValueError: If required columns are missing
    """
    excluded_batches = excluded_batches or []

    # Skip processing if TSV file already exists
    if Path(output_file).exists():
        try:
            with Path(output_file).open() as f:
                sample_count = sum(1 for line in f if line.strip())
        except OSError:
            logger.warning(f"TSV file exists but couldn't read it: {output_file}, regenerating")
        else:
            return sample_count

    try:
        if input_file.lower().endswith(".xlsx") or input_file.lower().endswith(".xls"):
            df = pd.read_excel(input_file, dtype=str)
        else:
            df = pd.read_csv(input_file, sep="\t", dtype=str)

        if excluded_batches and "batch" in df.columns:
            df = df[~df["batch"].isin(excluded_batches)]

        if input_file.lower().endswith(".xlsx") or input_file.lower().endswith(".xls"):
            lid_cols = df.filter(like="LID").columns
            barcode_cols = df.filter(like="BARCODE").columns

            if lid_cols.empty or barcode_cols.empty:
                logger.warning(f"Could not find LID or BARCODE columns in {input_file}")
                return 0

            lid_col = lid_cols[0]
            barcode_col = barcode_cols[0]

            result_df = pd.DataFrame({"lid": df[lid_col], "barcode": df[barcode_col]})

            # Remove rows with missing LID or barcode
            result_df = result_df.dropna(subset=["lid", "barcode"])

            numeric_mask = result_df["barcode"].notna() & result_df["barcode"].str.match(r"^\d+$")
            if numeric_mask.any():
                if barcode_len is None:
                    max_len = result_df.loc[numeric_mask, "barcode"].str.len().max()
                    barcode_len = max_len

                result_df.loc[numeric_mask, "barcode"] = (
                    result_df.loc[numeric_mask, "barcode"].astype(str).str.zfill(barcode_len)
                )

            result_df.to_csv(output_file, sep="\t", index=False, header=False)

            return len(result_df)
        df.to_csv(output_file, sep="\t", index=False)
        logger.success(f"Successfully exported data to {output_file}")
        logger.info(f"Total samples: {len(df)}")
        return len(df)

    except FileNotFoundError:
        logger.error(f"Input file not found: {input_file}")
        return 0
    except pd.errors.EmptyDataError:
        logger.error(f"Input file is empty: {input_file}")
        return 0
    except (OSError, pd.errors.ParserError, ValueError, KeyError) as e:
        logger.error(f"Error exporting TSV file: {e!s}")
        logger.exception("Traceback:")
        return 0


def process_excel_files(excluded_batches: list[str]) -> int:
    """
    Process Excel files in metadata directory and convert to TSV.

    Args:
        excluded_batches: List of batch names to exclude

    Returns:
        Number of samples exported
    """
    metadata_path = Path(METADATA_DIR)
    excel_files = list(metadata_path.glob("*.xlsx"))
    logger.info(f"Found {len(excel_files)} Excel files in metadata directory.")

    total_exported = 0
    for excel_file in excel_files:
        batch_name = excel_file.stem

        if batch_name in excluded_batches:
            continue

        tsv_path = metadata_path / f"{batch_name}.tsv"
        samples_exported = export_tsv_data(str(excel_file), str(tsv_path), excluded_batches)
        total_exported += samples_exported

    if total_exported > 0:
        logger.info(f"Total samples exported from Excel to TSV: {total_exported}")
    else:
        logger.info("No Excel files were processed.")

    return total_exported


def analyze_sample_discrepancies(merged_statistics: dict[str, dict], total_samples: int) -> None:
    """
    Analyze and report discrepancies between TSV and JSON data.

    TSV counts include rows and may contain the same sample ID in multiple batches,
    while merged JSON uses unique sample IDs as keys. Compare unique IDs in both
    directions instead of treating the row-count difference as a missing-sample
    count.

    Args:
        merged_statistics: Merged JSON statistics data
        total_samples: Number of rows in TSV file
    """
    json_sample_count = len(merged_statistics)
    tsv_sample_count = total_samples

    logger.info(f"Merged JSON file created: {OUTPUT_JSON}. Total unique samples in merged JSON: {json_sample_count}")
    logger.info(f"Merged TSV file created: {OUTPUT_FILE}. Total rows in merged TSV: {tsv_sample_count}")

    # Read the TSV once so row counts, duplicate rows, and unique-ID differences
    # are all derived from the same data.
    tsv_samples: set[str] = set()
    sample_batches: dict[str, set[str]] = {}
    sample_counts: Counter[str] = Counter()
    try:
        with OUTPUT_FILE.open(newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)  # Skip header
            for parts in reader:
                if len(parts) >= _MIN_PARTS and parts[1]:
                    batch_id, sample_id = parts[0], parts[1]
                    tsv_samples.add(sample_id)
                    sample_counts[sample_id] += 1
                    sample_batches.setdefault(sample_id, set()).add(batch_id)
    except FileNotFoundError:
        logger.error(f"TSV file {OUTPUT_FILE} not found")
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

    if missing_from_json:
        logger.info(f"Found {len(missing_from_json)} unique TSV sample IDs missing from JSON:")
        for sample_id in sorted(missing_from_json):
            batch_ids = ", ".join(sorted(sample_batches[sample_id]))
            logger.info(f"  - {sample_id} (batch: {batch_ids})")

    if missing_from_tsv:
        logger.info(f"Found {len(missing_from_tsv)} JSON sample IDs not present in TSV.")

    if not missing_from_json and not missing_from_tsv:
        logger.info("TSV and JSON contain the same unique sample IDs.")


def main() -> None:
    """
    Main function to orchestrate the data merging process.
    """
    import argparse  # noqa: PLC0415

    parser = argparse.ArgumentParser(description="Merge mtDNA data from multiple sources")
    parser.add_argument("excluded_batches", nargs="*", default=DEFAULT_EXCLUDED_BATCHES, help="Batch names to exclude")
    parser.add_argument("--force-refresh", action="store_true", help="Force refresh of cached directory scans")

    args = parser.parse_args()
    excluded_batches = args.excluded_batches
    force_refresh = args.force_refresh

    try:
        if excluded_batches != DEFAULT_EXCLUDED_BATCHES:
            logger.info(f"Excluding the following batches: {', '.join(excluded_batches)}")

        if force_refresh:
            logger.info("Force refresh enabled - ignoring cache")

        # Process Excel files to TSV
        process_excel_files(excluded_batches)

        # Create directory scanner for efficient processing - SINGLE SCAN
        logger.info("Performing comprehensive directory scan")
        scanner = DirectoryScanner(excluded_batches)
        scan_result = scanner.scan_all_directories(force_refresh=force_refresh)

        # Load successful samples and metadata using cached scan results
        logger.info("Processing successful samples from cached scan results")
        successful_samples = get_successful_samples(excluded_batches, scan_result)

        logger.info("Reading metadata files")
        metadata = read_metadata_files(excluded_batches)

        # Create merged outputs
        logger.info("Creating merged TSV file")
        samples_per_batch, total_samples = create_merged_tsv(successful_samples, metadata)

        # Print statistics using cached scan results
        batch_statistics(samples_per_batch, successful_samples, metadata, scan_result)

        # Merge JSON statistics
        merged_statistics, _batch_count = merge_statistics_json(excluded_batches, force_refresh=force_refresh)

        # Analyze discrepancies
        analyze_sample_discrepancies(merged_statistics, total_samples)

        logger.info("-" * 150)

    except ConfigurationError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    except OSError as e:
        logger.error(f"Unexpected error: {e}")
        logger.exception("Traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()
