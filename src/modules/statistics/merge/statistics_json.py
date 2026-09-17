"""Merge per-batch statistics JSON files into a manifest-validated merged export."""

import json
from datetime import datetime
from pathlib import Path

from loguru import logger

from src.modules.statistics.merge.models import MergePaths, SourceSignatureEntry


def _statistics_json_files(excluded_batches: set[str], regenerate_dir: Path) -> list[Path]:
    """Return batch statistics JSON files that should be included in the merge."""
    if not regenerate_dir.exists():
        return []

    return sorted(
        [
            batch_folder / "statistic_fullbatch.json"
            for batch_folder in regenerate_dir.iterdir()
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
    output_json: Path,
    output_manifest: Path,
    source_signature: dict[str, SourceSignatureEntry],
    excluded_batches: set[str],
    *,
    force_refresh: bool = False,
) -> tuple[dict[str, dict], int] | None:
    """Load the existing merged statistics JSON when all source inputs are unchanged."""
    if force_refresh or not output_json.exists() or not output_manifest.exists():
        return None

    try:
        with output_manifest.open() as f:
            manifest = json.load(f)

        if manifest.get("excluded_batches") != sorted(excluded_batches):
            return None
        if manifest.get("source_files") != source_signature:
            return None

        output_mtime_ns = output_json.stat().st_mtime_ns
        if any(source["mtime_ns"] > output_mtime_ns for source in source_signature.values()):
            return None

        with output_json.open() as f:
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
    excluded_batches: list[str],
    paths: MergePaths,
    *,
    force_refresh: bool = False,
) -> tuple[dict[str, dict], int]:
    """
    Find all statistic_fullbatch.json files in the regenerate folder,
    merge them into a single JSON file, and save it under the output directory.

    The output format uses sample_id as the key without combining with batch name.

    Args:
        excluded_batches: List of batch names to exclude
        paths: Merge input/output path contract
        force_refresh: Rebuild the merged export even when the manifest is fresh

    Returns:
        The merged statistics data and number of batches processed
    """
    merged_data: dict[str, dict] = {}
    excluded_batch_set = set(excluded_batches)
    batch_count = 0

    regenerate_dir = paths.regenerate_dir
    if not regenerate_dir.exists():
        logger.error(f"Regenerate directory {regenerate_dir} not found")
        return merged_data, 0

    json_files = _statistics_json_files(excluded_batch_set, regenerate_dir)
    try:
        source_signature = _statistics_source_signature(json_files)
    except OSError as e:
        logger.error(f"Error reading statistics source metadata: {e!s}")
        return merged_data, 0

    cached_statistics = _load_merged_statistics_if_fresh(
        paths.merged_statistics_json,
        paths.merged_statistics_manifest,
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

            for sample_id, sample_data in batch_data.items():
                sample_data["batch"] = batch_name
                merged_data[sample_id] = sample_data

            batch_count += 1

        except (OSError, json.JSONDecodeError, ValueError) as e:
            logger.error(f"Error processing statistics for batch {batch_name}: {e!s}")

    output_json = paths.merged_statistics_json
    output_json.parent.mkdir(parents=True, exist_ok=True)

    with output_json.open("w") as f:
        json.dump(merged_data, f, indent=2)

    with paths.merged_statistics_manifest.open("w") as f:
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
