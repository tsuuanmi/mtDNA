"""Shared JSON scan-cache persistence and source-freshness validation."""

import json
from datetime import datetime
from pathlib import Path
from typing import Any

from loguru import logger

CACHE_VERSION = 3
_ENTRY_KEYS = {"source_path", "source_type", "mtime_ns", "ctime_ns", "size"}


def load_cache(cache_file: Path) -> dict[str, dict[str, dict[str, Any]]]:
    """
    Load per-category scan caches from the JSON cache file.

    Args:
        cache_file: Cache file to load (missing or incompatible files yield an empty cache).

    Returns:
        Dictionary mapping scan categories to batch cache data.
    """
    cache: dict[str, dict[str, dict[str, Any]]] = {}
    if cache_file.exists():
        try:
            with cache_file.open() as f:
                data = json.load(f)

            if data.get("version") != CACHE_VERSION:
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


def save_cache(cache_file: Path, cache_data: dict[str, dict[str, dict[str, Any]]]) -> None:
    """
    Persist per-category scan caches to the JSON cache file.

    Args:
        cache_file: Cache file to write; parent directories are created.
        cache_data: Dictionary mapping scan categories to batch cache data.
    """
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    data = {"version": CACHE_VERSION, "last_scan": datetime.now().isoformat(), "scans": cache_data}  # noqa: DTZ005
    try:
        with cache_file.open("w") as f:
            json.dump(data, f, indent=2)
        cached_batches = sum(len(category_cache) for category_cache in cache_data.values())
        logger.debug(f"Saved cache with {cached_batches} batch entries")
    except OSError as e:
        logger.warning(f"Error saving cache file: {e}")


def source_entry(source_path: Path, payload: list[Any], *, payload_key: str) -> dict[str, Any]:
    """
    Build a cache entry for data derived from a source file or directory.

    Args:
        source_path: File or directory the payload was derived from.
        payload: Derived data to cache alongside the source metadata.
        payload_key: Cache key under which the payload is stored.
    """
    source_stat = source_path.stat()
    return {
        "source_path": str(source_path),
        "source_type": "file" if source_path.is_file() else "directory",
        "mtime_ns": source_stat.st_mtime_ns,
        "ctime_ns": source_stat.st_ctime_ns,
        "size": source_stat.st_size,
        payload_key: payload,
    }


def needs_rescan(
    source_path: Path,
    cached_data: dict[str, Any],
    *,
    payload_key: str,
    force_refresh: bool = False,
) -> bool:
    """
    Check if a cached source needs to be re-scanned based on source metadata.

    Args:
        source_path: Directory or file whose entries/data were cached.
        cached_data: Cached data for this source.
        payload_key: Cache key that must hold the payload in a valid entry.
        force_refresh: Force re-scanning regardless of cache state.

    Returns:
        True if source needs re-scanning
    """
    if force_refresh or not cached_data:
        return True

    try:
        source_stat = source_path.stat()
    except OSError:
        return True

    source_type = "file" if source_path.is_file() else "directory"
    if not (_ENTRY_KEYS | {payload_key}).issubset(cached_data):
        return True

    return (
        cached_data.get("source_path") != str(source_path)
        or cached_data.get("source_type") != source_type
        or cached_data.get("mtime_ns") != source_stat.st_mtime_ns
        or cached_data.get("ctime_ns") != source_stat.st_ctime_ns
        or cached_data.get("size") != source_stat.st_size
    )
