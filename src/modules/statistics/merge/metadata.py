"""Metadata TSV/Excel ingestion with concurrent NAS reads and a persistent source cache."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from loguru import logger

from src.modules.statistics.merge.cache import load_cache, needs_rescan, save_cache, source_entry
from src.modules.statistics.merge.config import MIN_TSV_PARTS
from src.modules.statistics.merge.models import MergePaths

# Concurrent worker cap: metadata reads are latency-bound on the NAS share,
# not bandwidth-bound, so parallel round trips dominate the speedup.
_METADATA_WORKERS = 32


def _read_tsv_pairs(tsv_file: Path) -> dict[str, str]:
    """Read one metadata TSV into a sample_id -> barcode mapping."""
    pairs: dict[str, str] = {}
    with tsv_file.open() as f:
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= MIN_TSV_PARTS:
                    pairs[parts[0]] = parts[1]
    return pairs


def read_metadata_files(
    excluded_batches: list[str],
    paths: MergePaths,
    *,
    force_refresh: bool = False,
) -> dict[str, dict[str, str]]:
    """
    Read all TSV files from the metadata directory and return a dictionary
    mapping batch name to a dict of sample_id -> barcode.

    Unchanged TSVs are served from the persistent scan cache; fresh files are
    read concurrently because the metadata directory lives on a high-latency
    NAS share.

    Args:
        excluded_batches: List of batch names to exclude
        paths: Merge input/output path contract
        force_refresh: Re-read every TSV even when a cache entry is fresh

    Returns:
        Dictionary mapping batch names to sample_id -> barcode mappings
    """
    metadata_path = paths.metadata_dir
    if not metadata_path.exists():
        logger.error(f"Metadata directory {metadata_path} not found")
        sys.exit(1)

    excluded = set(excluded_batches)
    tsv_files = sorted(metadata_path.glob("*.tsv"))
    included = [tsv_file for tsv_file in tsv_files if tsv_file.stem not in excluded]

    cache_data = load_cache(paths.scan_cache)
    category_cache = cache_data.setdefault("metadata", {})

    metadata: dict[str, dict[str, str]] = {}
    fresh_files: list[Path] = []
    cached_batches = 0

    for tsv_file in included:
        cached_batch = category_cache.get(tsv_file.stem, {})
        if not needs_rescan(tsv_file, cached_batch, payload_key="pairs", force_refresh=force_refresh):
            metadata[tsv_file.stem] = dict(cached_batch.get("pairs", []))
            cached_batches += 1
        else:
            fresh_files.append(tsv_file)

    if fresh_files:
        with ThreadPoolExecutor(max_workers=min(_METADATA_WORKERS, len(fresh_files))) as pool:
            for tsv_file, pairs in zip(fresh_files, pool.map(_read_tsv_pairs, fresh_files), strict=True):
                metadata[tsv_file.stem] = pairs
                category_cache[tsv_file.stem] = source_entry(
                    tsv_file,
                    [list(pair) for pair in sorted(pairs.items())],
                    payload_key="pairs",
                )
        save_cache(paths.scan_cache, cache_data)

    logger.info(f"Metadata read: {cached_batches} cached, {len(fresh_files)} read")
    return metadata


def convert_excel_to_tsv(excel_file: Path, tsv_path: Path) -> int:
    """
    Export an Excel workbook's LID/barcode columns to a headerless TSV.
    Numeric barcodes are zero-padded to the longest observed barcode so IDs
    preserve leading zeros.

    Args:
        excel_file: Path to the input Excel (.xlsx) file
        tsv_path: Path to the output TSV file

    Returns:
        Number of samples exported
    """
    try:
        df = pd.read_excel(excel_file, dtype=str)

        lid_cols = df.filter(like="LID").columns
        barcode_cols = df.filter(like="BARCODE").columns

        if lid_cols.empty or barcode_cols.empty:
            logger.warning(f"Could not find LID or BARCODE columns in {excel_file}")
            return 0

        lid_col = lid_cols[0]
        barcode_col = barcode_cols[0]

        result_df = pd.DataFrame({"lid": df[lid_col], "barcode": df[barcode_col]})

        # Remove rows with missing LID or barcode
        result_df = result_df.dropna(subset=["lid", "barcode"])

        numeric_mask = result_df["barcode"].notna() & result_df["barcode"].str.match(r"^\d+$")
        if numeric_mask.any():
            barcode_len = int(result_df.loc[numeric_mask, "barcode"].str.len().max())
            result_df.loc[numeric_mask, "barcode"] = (
                result_df.loc[numeric_mask, "barcode"].astype(str).str.zfill(barcode_len)
            )

        result_df.to_csv(tsv_path, sep="\t", index=False, header=False)
        logger.success(f"Successfully exported data to {tsv_path}")
        logger.info(f"Total samples: {len(result_df)}")
        return len(result_df)

    except FileNotFoundError:
        logger.error(f"Input file not found: {excel_file}")
        return 0
    except pd.errors.EmptyDataError:
        logger.error(f"Input file is empty: {excel_file}")
        return 0
    except (OSError, pd.errors.ParserError, ValueError, KeyError) as e:
        logger.error(f"Error exporting TSV file: {e!s}")
        logger.exception("Traceback:")
        return 0


def _batch_tsv_path(excel_file: Path) -> Path:
    """Return the metadata TSV path matching an Excel workbook."""
    return excel_file.parent / f"{excel_file.stem}.tsv"


def _tsv_readable(tsv_path: Path) -> bool:
    """Check whether a metadata TSV can be opened for reading."""
    try:
        with tsv_path.open():
            pass
    except OSError:
        return False
    return True


def process_excel_files(excluded_batches: list[str], paths: MergePaths) -> int:
    """
    Convert metadata Excel workbooks to TSV.

    Only batches whose TSV is missing or unreadable are converted; existing
    readable TSVs are left untouched. Readability is checked concurrently
    because the metadata directory lives on a high-latency NAS share.

    Args:
        excluded_batches: List of batch names to exclude
        paths: Merge input/output path contract

    Returns:
        Number of samples exported
    """
    metadata_path = paths.metadata_dir
    excel_files = list(metadata_path.glob("*.xlsx"))
    logger.info(f"Found {len(excel_files)} Excel files in metadata directory.")

    excluded = set(excluded_batches)
    candidates = [excel_file for excel_file in excel_files if excel_file.stem not in excluded]

    total_exported = 0
    converted = 0
    if candidates:
        with ThreadPoolExecutor(max_workers=min(_METADATA_WORKERS, len(candidates))) as pool:
            readable_flags = list(pool.map(_tsv_readable, [_batch_tsv_path(f) for f in candidates]))
        stale = [excel_file for excel_file, readable in zip(candidates, readable_flags, strict=True) if not readable]

        for excel_file in stale:
            total_exported += convert_excel_to_tsv(excel_file, _batch_tsv_path(excel_file))
            converted += 1

    skipped = len(candidates) - converted
    if converted:
        logger.info(
            f"Converted {converted} Excel workbooks to TSV ({total_exported} samples); "
            f"skipped {skipped} batches with existing TSVs",
        )
    else:
        logger.info(f"No Excel conversion needed; {skipped} batches already have TSVs")
    return total_exported
