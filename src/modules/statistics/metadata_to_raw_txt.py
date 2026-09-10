#!/usr/bin/env python
"""Create raw sample-ID TXT files from mtDNA metadata Excel workbooks."""

from __future__ import annotations

import argparse
import math
import re
from numbers import Real
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

import pandas as pd
from loguru import logger

DEFAULT_METADATA_DIR = Path("/mnt/nas/bca/mtDNA/science/data/metadata")
DEFAULT_RAW_DIR = Path("/mnt/nas/bca/mtDNA/science/data/raw")
SAMPLE_ID_COLUMNS = ("LID", "Sample ID", "Sample_ID", "SampleID")
EXCEL_SUFFIXES = {".xlsx", ".xls", ".xlsm"}
CONTROL_SAMPLE_PATTERN = "PC|NTC"


def normalize_header(value: object) -> str:
    """Normalize a spreadsheet header for case-insensitive matching."""
    return str(value).strip().lower().replace(" ", "").replace("_", "")


def format_sample_id(value: object) -> str | None:
    """Convert an Excel cell value to a raw sample ID string."""
    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, Real) and math.isnan(value):
        return None

    sample_id = str(value).strip()
    if not sample_id:
        return None

    return sample_id.removesuffix(".0") or None


def find_sample_id_column(columns: Iterable[object]) -> object | None:
    """Find the sample-ID column from a DataFrame's columns."""
    normalized_targets = {normalize_header(column) for column in SAMPLE_ID_COLUMNS}
    for column in columns:
        if normalize_header(column) in normalized_targets:
            return column
    return None


def is_control_sample(sample_id: str) -> bool:
    """Return whether a sample ID should be excluded as a control."""
    return re.search(CONTROL_SAMPLE_PATTERN, sample_id, flags=re.IGNORECASE) is not None


def extract_sample_ids(excel_path: Path) -> list[str]:
    """Extract unique non-control sample IDs from every sheet with a recognized ID column."""
    sample_ids: list[str] = []
    seen: set[str] = set()

    excel_file = pd.ExcelFile(excel_path)
    for sheet_name in excel_file.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        sample_id_column = find_sample_id_column(df.columns)
        if sample_id_column is None:
            logger.warning("{} [{}]: no sample ID column found", excel_path.name, sheet_name)
            continue

        for value in df[sample_id_column]:
            sample_id = format_sample_id(value)
            if sample_id is None or sample_id in seen or is_control_sample(sample_id):
                continue
            seen.add(sample_id)
            sample_ids.append(sample_id)

    return sample_ids


def iter_excel_files(metadata_dir: Path) -> Iterable[Path]:
    """Yield non-temporary Excel files in stable name order."""
    for path in sorted(metadata_dir.iterdir()):
        if path.name.startswith("~$"):
            continue
        if path.is_file() and path.suffix.lower() in EXCEL_SUFFIXES:
            yield path


def read_raw_sample_ids(output_path: Path) -> list[str] | None:
    """Read non-empty sample IDs from an existing raw TXT file."""
    if not output_path.exists():
        return None
    return [line.strip() for line in output_path.read_text().splitlines() if line.strip()]


def verify_output(output_path: Path, expected: list[str]) -> None:
    """Verify the written raw TXT exactly matches the extracted Excel sample IDs."""
    actual = read_raw_sample_ids(output_path)
    if actual != expected:
        message = f"Generated TXT does not match Excel LID samples: {output_path}"
        raise RuntimeError(message)


def write_raw_txt(excel_path: Path, raw_dir: Path, *, dry_run: bool) -> Path | None:
    """Write one raw TXT file for one metadata workbook when missing or mismatched."""
    output_path = raw_dir / f"{excel_path.stem}.txt"
    sample_ids = extract_sample_ids(excel_path)
    if not sample_ids:
        logger.warning("{}: no non-control sample IDs extracted; not writing TXT", excel_path.name)
        return None

    existing_sample_ids = read_raw_sample_ids(output_path)
    if existing_sample_ids == sample_ids:
        logger.info(
            "{}: OK, TXT matches Excel with {} non-control LID samples",
            output_path.name,
            len(sample_ids),
        )
        return None

    if existing_sample_ids is not None:
        logger.info(
            "{}: TXT sample count {} differs from Excel sample count {}; rewriting from Excel",
            output_path.name,
            len(existing_sample_ids),
            len(sample_ids),
        )
    else:
        logger.info("{}: TXT missing; writing {} non-control LID samples", output_path.name, len(sample_ids))

    if dry_run:
        logger.info("Would write {} sample IDs to {}", len(sample_ids), output_path)
        return output_path

    raw_dir.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(sample_ids) + "\n")
    verify_output(output_path, sample_ids)
    logger.success("Wrote and verified {} sample IDs in {}", len(sample_ids), output_path)
    return output_path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-m",
        "--metadata-dir",
        type=Path,
        default=DEFAULT_METADATA_DIR,
        help=f"Directory containing metadata Excel files (default: {DEFAULT_METADATA_DIR})",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Directory for generated raw TXT files (default: {DEFAULT_RAW_DIR})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show mismatched/missing TXT files without creating files",
    )
    return parser.parse_args()


def main() -> None:
    """Run the metadata-to-raw conversion."""
    args = parse_args()
    if not args.metadata_dir.is_dir():
        message = f"Metadata directory not found: {args.metadata_dir}"
        raise SystemExit(message)

    written = 0
    for excel_path in iter_excel_files(args.metadata_dir):
        if write_raw_txt(excel_path, args.output_dir, dry_run=args.dry_run) is not None:
            written += 1

    action = "Would write" if args.dry_run else "Wrote/updated"
    logger.info("{} {} raw TXT files", action, written)


if __name__ == "__main__":
    main()
