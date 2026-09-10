"""Convert Mutation Surveyor ``*_custom_report.txt`` files to Excel.

The TXT export is tab-delimited.  Each row is copied to Excel in the same
position so metadata rows, trace headers, and trace records stay unchanged.

Amino acid notation (e.g. ``p.Leu104LeuPro``) in variant cells is stripped
during conversion so that only nucleotide-level information and quality
scores are retained in the output.

The output filename follows the pattern
``{BATCH}_{anything}_df_custom_report.xlsx`` regardless of the input
suffix (``_CR_mutationSurveyor_CRSetting``, ``_default_custom_report``,
``_df_custom_report``, or ``_custom_report``).

Public API:
    - ``convert_txt_to_excel()`` — convert a TXT custom report to XLSX.
    - ``read_custom_report_txt()`` — parse a TXT custom report into rows.
    - ``is_txt_custom_report()`` — check whether a path looks like an MS TXT report.
"""

import argparse
import csv
import re
from pathlib import Path

from loguru import logger
from openpyxl import Workbook

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_SHEET_NAME = "Custom Report"
DEFAULT_ENCODING = "utf-8-sig"

# Filename suffixes that identify a Mutation Surveyor custom report.
# Ordered from longest to shortest so that ``_normalize_output_stem``
# strips the most-specific suffix first (``_df_custom_report`` must be
# matched before the shorter ``_custom_report``).
_CUSTOM_REPORT_SUFFIXES = (
    "_CR_mutationSurveyor_CRSetting",
    "_default_custom_report",
    "_df_custom_report",
    "_custom_report",
)

# Strip amino-acid notation (e.g. ",p.Leu104LeuPro") from variant cells.
_AMINO_ACID_RE = re.compile(r",p\.[^,]+")

__all__ = [
    "convert_txt_to_excel",
    "is_txt_custom_report",
    "read_custom_report_txt",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_amino_acid(value: str) -> str:
    """Remove amino acid notation from a variant cell value.

    Strips protein-change tokens like ``p.Leu104LeuPro`` that appear
    after a comma in variant strings, e.g.
    ``c.310T>TC,p.Leu104LeuPro,$7`` → ``c.310T>TC,$7``.
    """
    return _AMINO_ACID_RE.sub("", value)


def _normalize_output_stem(stem: str) -> str:
    """Normalize a filename stem to end with ``_df_custom_report``.

    Strips any recognised custom-report suffix from *stem* and appends
    ``_df_custom_report``, ensuring the output always follows the
    ``{BATCH}_{anything}_df_custom_report`` naming convention.
    """
    for suffix in _CUSTOM_REPORT_SUFFIXES:
        if stem.endswith(suffix):
            return stem[: -len(suffix)] + "_df_custom_report"
    # Fallback: if no suffix matched, append _df_custom_report
    return stem + "_df_custom_report"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def is_txt_custom_report(path: Path) -> bool:
    """Return True if *path* looks like a Mutation Surveyor TXT custom report.

    Checks for ``.txt`` extension and a filename suffix matching one of
    ``_CR_mutationSurveyor_CRSetting``, ``_default_custom_report``,
    ``_df_custom_report``, or ``_custom_report``.
    """
    path = Path(path)
    if path.suffix.lower() != ".txt":
        return False
    stem = path.stem
    return any(stem.endswith(suffix) for suffix in _CUSTOM_REPORT_SUFFIXES)


def read_custom_report_txt(txt_path: Path, encoding: str = DEFAULT_ENCODING) -> list[list[str]]:
    """Read a Mutation Surveyor TXT report into rows.

    A valid report must include a trace-table header row containing both
    ``Trace #`` and ``Sample Name``.

    Args:
        txt_path: Path to the TXT custom report file.
        encoding: Text encoding (default: ``utf-8-sig``).

    Returns:
        List of rows, where each row is a list of string cells.

    Raises:
        FileNotFoundError: If *txt_path* does not exist.
        ValueError: If the file is empty or lacks a trace-table header.
    """
    txt_path = Path(txt_path)
    if not txt_path.exists():
        msg = f"TXT report not found: {txt_path}"
        raise FileNotFoundError(msg)

    with txt_path.open(encoding=encoding, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        msg = f"TXT report is empty: {txt_path}"
        raise ValueError(msg)
    if not any("Trace #" in row and "Sample Name" in row for row in rows):
        msg = (
            "Input does not look like a Mutation Surveyor custom report: "
            "no row contains both 'Trace #' and 'Sample Name'."
        )
        raise ValueError(msg)

    return rows


def convert_txt_to_excel(
    txt_path: Path,
    output_path: Path | None = None,
    *,
    sheet_name: str = DEFAULT_SHEET_NAME,
    encoding: str = DEFAULT_ENCODING,
) -> Path:
    """Convert a TXT custom report to an XLSX custom report.

    Amino acid notation (``p.XXX``) in variant cells is stripped so only
    nucleotide-level information and quality scores are retained.

    The output filename follows the pattern
    ``{BATCH}_{anything}_df_custom_report.xlsx``.  When *output_path* is
    not given, it is derived from the input filename by normalising the
    stem to end with ``_df_custom_report`` and changing the extension to
    ``.xlsx``.

    Args:
        txt_path: Path to the TXT custom report.
        output_path: Optional explicit output path for the XLSX file.
        sheet_name: Name for the Excel worksheet (max 31 chars).
        encoding: Text encoding of the input file.

    Returns:
        Path to the written XLSX file.
    """
    txt_path = Path(txt_path)

    if output_path is not None:
        xlsx_path = Path(output_path)
    else:
        normalized_stem = _normalize_output_stem(txt_path.stem)
        xlsx_path = txt_path.with_name(normalized_stem + ".xlsx")

    logger.info("Converting TXT → XLSX: {} → {}", txt_path, xlsx_path)

    rows = read_custom_report_txt(txt_path, encoding=encoding)

    # Strip amino acid notation from every cell
    rows = [[_strip_amino_acid(cell) for cell in row] for row in rows]

    workbook = Workbook()
    worksheet = workbook.active
    if worksheet is None:
        msg = "Workbook.active returned None unexpectedly"
        raise RuntimeError(msg)
    worksheet.title = sheet_name[:31]

    for row in rows:
        worksheet.append([value or None for value in row])

    xlsx_path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(xlsx_path)

    logger.success("XLSX written: {}", xlsx_path)
    return xlsx_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    """CLI entry point: convert a Mutation Surveyor custom-report TXT to XLSX."""
    parser = argparse.ArgumentParser(description="Convert a Mutation Surveyor custom-report TXT to XLSX")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        dest="txt_path",
        help="Input *_custom_report.txt file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output .xlsx path (default: normalised stem with .xlsx extension)",
    )
    args = parser.parse_args(argv)

    convert_txt_to_excel(
        txt_path=Path(args.txt_path),
        output_path=Path(args.output) if args.output else None,
    )


if __name__ == "__main__":
    main()
