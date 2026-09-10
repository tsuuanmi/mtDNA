#!/usr/bin/env python3
import re
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from loguru import logger

MAX_SINGLE_LETTER_COLS = 26  # Excel columns A-Z
_NUM_COMPARISON_SIDES = 2  # pipeline/automate + sequencher/manual
_COMPARISON_LABELS = ("Pipeline", "Sequencher")
_VARIANT_COLUMNS = [
    "Variants (Pipeline)",
    "Variants (Sequencher)",
    "Variants Unique (Pipeline)",
    "Variants Unique (Sequencher)",
]
_SUFFIX_RE = re.compile(r"\(([^)]+)\)\s*$")


def arg_parser() -> tuple[Path, Path, list[str] | None]:
    parser = ArgumentParser(description="Merge multiple Excel files containing comparison results")
    parser.add_argument("-i", "--input_dir", type=Path, required=True)
    parser.add_argument("-o", "--output_file", type=Path, required=True)
    parser.add_argument(
        "-b",
        "--batches",
        type=str,
        default=None,
        help="Comma-separated batch IDs to include (defaults to all per-batch "
        "files found under input_dir). Use this to merge only the batches "
        "processed in the current run instead of every batch subfolder.",
    )
    args = parser.parse_args()
    batches: list[str] | None = None
    if args.batches:
        batches = [b.strip() for b in args.batches.split(",") if b.strip()]
    return args.input_dir, args.output_file, batches


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize comparison columns to stable ``Pipeline``/``Sequencher`` labels.

    Comparison outputs may name columns by the source tool (e.g. ``Variants
    (tracy)``, ``Variants (sequencher)``) or by the older generic ``A``/``B``
    labels. Different batches may use different automate tools (tracy vs
    blastn), so column names vary and ``pd.concat`` would create disjoint
    columns full of NaN. The first tool in column order is always the
    automate/pipeline side and the second the sequencher/manual side, so map by
    first-appearance order while keeping merged output headers compatible with
    the legacy ``Pipeline``/``Sequencher`` layout.
    """
    cols = list(df.columns)
    suffixes: list[str] = []
    for c in cols:
        m = _SUFFIX_RE.search(c)
        if m and m.group(1) not in suffixes:
            suffixes.append(m.group(1))
    if len(suffixes) < _NUM_COMPARISON_SIDES:
        return df
    mapping = {suffixes[i]: _COMPARISON_LABELS[i] for i in range(_NUM_COMPARISON_SIDES)}
    rename: dict[str, str] = {}
    for c in cols:
        m = _SUFFIX_RE.search(c)
        if m and m.group(1) in mapping:
            rename[c] = _SUFFIX_RE.sub(f"({mapping[m.group(1)]})", c)
    return df.rename(columns=rename)


def _load(input_dir: Path, batches: list[str] | None = None) -> pd.DataFrame:
    files = [
        str(p)
        for p in input_dir.rglob("*.xlsx")
        if not p.name.startswith("~$")
        and not p.name.startswith("all_batches")
        and (batches is None or p.parent.name in batches)
    ]
    if not files:
        scope = f"batches {batches}" if batches else "all batches"
        msg = f"No .xlsx files found in {input_dir} ({scope})"
        raise FileNotFoundError(msg)

    dfs = []
    for f in files:
        try:
            dfs.append(_normalize_columns(pd.read_excel(f)))
        except (OSError, pd.errors.ParserError) as e:
            logger.error(f"Error reading {f}: {e}")

    if not dfs:
        msg = "No valid Excel files to merge"
        raise ValueError(msg)

    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Merged {len(dfs)} Excel files")
    return df


def _normalize_separator(df: pd.DataFrame, separator: str) -> pd.DataFrame:
    sep = " " if separator == "space" else ", "
    for col in _VARIANT_COLUMNS:
        if col in df.columns:
            df[col] = (
                df[col]
                .astype(str)
                .str.replace(r",\s*", " ", regex=True)
                .str.strip()
                .str.replace(r"\s+", sep, regex=True)
            )
    return df


def _write_excel(df: pd.DataFrame, path: Path) -> None:
    col_widths = {}
    for col in df.columns:
        max_len = df[col].astype(str).map(lambda x: len(str(x))).max()
        col_widths[col] = min(max_len + 2, 50)

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="All Batches", index=False)
        ws = writer.sheets["All Batches"]
        for col, width in col_widths.items():
            ws.column_dimensions[col[0] if len(df.columns) <= MAX_SINGLE_LETTER_COLS else "A"].width = width

    logger.success(f"Wrote {path}")


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)
    logger.success(f"Wrote {path}")


def merge_excel_files(input_dir: str | Path, output_file: str | Path, batches: list[str] | None = None) -> None:
    """Merge all per-batch comparison .xlsx files under input_dir.

    Each batch's tool-parameterized column names are normalized to the legacy
    ``Pipeline``/``Sequencher`` labels so batches align when concatenated.
    Writes both ``.xlsx`` and ``.tsv`` to the output path (with appropriate
    suffix).

    When ``batches`` is given, only the per-batch subfolders whose name matches
    one of the listed batch IDs are merged; ``None`` (default) merges every
    per-batch file found under ``input_dir``.
    """
    output_file = Path(output_file)

    merged = _load(Path(input_dir), batches)
    merged["Sample ID"] = merged["Sample ID"].astype(str)
    merged = merged.sort_values(by=["Batch", "Sample ID"]).reset_index(drop=True)
    merged = _normalize_separator(merged, "space")

    _write_excel(merged, output_file.with_suffix(".xlsx"))
    _write_tsv(merged, output_file.with_suffix(".tsv"))


def main() -> None:
    input_dir, output_file, batches = arg_parser()
    logger.info(f"Starting to merge comparison files from: {input_dir}")
    if batches:
        logger.info(f"Filtering to batches: {', '.join(batches)}")
    merge_excel_files(input_dir, output_file, batches)


if __name__ == "__main__":
    main()
