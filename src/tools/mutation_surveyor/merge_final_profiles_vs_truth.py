#!/usr/bin/env python3
"""
merge_final_profiles_vs_truth.py

Merge *_final_profiles_vs_truth.xlsx from multiple batches into:
  1) all_batches_final_profiles_vs_truth.xlsx — all rows with Batch column filled
  2) HV1_AutoNoFlag_Y1F.xlsx  — LID_no_PC_NTC where
     HV1 Review Class == "Auto - no flag" AND HV1 Concordance in (Y1, F)
  3) HV23_AutoNoFlag_Y1F.xlsx — LID_no_PC_NTC where
     HV2-3 Review Class == "Auto - no flag" AND HV2-3 Concordance in (Y1, F)

Each input workbook has sheets: LID_no_PC_NTC, Summary, Discordance_detail.
The Batch column in LID_no_PC_NTC and Discordance_detail is None in individual
files; this script fills it from the batch directory name.
The Summary sheet gets a Batch column prepended.
"""

import argparse
import logging
from pathlib import Path
from typing import cast

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

logger = logging.getLogger(__name__)

BATCHES = [
    "MS_080426_002",
    "MS_080426_005",
    "MS_090426_003",
    "MS_100426_004",
    "MS_140426_001",
    "MS_170426_005",
    "MS_170426_006",
]

BASE_DIR = Path("/mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/mutation_surveyor")


def read_batch(batch: str, base_dir: Path) -> dict[str, pd.DataFrame]:
    """Read the three sheets from a batch's final_profiles_vs_truth workbook."""
    fpath = base_dir / batch / f"{batch}_final_profiles_vs_truth.xlsx"
    msg = f"Missing: {fpath}"
    if not fpath.exists():
        raise FileNotFoundError(msg)

    lid = pd.read_excel(fpath, sheet_name="LID_no_PC_NTC")
    summary = pd.read_excel(fpath, sheet_name="Summary")
    detail = pd.read_excel(fpath, sheet_name="Discordance_detail")

    lid["Batch"] = batch
    detail.insert(0, "Batch", batch)
    summary.insert(0, "Batch", batch)

    return {"LID_no_PC_NTC": lid, "Summary": summary, "Discordance_detail": detail}


def autosize_and_style(xlsx_path: Path) -> None:
    """Apply header styling, freeze panes, auto-filter, and auto-fit columns."""
    wb = load_workbook(xlsx_path)
    header_fill = PatternFill("solid", fgColor="D9EAF7")
    bold = Font(bold=True)

    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        for cell in ws[1]:
            cell.font = bold
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

        for col_cells in ws.columns:
            max_len = 0
            col_idx = col_cells[0].column
            if col_idx is None:
                continue  # skip columns without a valid index
            col_letter = get_column_letter(col_idx)
            for cell in col_cells[:2000]:
                val = "" if cell.value is None else str(cell.value)
                max_len = max(max_len, len(val))
            ws.column_dimensions[col_letter].width = min(max(max_len + 2, 10), 60)

    wb.save(xlsx_path)


def write_output(df_dict: dict[str, pd.DataFrame], output_path: Path) -> None:
    """Write a dict of {sheet_name: DataFrame} to an xlsx with styling."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name, df in df_dict.items():
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    autosize_and_style(output_path)
    logger.info("Written: %s (%d sheets)", output_path, len(df_dict))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Merge per-batch final_profiles_vs_truth workbooks into combined and filtered outputs."
    )
    parser.add_argument(
        "--base-dir",
        default=str(BASE_DIR),
        help=f"Base directory containing batch subdirectories (default: {BASE_DIR})",
    )
    parser.add_argument(
        "--batches",
        nargs="+",
        default=BATCHES,
        help=f"Batch IDs to merge (default: {BATCHES})",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path.home() / "workspaces" / "mtdna_raw" / "temp"),
        help="Output directory for merged files",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir)

    # ── Read and concatenate all batches ──
    all_lid: list[pd.DataFrame] = []
    all_summary: list[pd.DataFrame] = []
    all_detail: list[pd.DataFrame] = []

    for batch in args.batches:
        logger.info("Reading: %s", batch)
        sheets = read_batch(batch, base_dir)
        all_lid.append(sheets["LID_no_PC_NTC"])
        all_summary.append(sheets["Summary"])
        all_detail.append(sheets["Discordance_detail"])

    merged_lid = pd.concat(all_lid, ignore_index=True)
    merged_summary = pd.concat(all_summary, ignore_index=True)
    merged_detail = pd.concat(all_detail, ignore_index=True)

    # ── Output 1: all_batches_final_profiles_vs_truth.xlsx ──
    all_out = output_dir / "all_batches_final_profiles_vs_truth.xlsx"
    write_output(
        {
            "LID_no_PC_NTC": merged_lid,
            "Summary": merged_summary,
            "Discordance_detail": merged_detail,
        },
        all_out,
    )

    # ── Output 2: HV1_AutoNoFlag_Y1F.xlsx ──
    # Intersection: HV1 Review Class == "Auto - no flag" AND HV1 Concordance in (Y1, F)
    hv1_mask = (merged_lid["HV1 Review Class"] == "Auto - no flag") & merged_lid["HV1 Concordance"].isin(["Y1", "F"])
    hv1_filtered = cast("pd.DataFrame", merged_lid.loc[hv1_mask].copy())
    hv1_out = output_dir / "HV1_AutoNoFlag_Y1F.xlsx"
    write_output({"LID_no_PC_NTC": hv1_filtered}, hv1_out)

    # ── Output 3: HV23_AutoNoFlag_Y1F.xlsx ──
    # Intersection: HV2-3 Review Class == "Auto - no flag" AND HV2-3 Concordance in (Y1, F)
    hv23_mask = (merged_lid["HV2-3 Review Class"] == "Auto - no flag") & merged_lid["HV2-3 Concordance"].isin(
        ["Y1", "F"]
    )
    hv23_filtered = cast("pd.DataFrame", merged_lid.loc[hv23_mask].copy())
    hv23_out = output_dir / "HV23_AutoNoFlag_Y1F.xlsx"
    write_output({"LID_no_PC_NTC": hv23_filtered}, hv23_out)

    logger.info("Done. Merged %d batches.", len(args.batches))
    logger.info("  All rows:        %d", len(merged_lid))
    logger.info("  HV1 AutoNoFlag∩Y1F: %d", len(hv1_filtered))
    logger.info("  HV23 AutoNoFlag∩Y1F: %d", len(hv23_filtered))


if __name__ == "__main__":
    main()
