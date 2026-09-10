#!/usr/bin/env python3
"""
Add columns from STR_mtDNA_with_ids.tsv to ThongKeMauPhanTich.tsv
by mapping CCCD (in STR) with Identification Number (in ThongKe).
Then map mtDNA_barcode with Excel file to get mtDNA_batch and mtDNA_ID.
Finally, fill remaining gaps using DSTM_cleaned.tsv (barcode_xn and id_number joins).
"""

from pathlib import Path

import pandas as pd
from loguru import logger


def _clean_str(df: pd.DataFrame, cols: list[str]) -> None:
    """Replace float NaN / literal 'nan' with '' for given columns."""
    for col in cols:
        df[col] = df[col].astype(str).str.strip().replace("nan", "")


def main() -> None:  # noqa: PLR0915
    # File paths
    str_file = Path("data/STR_mtDNA_with_ids.tsv")
    thongke_file = Path("data/ThongKeMauPhanTich.tsv")
    excel_file = Path("data/02. Samples_Total_n_Rurun_From_mt28-CNVV (1).xlsx")
    dstm_file = Path("data/DSTM_cleaned.tsv")
    output_file = Path("data/ThongKeMauPhanTich_result.tsv")

    # Read sources
    df_str = pd.read_csv(str_file, sep="\t", dtype=str, low_memory=False, keep_default_na=False)
    df_thongke = pd.read_csv(thongke_file, sep="\t", dtype=str, keep_default_na=False)
    df_excel = pd.read_excel(excel_file, sheet_name="Total", dtype=str)
    df_excel = df_excel[["LID", "Barcode", "MẺ CHẠY"]].copy()

    df_dstm = pd.read_csv(dstm_file, sep="\t", dtype=str, low_memory=False, keep_default_na=False)
    df_dstm = df_dstm[["barcode_xn", "barcode_mtdna", "batch_str", "batch_mtdna", "id_number"]].copy()
    _clean_str(df_dstm, df_dstm.columns)  # type: ignore[reportArgumentType]

    logger.info(f"STR data shape: {df_str.shape}")
    logger.info(f"ThongKe data shape: {df_thongke.shape}")
    logger.info(f"Excel data shape: {df_excel.shape}")
    logger.info(f"DSTM data shape: {df_dstm.shape}")

    # ── Step 1: Map CCCD → STR columns ─────────────────────────────────────────
    df_str["CCCD"] = df_str["CCCD"].astype(str).str.strip()
    df_thongke["Identification Number"] = df_thongke["Identification Number"].astype(str).str.strip()

    str_cols = ["STR_barcode", "mtDNA_barcode", "STR_ID", "STR_batch"]
    df_str_clean = df_str[df_str["CCCD"] != ""].drop_duplicates(subset=["CCCD"], keep="first")  # type: ignore[reportCallIssue]
    logger.info(f"STR clean (unique CCCD): {len(df_str_clean)} rows")
    lookup = df_str_clean.set_index("CCCD")[str_cols].to_dict("index")  # type: ignore[reportCallIssue]

    for col in str_cols:
        df_thongke[col] = df_thongke["Identification Number"].map(lambda x, col=col: lookup.get(x, {}).get(col) or "")

    # ── Step 2: Map mtDNA_barcode → Excel → mtDNA_batch, mtDNA_ID ─────────────
    df_excel["Barcode"] = df_excel["Barcode"].astype(str).str.strip()  # type: ignore[reportAttributeAccessIssue]
    df_thongke["mtDNA_barcode"] = df_thongke["mtDNA_barcode"].astype(str).str.strip()

    df_excel_clean = df_excel.drop_duplicates(subset=["Barcode"], keep="first")  # type: ignore[reportCallIssue]
    logger.info(f"Excel clean (unique Barcode): {len(df_excel_clean)} rows")
    excel_lookup = df_excel_clean.set_index("Barcode")[["MẺ CHẠY", "LID"]].to_dict("index")  # type: ignore[reportCallIssue]

    df_thongke["mtDNA_batch"] = df_thongke["mtDNA_barcode"].map(lambda x: excel_lookup.get(x, {}).get("MẺ CHẠY") or "")
    df_thongke["mtDNA_ID"] = df_thongke["mtDNA_barcode"].map(lambda x: excel_lookup.get(x, {}).get("LID") or "")

    # ── Step 3: Fill gaps using DSTM_cleaned ──────────────────────────────────
    # STR_barcode and mtDNA_barcode are equivalent in the result, so use
    # barcode_xn (unified lab barcode) as the primary join key.  Falls back
    # to id_number for rows that barcode_xn cannot match.
    df_thongke["STR_barcode"] = df_thongke["STR_barcode"].astype(str).str.strip()
    df_thongke["mtDNA_barcode"] = df_thongke["mtDNA_barcode"].astype(str).str.strip()

    # --- Stage 1: barcode_xn join (primary) ---
    df_dstm_xn = df_dstm[df_dstm["barcode_xn"] != ""].drop_duplicates(subset=["barcode_xn"], keep="first")  # type: ignore[reportCallIssue,reportAttributeAccessIssue]
    xn_merge = (
        df_thongke[["STR_barcode"]]
        .merge(
            df_dstm_xn[["barcode_xn", "batch_str", "batch_mtdna"]],
            left_on="STR_barcode",
            right_on="barcode_xn",
            how="left",
        )
        .drop(columns=["STR_barcode"])
    )
    xn_merge = xn_merge.fillna("")

    # --- Stage 2: id_number fallback ---
    df_dstm_id = df_dstm[df_dstm["id_number"] != ""].drop_duplicates(subset=["id_number"], keep="first")  # type: ignore[reportCallIssue,reportAttributeAccessIssue]
    id_merge = (
        df_thongke[["Identification Number"]]
        .merge(
            df_dstm_id[["id_number", "barcode_xn", "batch_str", "batch_mtdna"]],
            left_on="Identification Number",
            right_on="id_number",
            how="left",
        )
        .drop(columns=["id_number", "Identification Number"])
    )
    id_merge = id_merge.fillna("")

    # Fill blanks: existing > barcode_xn DSTM > id_number DSTM
    # STR_barcode and mtDNA_barcode are equivalent — both filled from barcode_xn
    for thongke_col, xn_col, id_col in [
        ("STR_barcode", "barcode_xn", "barcode_xn"),
        ("mtDNA_barcode", "barcode_xn", "barcode_xn"),
        ("STR_batch", "batch_str", "batch_str"),
        ("mtDNA_batch", "batch_mtdna", "batch_mtdna"),
    ]:
        existing = df_thongke[thongke_col].astype(str).str.strip()
        from_xn = xn_merge[xn_col]
        from_id = id_merge[id_col]
        df_thongke[thongke_col] = existing.where(existing != "", from_xn)
        df_thongke[thongke_col] = df_thongke[thongke_col].where(df_thongke[thongke_col] != "", from_id)

    # ── Statistics ────────────────────────────────────────────────────────────
    total = len(df_thongke)
    matched_cccd = (df_thongke["STR_barcode"] != "").sum()
    filled_str_batch = (df_thongke["STR_batch"] != "").sum()
    filled_mtdna_barcode = (df_thongke["mtDNA_barcode"] != "").sum()
    has_real_batch = (df_thongke["mtDNA_batch"] != "").sum()
    no_barcode = total - filled_mtdna_barcode

    logger.info("\n--- Fill Summary ---")
    logger.info(f"STR_barcode     : {matched_cccd}/{total} ({matched_cccd / total * 100:.1f}%)")
    logger.info(
        f"mtDNA_barcode   : {filled_mtdna_barcode}/{total} "
        f"({filled_mtdna_barcode / total * 100:.1f}%)  | {no_barcode} no barcode",
    )
    logger.info(f"STR_batch       : {filled_str_batch}/{total} ({filled_str_batch / total * 100:.1f}%)")
    logger.info(f"mtDNA_batch     : {has_real_batch}/{total} ({has_real_batch / total * 100:.1f}%)")

    batch_counts = df_thongke["mtDNA_batch"].value_counts()
    batch_counts = batch_counts[batch_counts.index != ""].head(10)  # type: ignore[reportAttributeAccessIssue]
    logger.info("\n--- Top 10 mtDNA_batch ---")
    for batch, count in batch_counts.items():
        logger.info(f"  {batch:<30s} {count}  ({count / total * 100:.1f}%)")
    no_batch = total - has_real_batch
    if no_batch:
        logger.info(f"  {'(no batch)':<30s} {no_batch}  ({no_batch / total * 100:.1f}%)")

    # ── Save ─────────────────────────────────────────────────────────────────
    final_cols = [
        "Identification Number",
        "Has STR",
        "Has mtDNA",
        "status",
        "address",
        "STR_barcode",
        "mtDNA_barcode",
        "STR_ID",
        "mtDNA_ID",
        "STR_batch",
        "mtDNA_batch",
    ]
    df_thongke = df_thongke[final_cols]
    df_thongke.to_csv(output_file, sep="\t", index=False)
    logger.info(f"\nResult saved to: {output_file}")
    logger.info(f"Result columns: {df_thongke.columns.tolist()}")


if __name__ == "__main__":
    main()
