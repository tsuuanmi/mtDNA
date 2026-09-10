#!/usr/bin/env python3
"""
ms_final_profiles_merge.py

Merge ms_hv23_merge.py and ms_hv1_merge.py LID sheets into one final
MS mtDNA profile workbook.

Input workbooks:
    - *_HV23_merged.xlsx, sheet: LID
    - *_HV1_merged.xlsx, sheet: LID

Output workbook:
    - *_final_profiles.xlsx

This script does not use a truthset. Truthset concordance is performed later by
ms_final_profiles_vs_truth.py.
"""

import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Merge HV2-3 and HV1 LID-level MS workbooks into final profiles")
    p.add_argument("--hv23", required=True, help="Input *_HV23_merged.xlsx from ms_hv23_merge.py")
    p.add_argument("--hv1", required=True, help="Input *_HV1_merged.xlsx from ms_hv1_merge.py")
    p.add_argument("--output", required=True, help="Output *_final_profiles.xlsx workbook")
    return p.parse_args()


def safe_str(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def token_sort_key(tok):
    m = re.match(r"^(\d+(?:\.\d+)?)", safe_str(tok))
    if not m:
        return (10**9, tok)
    return (float(m.group(1)), tok)


def merge_ranges(hv23_range, hv1_range):
    vals = []
    if safe_str(hv23_range):
        vals.append(safe_str(hv23_range))
    if safe_str(hv1_range):
        vals.append(safe_str(hv1_range))
    return " ".join(vals)


def merge_variants(hv23_vars, hv1_vars):
    toks = []

    for block in [hv23_vars, hv1_vars]:
        if safe_str(block):
            toks.extend(safe_str(block).split())

    toks = sorted(set(toks), key=token_sort_key)
    return " ".join(toks)


def combined_review(row):
    hv23_review = safe_str(row.get("HV2-3 Review Class", ""))
    hv1_review = safe_str(row.get("HV1 Review Class", ""))

    if hv23_review == "Review" or hv1_review == "Review":
        return "Review"

    return "Auto - no flag"


def combined_flags(row):
    vals = []

    for col in [
        "HV2-3 Flag list",
        "HV1 Review-triggering Flag list",
    ]:
        s = safe_str(row.get(col, ""))
        if s:
            vals.append(s)

    return "; ".join(vals)


def pct(v, total):
    if total == 0:
        return ""
    return round((100.0 * v / total), 2)


def build_summary(df):
    rows = []

    total = len(df)

    rows.append({"Metric": "Total LID", "Value": total, "%": 100.0 if total else ""})

    def add_counter(prefix, series):
        vc = series.fillna("").astype(str).value_counts(dropna=False)

        for k, v in vc.items():
            label = k if safe_str(k) else "(blank)"
            rows.append({"Metric": f"{prefix} - {label}", "Value": int(v), "%": pct(v, total)})

    add_counter("Combined Review", df["Combined Review Class"])

    if "HV2-3 Review Class" in df.columns:
        add_counter("HV2-3 Review", df["HV2-3 Review Class"])

    if "HV1 Review Class" in df.columns:
        add_counter("HV1 Review", df["HV1 Review Class"])

    flag_counter = Counter()

    for val in df["Combined Flag list"]:
        s = safe_str(val)
        if not s:
            continue

        for f in s.split("; "):
            if safe_str(f):
                flag_counter[f] += 1

    for k, v in sorted(flag_counter.items(), key=lambda x: (-x[1], x[0])):
        rows.append({"Metric": f"Flag - {k}", "Value": v, "%": pct(v, total)})

    return pd.DataFrame(rows)


def main():
    args = parse_args()

    hv23 = pd.read_excel(args.hv23, sheet_name="LID")
    hv1 = pd.read_excel(args.hv1, sheet_name="LID")

    hv23.columns = [str(c).strip() for c in hv23.columns]
    hv1.columns = [str(c).strip() for c in hv1.columns]

    # Backward/forward compatibility: earlier notes/scripts sometimes used
    # an underscore in this column name, while ms_hv23_merge.py writes the
    # hyphenated form. Normalize to the hyphenated column used in the current
    # pipeline.
    if "HV3R PHP 341_437 Sites" in hv23.columns and "HV3R PHP 341-437 Sites" not in hv23.columns:
        hv23 = hv23.rename(columns={"HV3R PHP 341_437 Sites": "HV3R PHP 341-437 Sites"})

    hv23_keep = [
        "Batch",
        "Well",
        "LID",
        "Range by HV2-3 Regions Adjusted",
        "Variants by HV2-3 Regions Adjusted",
        "HV2-3 Pair Status",
        "HV3R PHP 341-437 Sites",
        "HV2-3 Flag list",
        "HV2-3 Review Class",
    ]

    hv1_keep = [
        "Batch",
        "Well",
        "LID",
        "HV1 Pair Status",
        "Range by HV1 Region Adjusted",
        "Variants by HV1 Region Adjusted",
        "Mismatch by HV1 Region Adjusted",
        "HV1 Flag list",
        "HV1 Review-triggering Flag list",
        "HV1 Review Class",
    ]

    hv23_keep = [c for c in hv23_keep if c in hv23.columns]
    hv1_keep = [c for c in hv1_keep if c in hv1.columns]

    hv23 = hv23[hv23_keep].copy()
    hv1 = hv1[hv1_keep].copy()

    merged = pd.merge(hv23, hv1, on=["Batch", "LID"], how="outer", suffixes=("_HV23", "_HV1"))

    merged["Well"] = merged["Well_HV23"].fillna(merged["Well_HV1"])

    merged["Sample profile - Range"] = merged.apply(
        lambda r: merge_ranges(r.get("Range by HV2-3 Regions Adjusted", ""), r.get("Range by HV1 Region Adjusted", "")),
        axis=1,
    )

    merged["Sample profile - Variants"] = merged.apply(
        lambda r: merge_variants(
            r.get("Variants by HV2-3 Regions Adjusted", ""),
            r.get("Variants by HV1 Region Adjusted", ""),
        ),
        axis=1,
    )

    merged["Combined Flag list"] = merged.apply(combined_flags, axis=1)
    merged["Combined Review Class"] = merged.apply(combined_review, axis=1)

    final_cols = [
        "Batch",
        "Well",
        "LID",
        "Sample profile - Range",
        "Sample profile - Variants",
        "Range by HV2-3 Regions Adjusted",
        "Variants by HV2-3 Regions Adjusted",
        "HV2-3 Pair Status",
        "HV3R PHP 341-437 Sites",
        "HV2-3 Flag list",
        "HV2-3 Review Class",
        "HV1 Pair Status",
        "Range by HV1 Region Adjusted",
        "Variants by HV1 Region Adjusted",
        "Mismatch by HV1 Region Adjusted",
        "HV1 Flag list",
        "HV1 Review-triggering Flag list",
        "HV1 Review Class",
        "Combined Flag list",
        "Combined Review Class",
    ]

    final_cols = [c for c in final_cols if c in merged.columns]

    merged = merged[final_cols].copy()

    no_pc_ntc = merged[~merged["LID"].astype(str).str.upper().str.startswith(("PC", "NTC"))].copy()

    summary_all = build_summary(merged)
    summary_no_pc_ntc = build_summary(no_pc_ntc)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        no_pc_ntc.to_excel(writer, sheet_name="LID_no_PC_NTC", index=False)
        summary_no_pc_ntc.to_excel(writer, sheet_name="Summary_no_PC_NTC", index=False)

        merged.to_excel(writer, sheet_name="LID", index=False)
        summary_all.to_excel(writer, sheet_name="Summary", index=False)

        for sheet_name, ws in writer.sheets.items():
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions

            for col in ws.columns:
                max_len = 0
                letter = col[0].column_letter

                for cell in col:
                    value = "" if cell.value is None else str(cell.value)
                    max_len = max(max_len, len(value))

                ws.column_dimensions[letter].width = min(max_len + 2, 50)

    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
