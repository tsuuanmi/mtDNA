#!/usr/bin/env python3
"""
ms_control_QC.py

Control QC module for MS mtDNA post-processing.

Purpose
-------
This script is an optional pipeline step after ms_trace_process.py has generated MS_trace_processed.xlsx.
It checks PC and NTC controls at a profile level using the already-normalized
`Variant Converted 2` values from the MS_trace_processed output.

Expected control naming in MS_trace_processed LID column
-----------------------------------------------
PC controls:
    PC1_mtDNA_247547
    PC2_mtDNA_247547

    Pattern: ^PC\\d+_(.+)$
    The part after PC<number>_ is treated as the reference/control source LID
    and is looked up in the control reference file.

NTC controls:
    NTC1
    NTC2

    Pattern: ^NTC\\d+$

Analysis ranges used for profile construction
---------------------------------------------
Only variants in the defined mtDNA analysis ranges are used:
    HV2-3: 73-340, 438-576
    HV1  : 16024-16365

Control reference file
----------------------
Excel file with one row per PC reference LID.
Accepted column names are intentionally flexible.

Reference key column, first found among:
    Reference_LID, LID, Control_LID, Control_ID, PC_LID

Expected profile column, first found among:
    Expected_Profile, Profile, Expected_Variants, Variants

Example:
    Reference_LID       Expected_Profile
    mtDNA_247547        73G 263G 309.1C 315.1C 16093C ...

Behavior
--------
- PC observed profile = sorted unique union of Variant Converted 2 from its 4 traces,
  filtered to the defined analysis ranges.
- PC passes only if observed profile exactly matches expected profile.
- NTC observed profile must be empty.
- Each control LID must have exactly one trace for each expected primer:
  HV1F, HV1R, HV2F, HV3R.
- If any control fails, the script exits with code 1.
- The script always writes an Excel QC report.
"""

import argparse
import os
import re
from collections.abc import Iterable, Sequence

import pandas as pd

EXPECTED_PRIMERS = ["HV1F", "HV1R", "HV2F", "HV3R"]
ANALYSIS_RANGES = [
    (73.0, 340.999999),
    (438.0, 576.999999),
    (16024.0, 16365.999999),
]

REF_KEY_COLUMNS = ["Reference_LID", "LID", "Control_LID", "Control_ID", "PC_LID"]
REF_PROFILE_COLUMNS = ["Expected_Profile", "Profile", "Expected_Variants", "Variants"]

PC_RE = re.compile(r"^PC\d+_(.+)$", re.IGNORECASE)
NTC_RE = re.compile(r"^NTC\d+$", re.IGNORECASE)
POS_RE = re.compile(r"^(\d+(?:\.\d+)?)(.+)$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Optional PC/NTC control QC for MS_trace_processed output")
    p.add_argument(
        "--trace-processed",
        required=False,
        help="Input {BATCH}_MS_trace_processed.xlsx from ms_trace_process.py",
    )
    p.add_argument("--compare-qc", required=False, help="Deprecated alias for --trace-processed")
    p.add_argument("--control-ref", required=True, help="Global PC reference profile Excel file")
    p.add_argument("--out", required=True, help="Output control QC report xlsx")
    p.add_argument(
        "--allow-missing-ntc",
        action="store_true",
        help="Do not fail if no NTC controls are found. Default: fail.",
    )
    p.add_argument(
        "--allow-missing-pc",
        action="store_true",
        help="Do not fail if no PC controls are found. Default: fail.",
    )
    args = p.parse_args()
    if not args.trace_processed and not args.compare_qc:
        p.error("one of --trace-processed or deprecated --compare-qc is required")
    if args.trace_processed and args.compare_qc:
        p.error("use only one of --trace-processed or --compare-qc")
    args.trace_processed = args.trace_processed or args.compare_qc
    return args


def clean_str(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def find_first_column(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    column_map = {str(c).strip().lower(): c for c in columns}
    for candidate in candidates:
        c = column_map.get(candidate.lower())
        if c is not None:
            return c
    return None


def token_position(token: str) -> float | None:
    token = clean_str(token)
    m = POS_RE.match(token)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def in_analysis_range(pos: float) -> bool:
    return any(start <= pos <= end for start, end in ANALYSIS_RANGES)


def variant_sort_key(token: str) -> tuple[float, str]:
    pos = token_position(token)
    if pos is None:
        return (10**12, token)
    return (pos, token)


def normalize_profile(profile_text: str, filter_ranges: bool = True) -> list[str]:
    """Convert profile text into sorted unique tokens, optionally filtered by analysis ranges."""
    tokens: set[str] = set()
    for raw in clean_str(profile_text).replace(";", " ").replace(",", " ").split():
        tok = clean_str(raw)
        if not tok:
            continue
        pos = token_position(tok)
        if filter_ranges:
            if pos is None or not in_analysis_range(pos):
                continue
        tokens.add(tok)
    return sorted(tokens, key=variant_sort_key)


def profile_to_string(tokens: Iterable[str]) -> str:
    return " ".join(sorted(set(tokens), key=variant_sort_key))


def read_reference(path: str) -> dict[str, list[str]]:
    try:
        ref = pd.read_excel(path)
    except Exception as exc:
        raise SystemExit(f"ERROR: Could not read control reference file: {path}\n{exc}")

    key_col = find_first_column(ref.columns, REF_KEY_COLUMNS)
    profile_col = find_first_column(ref.columns, REF_PROFILE_COLUMNS)

    if key_col is None or profile_col is None:
        raise SystemExit(
            "ERROR: Control reference file must contain one key column and one profile column.\n"
            f"Accepted key columns: {', '.join(REF_KEY_COLUMNS)}\n"
            f"Accepted profile columns: {', '.join(REF_PROFILE_COLUMNS)}"
        )

    reference: dict[str, list[str]] = {}
    duplicates: list[str] = []

    for _, row in ref.iterrows():
        key = clean_str(row.get(key_col))
        if not key:
            continue
        if key in reference:
            duplicates.append(key)
            continue
        reference[key] = normalize_profile(row.get(profile_col, ""), filter_ranges=True)

    if duplicates:
        dup_str = ", ".join(sorted(set(duplicates)))
        raise SystemExit(f"ERROR: Duplicate reference IDs in control reference file: {dup_str}")

    if not reference:
        raise SystemExit("ERROR: No usable PC reference profiles found in control reference file.")

    return reference


def detect_control_type(lid: str) -> tuple[str, str]:
    """Return (control_type, reference_lid). control_type is PC, NTC, or empty."""
    lid = clean_str(lid)
    m = PC_RE.match(lid)
    if m:
        return "PC", clean_str(m.group(1))
    if NTC_RE.match(lid):
        return "NTC", ""
    return "", ""


def check_required_columns(df: pd.DataFrame, path: str) -> None:
    required = ["LID", "Primer", "Variant Converted 2"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: Missing required columns in {path}: {', '.join(missing)}\n"
            "Expected input is MS_trace_processed output from ms_trace_process.py."
        )


def get_control_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    parsed = out["LID"].apply(lambda x: detect_control_type(clean_str(x)))
    out["Control_Type"] = parsed.apply(lambda x: x[0])
    out["Reference_LID"] = parsed.apply(lambda x: x[1])
    return out[out["Control_Type"].isin(["PC", "NTC"])].copy()


def check_primer_inventory(control_lid: str, sub: pd.DataFrame) -> tuple[bool, list[str], str]:
    issues: list[str] = []
    counts = sub["Primer"].astype(str).str.upper().value_counts().to_dict()

    for primer in EXPECTED_PRIMERS:
        n = int(counts.get(primer, 0))
        if n == 0:
            issues.append(f"Missing {primer}")
        elif n > 1:
            issues.append(f"Duplicate {primer} ({n} traces)")

    unexpected = sorted(p for p in counts.keys() if p and p not in EXPECTED_PRIMERS)
    for primer in unexpected:
        issues.append(f"Unexpected primer {primer} ({counts[primer]} traces)")

    status = not issues
    detail = "; ".join(issues)
    return status, issues, detail


def build_observed_profile(sub: pd.DataFrame) -> list[str]:
    tokens: set[str] = set()
    for value in sub["Variant Converted 2"].fillna("").astype(str):
        tokens.update(normalize_profile(value, filter_ranges=True))
    return sorted(tokens, key=variant_sort_key)


def compare_profiles(observed: list[str], expected: list[str]) -> tuple[list[str], list[str]]:
    observed_set = set(observed)
    expected_set = set(expected)
    missing = sorted(expected_set - observed_set, key=variant_sort_key)
    extra = sorted(observed_set - expected_set, key=variant_sort_key)
    return missing, extra


def make_report(
    df_controls: pd.DataFrame,
    reference: dict[str, list[str]],
    allow_missing_pc: bool,
    allow_missing_ntc: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, bool]:
    details: list[dict] = []
    profiles: list[dict] = []

    pc_lids = sorted(df_controls.loc[df_controls["Control_Type"] == "PC", "LID"].dropna().astype(str).unique())
    ntc_lids = sorted(df_controls.loc[df_controls["Control_Type"] == "NTC", "LID"].dropna().astype(str).unique())

    if not pc_lids and not allow_missing_pc:
        details.append(
            {
                "Control_LID": "",
                "Control_Type": "PC",
                "Check": "Presence",
                "Status": "FAIL",
                "Issue": "No PC controls found",
                "Reference_LID": "",
                "Observed_Profile": "",
                "Expected_Profile": "",
                "Missing_Variants": "",
                "Extra_Variants": "",
            }
        )

    if not ntc_lids and not allow_missing_ntc:
        details.append(
            {
                "Control_LID": "",
                "Control_Type": "NTC",
                "Check": "Presence",
                "Status": "FAIL",
                "Issue": "No NTC controls found",
                "Reference_LID": "",
                "Observed_Profile": "",
                "Expected_Profile": "",
                "Missing_Variants": "",
                "Extra_Variants": "",
            }
        )

    # Evaluate each control independently.
    for control_lid, sub in df_controls.groupby("LID", dropna=False):
        control_lid = clean_str(control_lid)
        control_type, reference_lid = detect_control_type(control_lid)
        inventory_ok, inventory_issues, inventory_detail = check_primer_inventory(control_lid, sub)
        observed = build_observed_profile(sub)
        observed_text = profile_to_string(observed)

        if control_type == "NTC":
            expected: list[str] = []
            missing, extra = [], observed
            profile_ok = len(observed) == 0
            status = "PASS" if inventory_ok and profile_ok else "FAIL"
            issue_parts = []
            if not inventory_ok:
                issue_parts.append(inventory_detail)
            if not profile_ok:
                issue_parts.append("NTC has variants in analysis range")

            details.append(
                {
                    "Control_LID": control_lid,
                    "Control_Type": "NTC",
                    "Check": "NTC_Profile",
                    "Status": status,
                    "Issue": "; ".join(issue_parts),
                    "Reference_LID": "",
                    "Observed_Profile": observed_text,
                    "Expected_Profile": "",
                    "Missing_Variants": "",
                    "Extra_Variants": profile_to_string(extra),
                }
            )
            profiles.append(
                {
                    "Control_LID": control_lid,
                    "Control_Type": "NTC",
                    "Reference_LID": "",
                    "Observed_Profile": observed_text,
                    "Expected_Profile": "",
                }
            )
            continue

        if control_type == "PC":
            expected = reference.get(reference_lid)
            if expected is None:
                status = "FAIL"
                issue_parts = [f"Reference profile not found for {reference_lid}"]
                if not inventory_ok:
                    issue_parts.insert(0, inventory_detail)
                details.append(
                    {
                        "Control_LID": control_lid,
                        "Control_Type": "PC",
                        "Check": "PC_Profile",
                        "Status": status,
                        "Issue": "; ".join(issue_parts),
                        "Reference_LID": reference_lid,
                        "Observed_Profile": observed_text,
                        "Expected_Profile": "",
                        "Missing_Variants": "",
                        "Extra_Variants": observed_text,
                    }
                )
                profiles.append(
                    {
                        "Control_LID": control_lid,
                        "Control_Type": "PC",
                        "Reference_LID": reference_lid,
                        "Observed_Profile": observed_text,
                        "Expected_Profile": "",
                    }
                )
                continue

            missing, extra = compare_profiles(observed, expected)
            profile_ok = len(missing) == 0 and len(extra) == 0
            status = "PASS" if inventory_ok and profile_ok else "FAIL"
            issue_parts = []
            if not inventory_ok:
                issue_parts.append(inventory_detail)
            if missing:
                issue_parts.append("Missing expected variants")
            if extra:
                issue_parts.append("Extra observed variants")

            details.append(
                {
                    "Control_LID": control_lid,
                    "Control_Type": "PC",
                    "Check": "PC_Profile",
                    "Status": status,
                    "Issue": "; ".join(issue_parts),
                    "Reference_LID": reference_lid,
                    "Observed_Profile": observed_text,
                    "Expected_Profile": profile_to_string(expected),
                    "Missing_Variants": profile_to_string(missing),
                    "Extra_Variants": profile_to_string(extra),
                }
            )
            profiles.append(
                {
                    "Control_LID": control_lid,
                    "Control_Type": "PC",
                    "Reference_LID": reference_lid,
                    "Observed_Profile": observed_text,
                    "Expected_Profile": profile_to_string(expected),
                }
            )

    # PC replicate concordance check: all PCs referencing the same Reference_LID should match.
    pc_profile_df = pd.DataFrame([p for p in profiles if p["Control_Type"] == "PC"])
    if not pc_profile_df.empty:
        for reference_lid, sub in pc_profile_df.groupby("Reference_LID", dropna=False):
            observed_profiles = sorted(set(sub["Observed_Profile"].astype(str)))
            if len(observed_profiles) > 1:
                details.append(
                    {
                        "Control_LID": "; ".join(sorted(sub["Control_LID"].astype(str))),
                        "Control_Type": "PC",
                        "Check": "PC_Replicate_Concordance",
                        "Status": "FAIL",
                        "Issue": "PC replicate profiles are discordant",
                        "Reference_LID": reference_lid,
                        "Observed_Profile": " || ".join(observed_profiles),
                        "Expected_Profile": "",
                        "Missing_Variants": "",
                        "Extra_Variants": "",
                    }
                )

    details_df = pd.DataFrame(details)
    profiles_df = pd.DataFrame(profiles)

    fail_count = int((details_df["Status"] == "FAIL").sum()) if not details_df.empty else 0
    pass_count = int((details_df["Status"] == "PASS").sum()) if not details_df.empty else 0
    overall_status = "FAIL" if fail_count > 0 else "PASS"

    summary_rows = [
        {"Metric": "Control_QC_Status", "Value": overall_status},
        {"Metric": "PC_controls_found", "Value": len(pc_lids)},
        {"Metric": "NTC_controls_found", "Value": len(ntc_lids)},
        {"Metric": "PASS_checks", "Value": pass_count},
        {"Metric": "FAIL_checks", "Value": fail_count},
        {"Metric": "Analysis_ranges", "Value": "73-340;438-576;16024-16365"},
        {"Metric": "Expected_primers", "Value": ",".join(EXPECTED_PRIMERS)},
    ]
    summary_df = pd.DataFrame(summary_rows)
    return summary_df, details_df, profiles_df, overall_status == "PASS"


def write_report(
    out_path: str,
    summary_df: pd.DataFrame,
    details_df: pd.DataFrame,
    profiles_df: pd.DataFrame,
    reference: dict[str, list[str]],
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    ref_df = pd.DataFrame(
        [
            {"Reference_LID": key, "Expected_Profile": profile_to_string(value)}
            for key, value in sorted(reference.items())
        ]
    )
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        details_df.to_excel(writer, sheet_name="Details", index=False)
        profiles_df.to_excel(writer, sheet_name="Observed_Profiles", index=False)
        ref_df.to_excel(writer, sheet_name="Reference_Profiles", index=False)


def main() -> None:
    args = parse_args()

    try:
        df = pd.read_excel(args.trace_processed)
    except Exception as exc:
        raise SystemExit(f"ERROR: Could not read MS_trace_processed file: {args.trace_processed}\n{exc}")

    check_required_columns(df, args.trace_processed)
    reference = read_reference(args.control_ref)
    df_controls = get_control_rows(df)

    summary_df, details_df, profiles_df, ok = make_report(
        df_controls=df_controls,
        reference=reference,
        allow_missing_pc=args.allow_missing_pc,
        allow_missing_ntc=args.allow_missing_ntc,
    )

    write_report(args.out, summary_df, details_df, profiles_df, reference)

    status = "PASS" if ok else "FAIL"
    print("====================================================")
    print("CONTROL QC - PC / NTC")
    print("====================================================")
    print(f"Input MS_trace_processed : {args.trace_processed}")
    print(f"Control reference: {args.control_ref}")
    print(f"Report           : {args.out}")
    print(f"Status           : {status}")
    if not ok and not details_df.empty:
        failed = details_df[details_df["Status"] == "FAIL"]
        print("Failed checks:")
        for _, row in failed.iterrows():
            ctl = clean_str(row.get("Control_LID")) or "<batch>"
            check = clean_str(row.get("Check"))
            issue = clean_str(row.get("Issue"))
            print(f"  - {ctl} | {check} | {issue}")
    print("====================================================")

    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
