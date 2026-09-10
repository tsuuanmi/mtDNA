#!/usr/bin/env python3
"""
merge_files.py

Generic Excel merge utility for the MS post-processing pipeline.

It accepts one or more input files or glob patterns and concatenates all sheets'
first worksheets row-wise. It adds two helper columns:
  - Batch
  - Source_File

Batch is inferred from common pipeline suffixes, especially
`_MS_trace_processed.xlsx`.
"""

import argparse
import glob
import os
import re

import pandas as pd

KNOWN_SUFFIXES = [
    "_MS_trace_processed.xlsx",
    "_control_QC.xlsx",
    "_HV23_merged.xlsx",
    "_HV1_merged.xlsx",
    "_final_profiles.xlsx",
    "_final_profiles_vs_truth.xlsx",
    # Backward compatibility for older files.
    "_compare_QC.xlsx",
    "_df_compare_QC.xlsx",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Merge Excel files and add Batch/Source_File columns")
    parser.add_argument("--input", nargs="+", required=True, help="Input files or glob patterns")
    parser.add_argument("--output", required=True, help="Output Excel file path")
    parser.add_argument(
        "--batch-suffix",
        default=None,
        help="Optional suffix to remove from filenames when deriving Batch, e.g. _MS_trace_processed.xlsx",
    )
    return parser.parse_args()


def expand_inputs(patterns):
    files = []
    missing_patterns = []
    for pattern in patterns:
        expanded = glob.glob(pattern)
        if expanded:
            files.extend(expanded)
        elif os.path.exists(pattern):
            files.append(pattern)
        else:
            missing_patterns.append(pattern)
    return sorted(set(files)), missing_patterns


def infer_batch(filename, batch_suffix=None):
    name = os.path.basename(filename)
    if batch_suffix and name.endswith(batch_suffix):
        return name[: -len(batch_suffix)]
    for suffix in KNOWN_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return re.sub(r"\.xlsx$", "", name, flags=re.IGNORECASE)


def main():
    args = parse_args()
    files, missing_patterns = expand_inputs(args.input)

    if missing_patterns:
        print("WARNING: These input patterns/files did not match existing files:")
        for pattern in missing_patterns:
            print(f"  - {pattern}")

    if not files:
        print("NO INPUT FILE FOUND")
        raise SystemExit(1)

    dfs = []
    for f in files:
        df = pd.read_excel(f)
        filename = os.path.basename(f)
        df["Batch"] = infer_batch(filename, args.batch_suffix)
        df["Source_File"] = filename
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    combined.to_excel(args.output, index=False)

    print(f"DONE -> {args.output} ({len(files)} files)")


if __name__ == "__main__":
    main()
