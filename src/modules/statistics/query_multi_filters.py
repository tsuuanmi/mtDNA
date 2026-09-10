#!/usr/bin/env python
"""
Run variant filter queries against a merged_statistics.json file.

For each query:
  - A sample MUST contain ALL required variants.
  - A sample must have NO OTHER variants within the specified exclusion window(s).
  - Exclusion windows may span multiple ranges; their union is checked.
  - The 'Variants' column shows all variants across all regions.

Each QUERIES entry: (output_filename, required_variants_set, exclusion_windows)
  - required_variants_set : set of (int_pos, seq_str)
  - exclusion_windows     : list of (start, end) — union of all windows is checked
"""

import csv
import json
from pathlib import Path
from typing import Any

from loguru import logger

from src.core.variants import pos_base, pos_sort_key

INPUT_JSON = "/mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json"

QUERIES = [
    (
        "NBBA0129.tsv",
        {(16136, "C"), (16183, "C"), (16189, "C"), (16260, "T"), (16287, "T"), (16325, "C")},
        [(16042, 16193), (16240, 16365)],
    ),
    (
        "NBAA0316.tsv",
        {(16140, "C"), (16187, "T"), (16189, "C"), (16256, "T"), (16266, "G")},
        [(16024, 16365)],
    ),
    (
        "NBAB0428.tsv",
        {(16140, "C"), (16182, "C"), (16183, "C"), (16189, "C"), (16261, "T"), (16266, "A")},
        [(16038, 16193), (16242, 16365)],
    ),
    (
        "NBBC1011.tsv",
        {(16223, "T"), (16287, "T"), (16319, "A"), (16362, "C")},
        [(16038, 16365)],
    ),
    (
        "NBAB0438.tsv",
        {(16051, "G"), (16189, "C"), (16269, "G"), (16292, "T"), (16300, "G"), (16304, "C")},
        [(16043, 16193), (16238, 16365)],
    ),
    (
        "NBBA0418.tsv",
        {(16108, "T"), (16129, "A"), (16162, "G"), (16172, "C"), (16304, "C"), (16357, "C")},
        [(16024, 16365)],
    ),
    (
        "NBBC1026.tsv",
        {(16086, "C"), (16129, "A"), (16209, "C"), (16223, "T"), (16260, "T"), (16272, "G")},
        [(16024, 16365)],
    ),
]

FIELDS = ["Sample_ID", "Batch", "Analyzed_Range", "Variants"]


def variants_in_windows(sample: dict[str, Any], windows: list[tuple[int, int]]) -> set[tuple[Any, str]]:
    """Return (pos, seq) pairs for all variants whose position falls in any window."""
    result: set[tuple[Any, str]] = set()
    for vtype in ("snps", "deletions", "insertions"):
        for v in sample.get("variants", {}).get(vtype, []):
            pos = v.get("pos")
            if pos is not None and any(s <= pos_base(pos) <= e for s, e in windows):
                result.add((pos_base(pos), v.get("seq", "")))
    return result


def format_analyzed_range(intervals_dict: dict[str, Any]) -> str:
    flat = sorted(
        (iv for ivs in intervals_dict.values() for iv in ivs),
        key=lambda x: x[0],
    )
    return ", ".join(f"{s}-{e}" for s, e in flat) or "None"


def format_variants(sample: dict[str, Any]) -> str:
    tokens = []
    for vtype in ("snps", "insertions", "deletions"):
        for v in sample.get("variants", {}).get(vtype, []):
            pos, seq = v.get("pos", ""), v.get("seq", "")
            token = f"{pos}DEL" if seq == "-" else f"{pos}{seq}"
            tokens.append((str(pos), token))
    tokens.sort(key=lambda x: pos_sort_key(x[0]))
    return " ".join(t for _, t in tokens) or "None"


logger.info("Loading merged_statistics.json …")
with Path(INPUT_JSON).open(encoding="utf-8") as fh:
    data = json.load(fh)
logger.info(f"  Total samples: {len(data):,}\n")

for filename, required, windows in QUERIES:
    rows = []
    for sid, sample in data.items():
        in_window = variants_in_windows(sample, windows)
        if not required.issubset(in_window) or (in_window - required):
            continue

        rows.append(
            {
                "Sample_ID": sid,
                "Batch": sample.get("batch", ""),
                "Analyzed_Range": format_analyzed_range(sample.get("intervals", {})),
                "Variants": format_variants(sample),
            },
        )

    rows.sort(key=lambda r: (r["Batch"], r["Sample_ID"]))

    with Path(filename).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    windows_str = " + ".join(f"{s}-{e}" for s, e in windows)
    pct = len(rows) / len(data) * 100
    logger.info(f"  {filename:20s}  {len(rows):>5,} samples  ({pct:.4f}%)  window: {windows_str}")

logger.info("\nDone.")
