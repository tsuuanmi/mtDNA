#!/usr/bin/env python3
"""
ms_hv23_merge.py

HV2-3 processing from:
    *_MS_trace_processed.xlsx
which is the output of ms_trace_process.py.

This version provides TWO HV2-3 views:

1) Base HV2-3 merge
   - Range by HV2-3 Regions:
       derived from MS Range of HV2F + HV3R, then trimmed to biological HV2-3
       regions only: 73-340 and 438-576
   - Variants by HV2-3 Regions:
       derived from Variant Converted 2 of HV2F + HV3R
   - duplicate variants are kept once
   - variants outside Range by HV2-3 Regions are removed
   - Mismatch by HV2-3 Regions:
       flags positions that disagree between 2F and 3R within overlap
   - Concordance_auto_2 (by HV2-3 Regions):
       optional: compare merged MS vs sample profile HV2-3 Variants using Y/Y1/Y2/F logic
       if truth columns are present; otherwise leave blank

2) Adjusted HV2-3 view
   - HV2F:
       keep only Variant Converted 2 in 73-303
       adjusted range is restricted to 73-303
   - HV3R:
       keep only Variant Converted 2 in 303-340 and 438-576
       adjusted range is restricted to those regions
   - merge adjusted traces by LID
   - Range by HV2-3 Regions Adjusted:
       use range of 2F and 3R, separated by spaces if not contiguous
   - Variants by HV2-3 Regions Adjusted:
       merged from adjusted trace variants
   - Mismatch by HV2-3 Regions Adjusted:
       mismatch positions within overlap between adjusted 2F and 3R
   - Concordance_auto_2 (by HV2-3 Regions Adjusted):
       optional: compare adjusted merged MS vs sample profile HV2-3 Variants if present

Fixes vs v18 patch
-------------------
- Use strict allele identity for trace-to-trace mismatch QC; IUPAC/PHP compatibility no longer bypasses trace mismatch.
- Add selected adjusted-trace same-position multi-call flagging, e.g. HV2-3 same-position multi-call: HV2F 152A/152C.

Fixes vs v17
------------
- Accept MS-only trace-processed input from ms_trace_process.py.
- Truthset columns are optional. If HV2-3 Variants is absent, truth status, concordance, and discordance outputs are left blank instead of raising KeyError.
- Core HV2-3 merge, adjusted range, mismatch, flag, review, and severity logic remains MS-driven.

Fixes vs v14
------------
- Add HV3R source-trace PHP QC for positions 341-437:
  >=1 PHP sends LID to Review/MEDIUM; >=3 PHP escalates to HIGH.

Fixes vs v9
-----------
- Select one HV2F and one HV3R deterministically per LID
- Add pair/truth status columns
- Use first non-empty truth instead of blind g.iloc[0]
- Do not return "Y" when there is no comparable data
- Do not merge merely adjacent ranges by default
- Preserve biological target boundaries when building merged ranges
"""

import argparse
import re
from pathlib import Path

import pandas as pd

IUPAC = {
    "W": set("AT"),
    "R": set("AG"),
    "M": set("AC"),
    "K": set("GT"),
    "Y": set("CT"),
    "S": set("CG"),
    "D": set("AGT"),
    "H": set("ACT"),
    "V": set("ACG"),
    "B": set("CGT"),
}

HV23_TARGET = [(73.0, 340.0), (438.0, 576.0)]
HV2F_ADJUSTED = [(73.0, 303.0)]
HV3R_ADJUSTED = [(303.0, 340.0), (438.0, 576.0)]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input *_MS_trace_processed.xlsx from ms_trace_process.py")
    p.add_argument("--output", required=True, help="Output xlsx with HV2-3 merged base + adjusted columns")
    return p.parse_args()


def safe_str(x):
    return "" if pd.isna(x) else str(x).strip()


def normalize_text(x):
    return " ".join(safe_str(x).split())


def is_ignored_hv23_length_variant(pos, allele):
    """
    Ignore length-variant insertions during HV2-3 merged concordance comparison:
    - 309.xC
    - 573.xC
    """
    a = safe_str(allele).replace("_HET", "_het")
    if 309 < float(pos) < 310 and a == "C":
        return True
    if 573 < float(pos) < 574 and a == "C":
        return True
    return False


def parse_range_pair(text):
    m = re.match(r"^(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)$", safe_str(text))
    if not m:
        return None
    a = float(m.group(1))
    b = float(m.group(2))
    return (min(a, b), max(a, b))


def parse_multi_ranges(text):
    out = []
    for a, b in re.findall(r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)", safe_str(text)):
        a = float(a)
        b = float(b)
        out.append((min(a, b), max(a, b)))
    return out


def intersect_ranges(r, allowed):
    if r is None:
        return []
    out = []
    for lo, hi in allowed:
        a = max(r[0], lo)
        b = min(r[1], hi)
        if a <= b:
            out.append((a, b))
    return out


def overlap_pair(r1, r2):
    if r1 is None or r2 is None:
        return None
    lo = max(r1[0], r2[0])
    hi = min(r1[1], r2[1])
    return (lo, hi) if lo <= hi else None


def normalize_ranges(ranges, merge_adjacent=False):
    if not ranges:
        return []
    ranges = sorted(ranges)
    merged = [ranges[0]]
    for a, b in ranges[1:]:
        pa, pb = merged[-1]
        cond = a <= pb if not merge_adjacent else a <= pb + 1
        if cond:
            merged[-1] = (pa, max(pb, b))
        else:
            merged.append((a, b))
    return merged


def format_num(x):
    return str(int(float(x))) if float(x).is_integer() else str(x)


def format_ranges(ranges):
    return " ".join(f"{format_num(a)}-{format_num(b)}" for a, b in ranges)


def parse_variant_tokens(text):
    out = []
    if pd.isna(text) or str(text).strip() == "":
        return out
    for tok in re.split(r"[,\s]+", str(text).strip()):
        if not tok:
            continue
        m = re.match(r"^(\d+(?:\.\d+)?)([A-Za-z_]+)$", tok)
        if not m:
            continue
        pos = float(m.group(1))
        allele = m.group(2).replace("_HET", "_het")
        out.append((pos, allele, f"{format_num(pos)}{allele}"))
    return out


def in_any_range(pos, ranges):
    return any(lo <= pos <= hi for lo, hi in ranges)


def token_sort_key(tok):
    m = re.match(r"^(\d+(?:\.\d+)?)([A-Za-z_]+)$", safe_str(tok))
    if not m:
        return (10**9, tok)
    return (float(m.group(1)), m.group(2))


def merge_unique_tokens(texts, allowed_ranges):
    seen = {}
    for txt in texts:
        for pos, allele, tok in parse_variant_tokens(txt):
            if in_any_range(pos, allowed_ranges):
                key = (pos, allele.upper())
                if key not in seen:
                    seen[key] = tok
    return " ".join(sorted(seen.values(), key=token_sort_key))


def alleles_compatible(a1, a2):
    a1 = safe_str(a1).replace("_HET", "_het")
    a2 = safe_str(a2).replace("_HET", "_het")
    c1 = a1.replace("_het", "")
    c2 = a2.replace("_het", "")
    if c1 == c2:
        return True
    s1 = IUPAC.get(c1)
    s2 = IUPAC.get(c2)
    if s1 is not None and s2 is None and c2 in s1:
        return True
    if s2 is not None and s1 is None and c1 in s2:
        return True
    if s1 is not None and s2 is not None and len(s1.intersection(s2)) > 0:
        return True
    return False


def alleles_identical_for_trace_mismatch(a1, a2):
    """
    Strict trace-to-trace QC comparison.

    Unlike concordance, trace mismatch detection should not use IUPAC/PHP
    compatibility. Two selected traces pass at an overlap position only when
    the normalized allele string is exactly identical.
    """
    a1 = safe_str(a1).replace("_HET", "_het").strip()
    a2 = safe_str(a2).replace("_HET", "_het").strip()
    return a1 == a2


def compatibility_label(ms_allele, truth_allele):
    ms = safe_str(ms_allele).replace("_HET", "_het")
    truth = safe_str(truth_allele).replace("_HET", "_het")
    ms_has_het = "_het" in ms.lower()
    truth_has_het = "_het" in truth.lower()
    ms_core = ms.replace("_het", "")
    truth_core = truth.replace("_het", "")
    if ms_core == truth_core:
        return "Y2" if (ms_has_het or truth_has_het) else "Y"
    ms_set = IUPAC.get(ms_core)
    truth_set = IUPAC.get(truth_core)
    if ms_set is not None and truth_set is None and truth_core in ms_set:
        return "Y2" if (ms_has_het or truth_has_het) else "Y1"
    if truth_set is not None and ms_set is None and ms_core in truth_set:
        return "Y2" if (ms_has_het or truth_has_het) else "Y1"
    if ms_set is not None and truth_set is not None and len(ms_set.intersection(truth_set)) > 0:
        return "Y2" if (ms_has_het or truth_has_het) else "Y1"
    return "F"


def first_nonempty(series):
    vals = [normalize_text(x) for x in series]
    vals = [x for x in vals if x != ""]
    return vals[0] if vals else ""


def consistent_nonempty(series):
    vals = [normalize_text(x) for x in series]
    vals = [x for x in vals if x != ""]
    return len(set(vals)) <= 1


def has_truth_columns(df):
    return "HV2-3 Variants" in df.columns


def get_truth_hv23_for_group(g):
    if "HV2-3 Variants" not in g.columns:
        return "", ""
    truth_status = "OK" if consistent_nonempty(g["HV2-3 Variants"]) else "Inconsistent HV2-3 truth"
    truth_hv23 = first_nonempty(g["HV2-3 Variants"])
    return truth_status, truth_hv23


def select_trace_pair(group_df):
    out = {}
    for primer in ["HV2F", "HV3R"]:
        sub = group_df[group_df["Primer"].astype(str).str.upper() == primer].copy()
        if sub.empty:
            out[primer] = None
            continue

        def span(row):
            r = parse_range_pair(row.get("MS Range", ""))
            return -1 if r is None else (r[1] - r[0])

        def trace_num(row):
            v = safe_str(row.get("Trace #", ""))
            try:
                return int(float(v))
            except Exception:
                return 10**9

        sub["_span"] = sub.apply(span, axis=1)
        sub["_trace_num"] = sub.apply(trace_num, axis=1)
        sub = sub.sort_values(["_span", "_trace_num", "Sample Name"], ascending=[False, True, True])
        out[primer] = sub.iloc[0].to_dict()
    return out


def build_region_ranges_from_rows_base(rows):
    rows = [r for r in rows if r is not None]
    rows = sorted(rows, key=lambda r: 0 if safe_str(r.get("Primer")).upper() == "HV2F" else 1)
    trimmed_by_row = []
    for r in rows:
        rr = parse_range_pair(r.get("MS Range", ""))
        trimmed_by_row.append(intersect_ranges(rr, HV23_TARGET))
    flat = [x for sub in trimmed_by_row for x in sub]
    if len(rows) >= 2:
        left = trimmed_by_row[0]
        right = trimmed_by_row[1]
        if left and right:
            first_seg = left[0]
            last_seg = right[-1]
            ov = overlap_pair(first_seg, last_seg)
            if ov is not None:
                span = (first_seg[0], last_seg[1])
                return intersect_ranges(span, HV23_TARGET)
    return normalize_ranges(flat, merge_adjacent=False)


def build_region_ranges_from_rows_adjusted(rows):
    rows = [r for r in rows if r is not None]
    rows = sorted(rows, key=lambda r: 0 if safe_str(r.get("Primer")).upper() == "HV2F" else 1)

    left = parse_multi_ranges(rows[0].get("HV2-3 Trace Range Adjusted", "")) if len(rows) >= 1 else []
    right = parse_multi_ranges(rows[1].get("HV2-3 Trace Range Adjusted", "")) if len(rows) >= 2 else []

    flat = left + right

    if left and right:
        first_seg = left[0]
        last_seg = right[-1]
        ov = overlap_pair(first_seg, last_seg)
        if ov is not None:
            span = (first_seg[0], last_seg[1])
            return intersect_ranges(span, HV23_TARGET)

    return normalize_ranges(flat, merge_adjacent=False)


def safe_numeric_sort(values):
    def _key(x):
        try:
            return float(x)
        except Exception:
            return 10**9

    return sorted(values, key=_key)


def variant_allele_sets_in_ranges(text, ranges):
    out = {}
    for pos, allele, _ in parse_variant_tokens(text):
        if not in_any_range(pos, ranges):
            continue
        norm_allele = safe_str(allele).replace("_HET", "_het").strip()
        out.setdefault(pos, set()).add(norm_allele)
    return out


def allele_sets_identical_for_trace_mismatch(alleles1, alleles2):
    if alleles1 is None or alleles2 is None:
        return False
    if len(alleles1) != len(alleles2):
        return False
    for a in alleles1:
        if not any(alleles_identical_for_trace_mismatch(a, b) for b in alleles2):
            return False
    return True


def mismatch_positions_between_two_rows(row1, row2, variant_key, range_ranges1, range_ranges2):
    if row1 is None or row2 is None:
        return []
    overlap_ranges = []
    for a in range_ranges1:
        for b in range_ranges2:
            ov = overlap_pair(a, b)
            if ov is not None:
                overlap_ranges.append(ov)
    overlap_ranges = normalize_ranges(overlap_ranges, merge_adjacent=False)

    if not overlap_ranges:
        return []

    d1 = variant_allele_sets_in_ranges(row1.get(variant_key, ""), overlap_ranges)
    d2 = variant_allele_sets_in_ranges(row2.get(variant_key, ""), overlap_ranges)

    positions = sorted(set(d1.keys()).union(d2.keys()))
    mismatches = []
    for pos in positions:
        a1 = d1.get(pos)
        a2 = d2.get(pos)
        if a1 is None or a2 is None or not allele_sets_identical_for_trace_mismatch(a1, a2):
            mismatches.append(format_num(pos))
    return mismatches


def concordance_hv23(merged_ranges, merged_variants_text, truth_text):
    """
    Compare merged HV2-3 variants vs truth within merged_ranges only.

    Ignore the following length-variant insertions during comparison:
    - 309.xC
    - 573.xC
    """
    if not merged_ranges:
        return ""

    ms_all = {
        pos: allele for pos, allele, _ in parse_variant_tokens(merged_variants_text) if in_any_range(pos, merged_ranges)
    }
    truth_all = {pos: allele for pos, allele, _ in parse_variant_tokens(truth_text) if in_any_range(pos, merged_ranges)}

    ms = {pos: allele for pos, allele in ms_all.items() if not is_ignored_hv23_length_variant(pos, allele)}
    truth = {pos: allele for pos, allele in truth_all.items() if not is_ignored_hv23_length_variant(pos, allele)}

    if not ms and not truth:
        return ""

    status = "Y"
    for pos, truth_allele in truth.items():
        if pos not in ms:
            return "F"
        verdict = compatibility_label(ms[pos], truth_allele)
        if verdict == "F":
            return "F"
        if verdict == "Y2":
            status = "Y2"
        elif verdict == "Y1" and status == "Y":
            status = "Y1"

    for pos in ms:
        if pos not in truth:
            return "F"

    return status


# =========================
# HV2-3 REVIEW FLAGS / DISCORDANCE
# =========================
PHP_CODES = set("RMSKYWVHBDN")


def covers_segment(ranges, target):
    lo, hi = target
    return any(a <= lo and b >= hi for a, b in ranges)


def is_common_hv23_insertion(pos, allele):
    a = safe_str(allele).replace("_HET", "_het")
    if 309 < float(pos) < 310 and a == "C":
        return True
    if 315 < float(pos) < 316 and a == "C":
        return True
    if 573 < float(pos) < 574 and a == "C":
        return True
    if 452 < float(pos) < 455 and a == "T":
        return True
    return False


def is_common_hv23_deletion(pos, allele):
    a = safe_str(allele).upper()
    if "DEL" not in a:
        return False
    if float(pos) == 249:
        return True
    if 513 <= float(pos) <= 525:
        return True
    return False


def hv3r_php_341_437_from_selected_trace(hv3r_row):
    """
    Source-trace QC for selected HV3R.

    Rationale:
    Some poor-quality HV3R traces show PHP/IUPAC calls in 341-437.
    These sites may be outside the adjusted merged profile and therefore would
    otherwise disappear from the merged-variant based Review logic.
    """
    if hv3r_row is None:
        return 0, ""

    sites = []
    for pos, allele, tok in parse_variant_tokens(hv3r_row.get("Variant Converted 2", "")):
        if 341 <= pos <= 437 and allele.upper() in PHP_CODES:
            sites.append(tok)

    # Preserve numeric order while avoiding duplicate calls at the same site.
    seen = {}
    for tok in sites:
        seen[tok] = tok
    ordered = sorted(seen.values(), key=token_sort_key)
    return len(ordered), " ".join(ordered)


def find_same_position_multicalls(trace_row, variant_key="HV2-3 Trace Variants Adjusted"):
    """
    Detect internal multi-call positions in one selected adjusted trace.

    A multi-call exists when the same selected trace has more than one
    distinct normalized allele at the same position. Repeated identical calls
    are treated as duplicates, not as multi-calls.
    """
    if trace_row is None:
        return []

    by_pos = {}
    for pos, allele, tok in parse_variant_tokens(trace_row.get(variant_key, "")):
        norm_allele = safe_str(allele).replace("_HET", "_het").strip()
        by_pos.setdefault(pos, {})[norm_allele.upper()] = f"{format_num(pos)}{norm_allele}"

    primer = safe_str(trace_row.get("Primer", "")).upper()
    out = []
    for pos in sorted(by_pos):
        tokens = sorted(by_pos[pos].values(), key=token_sort_key)
        if len(tokens) > 1:
            out.append(f"{primer} " + "/".join(tokens))
    return out


def format_same_position_multicall_text(items):
    seen = []
    for item in items:
        item = safe_str(item)
        if item and item not in seen:
            seen.append(item)
    return ", ".join(seen)


def build_hv23_flags(
    pair_status,
    truth_status,
    adjusted_ranges,
    adjusted_variants,
    adjusted_mismatch,
    hv3r_php_341_437_count=0,
    same_position_multicall_text="",
):
    flags = []
    toks = parse_variant_tokens(adjusted_variants)

    required_segments = [(73.0, 303.0), (303.0, 340.0), (438.0, 576.0)]
    if not all(covers_segment(adjusted_ranges, seg) for seg in required_segments):
        flags.append("No full region")

    if safe_str(pair_status) != "OK":
        flags.append("Pair status not OK")

    if safe_str(truth_status) not in {"", "OK"}:
        flags.append("Truth status not OK")

    if safe_str(adjusted_mismatch) != "":
        flags.append("Trace mismatch")

    if safe_str(same_position_multicall_text):
        flags.append(f"HV2-3 same-position multi-call: {safe_str(same_position_multicall_text)}")

    if hv3r_php_341_437_count >= 2:
        flags.append("HV3R PHP 341-437")

    if not any(abs(pos - 315.1) < 1e-9 and allele.upper() == "C" for pos, allele, _ in toks):
        flags.append("No 315.1C")

    if any(allele.upper() in PHP_CODES for _, allele, _ in toks):
        flags.append("PHP")

    if "het" in safe_str(adjusted_variants).lower():
        flags.append("het")

    if any(int(float(pos)) == 310 for pos, _, _ in toks):
        flags.append("310 variant")

    if any(((303 <= pos <= 309) or (311 <= pos <= 315)) and allele.upper() != "C" for pos, allele, _ in toks):
        flags.append("non-C in HV2 polyC")

    if any((568 <= pos <= 573) and allele.upper() != "C" for pos, allele, _ in toks):
        flags.append("non-C in HV3 polyC")

    for pos, allele, _ in toks:
        if not float(pos).is_integer():
            if not is_common_hv23_insertion(pos, allele):
                flags.append("unexpected insertion")
                break

    for pos, allele, _ in toks:
        if "DEL" in allele.upper():
            if not is_common_hv23_deletion(pos, allele):
                flags.append("unexpected deletion")
                break

    return "; ".join(sorted(set(flags)))


def hv23_review_class(flag_list):
    return "Auto - no flag" if safe_str(flag_list) == "" else "Review"


def variant_map_in_ranges(text, ranges, ignore_length_variants=False):
    out = {}
    for pos, allele, tok in parse_variant_tokens(text):
        if not in_any_range(pos, ranges):
            continue
        if ignore_length_variants and is_ignored_hv23_length_variant(pos, allele):
            continue
        out[pos] = (allele, tok)
    return out


def discordance_reason_hv23(merged_ranges, ms_text, truth_text):
    if not merged_ranges:
        return ""

    ms = variant_map_in_ranges(ms_text, merged_ranges, ignore_length_variants=True)
    truth = variant_map_in_ranges(truth_text, merged_ranges, ignore_length_variants=True)

    if not ms and not truth:
        return ""

    missing = []
    extra = []
    mismatch = []

    for pos in sorted(set(truth) - set(ms)):
        missing.append(truth[pos][1])

    for pos in sorted(set(ms) - set(truth)):
        extra.append(ms[pos][1])

    for pos in sorted(set(ms) & set(truth)):
        ms_allele, ms_tok = ms[pos]
        truth_allele, truth_tok = truth[pos]
        verdict = compatibility_label(ms_allele, truth_allele)
        if verdict == "F":
            mismatch.append(f"{ms_tok}!={truth_tok}")

    parts = []
    if missing:
        parts.append("Missing [" + ", ".join(missing) + "]")
    if extra:
        parts.append("Extra [" + ", ".join(extra) + "]")
    if mismatch:
        parts.append("Mismatch [" + ", ".join(mismatch) + "]")

    return "; ".join(parts)


def hv23_severity(flag_list, concordance, discordance_reason, hv3r_php_341_437_count=0):
    flags = set([x.strip() for x in safe_str(flag_list).split(";") if x.strip()])
    conc = safe_str(concordance)

    if safe_str(flag_list) == "" and conc in {"Y", "Y1", "Y2"}:
        return "AUTO"

    if hv3r_php_341_437_count >= 3:
        return "HIGH"

    high_flags = {
        "No full region",
        "Pair status not OK",
        "Trace mismatch",
        "unexpected deletion",
        "unexpected insertion",
    }
    medium_flags = {
        "PHP",
        "het",
        "Truth status not OK",
        "310 variant",
        "non-C in HV2 polyC",
        "non-C in HV3 polyC",
        "HV3R PHP 341-437",
    }
    low_flags = {"No 315.1C"}

    if conc == "F" and safe_str(discordance_reason):
        return "HIGH"
    if flags.intersection(high_flags):
        return "HIGH"
    if hv3r_php_341_437_count >= 2:
        return "MEDIUM"
    if flags.intersection(medium_flags):
        return "MEDIUM"
    if any(flag.startswith("HV2-3 same-position multi-call:") for flag in flags):
        return "MEDIUM"
    if flags.intersection(low_flags):
        return "LOW"
    if conc == "F":
        return "HIGH"
    return "LOW" if flags else "AUTO"


def trim_trace_adjusted(primer, ms_range_text, variant_converted2_text):
    primer = safe_str(primer).upper()
    if primer == "HV2F":
        allowed = HV2F_ADJUSTED
    elif primer == "HV3R":
        allowed = HV3R_ADJUSTED
    else:
        return "", ""
    rr = parse_range_pair(ms_range_text)
    trimmed_ranges = intersect_ranges(rr, allowed)
    kept_map = {}
    for pos, allele, tok in parse_variant_tokens(variant_converted2_text):
        if in_any_range(pos, trimmed_ranges):
            key = (pos, allele.upper())
            if key not in kept_map:
                kept_map[key] = tok
    kept = sorted(kept_map.values(), key=token_sort_key)
    return format_ranges(trimmed_ranges), " ".join(kept)


def main():
    args = parse_args()
    df = pd.read_excel(args.input)
    df.columns = [str(c).strip() for c in df.columns]

    df = df[df["Primer"].astype(str).str.upper().isin(["HV2F", "HV3R"])].copy()

    df["HV2-3 Pair Status"] = ""
    df["HV2-3 Truth Status"] = ""
    df["Range by HV2-3 Regions"] = ""
    df["Variants by HV2-3 Regions"] = ""
    df["Mismatch by HV2-3 Regions"] = ""
    df["Concordance_auto_2 (by HV2-3 Regions)"] = ""

    df["HV2-3 Trace Range Adjusted"] = ""
    df["HV2-3 Trace Variants Adjusted"] = ""
    df["Range by HV2-3 Regions Adjusted"] = ""
    df["Variants by HV2-3 Regions Adjusted"] = ""
    df["Mismatch by HV2-3 Regions Adjusted"] = ""
    df["Concordance_auto_2 (by HV2-3 Regions Adjusted)"] = ""

    for idx, row in df.iterrows():
        adj_range, adj_vars = trim_trace_adjusted(
            row.get("Primer", ""), row.get("MS Range", ""), row.get("Variant Converted 2", "")
        )
        df.at[idx, "HV2-3 Trace Range Adjusted"] = adj_range
        df.at[idx, "HV2-3 Trace Variants Adjusted"] = adj_vars

    total_lid = 0
    ok_pairs = 0
    missing_hv2f = 0
    missing_hv3r = 0
    missing_both = 0
    truth_inconsistent = 0

    for lid, g in df.groupby("LID", dropna=False):
        total_lid += 1
        pair = select_trace_pair(g)
        hv2f = pair["HV2F"]
        hv3r = pair["HV3R"]
        pair_rows = [hv2f, hv3r]

        if hv2f is not None and hv3r is not None:
            pair_status = "OK"
            ok_pairs += 1
        elif hv2f is not None:
            pair_status = "Missing HV3R"
            missing_hv3r += 1
        elif hv3r is not None:
            pair_status = "Missing HV2F"
            missing_hv2f += 1
        else:
            pair_status = "Missing HV2F and HV3R"
            missing_both += 1

        truth_status, truth_hv23 = get_truth_hv23_for_group(g)
        if truth_status not in {"", "OK"}:
            truth_inconsistent += 1

        # Base merge
        base_ranges = build_region_ranges_from_rows_base(pair_rows)
        base_ranges_text = format_ranges(base_ranges)
        base_variants = merge_unique_tokens(
            [r.get("Variant Converted 2", "") for r in pair_rows if r is not None], base_ranges
        )

        base_range1 = (
            intersect_ranges(parse_range_pair(hv2f.get("MS Range", "")), HV23_TARGET) if hv2f is not None else []
        )
        base_range2 = (
            intersect_ranges(parse_range_pair(hv3r.get("MS Range", "")), HV23_TARGET) if hv3r is not None else []
        )
        base_mismatch = ""
        if hv2f is not None and hv3r is not None:
            mism = mismatch_positions_between_two_rows(
                hv2f, hv3r, variant_key="Variant Converted 2", range_ranges1=base_range1, range_ranges2=base_range2
            )
            base_mismatch = " ".join(safe_numeric_sort(set(mism))) if mism else ""

        base_conc = concordance_hv23(base_ranges, base_variants, truth_hv23) if truth_hv23 else ""

        # Adjusted merge
        adj_ranges = build_region_ranges_from_rows_adjusted(pair_rows)
        adj_ranges_text = format_ranges(adj_ranges)
        adj_variants = merge_unique_tokens(
            [r.get("HV2-3 Trace Variants Adjusted", "") for r in pair_rows if r is not None], adj_ranges
        )

        adj_range1 = parse_multi_ranges(hv2f.get("HV2-3 Trace Range Adjusted", "")) if hv2f is not None else []
        adj_range2 = parse_multi_ranges(hv3r.get("HV2-3 Trace Range Adjusted", "")) if hv3r is not None else []
        adj_mismatch = ""
        if hv2f is not None and hv3r is not None:
            mism = mismatch_positions_between_two_rows(
                hv2f,
                hv3r,
                variant_key="HV2-3 Trace Variants Adjusted",
                range_ranges1=adj_range1,
                range_ranges2=adj_range2,
            )
            adj_mismatch = " ".join(safe_numeric_sort(set(mism))) if mism else ""

        adj_conc = concordance_hv23(adj_ranges, adj_variants, truth_hv23) if truth_hv23 else ""

        for i in g.index:
            df.at[i, "HV2-3 Pair Status"] = pair_status
            df.at[i, "HV2-3 Truth Status"] = truth_status
            df.at[i, "Range by HV2-3 Regions"] = base_ranges_text
            df.at[i, "Variants by HV2-3 Regions"] = base_variants
            df.at[i, "Mismatch by HV2-3 Regions"] = base_mismatch
            df.at[i, "Concordance_auto_2 (by HV2-3 Regions)"] = base_conc

            df.at[i, "Range by HV2-3 Regions Adjusted"] = adj_ranges_text
            df.at[i, "Variants by HV2-3 Regions Adjusted"] = adj_variants
            df.at[i, "Mismatch by HV2-3 Regions Adjusted"] = adj_mismatch
            df.at[i, "Concordance_auto_2 (by HV2-3 Regions Adjusted)"] = adj_conc

    final_cols = [
        "Batch",
        "Well",
        "LID",
        "Primer",
        "Trace #",
        "Sample Name",
        "# of Mutation",
        "Read Start",
        "Read End",
        "MS Range",
        "Variant Summary",
        "Variant Converted",
        "Variant Converted 2",
        "HV2-3 Pair Status",
        "HV2-3 Truth Status",
        "HV2-3 Trace Range Adjusted",
        "HV2-3 Trace Variants Adjusted",
        "Range by HV2-3 Regions",
        "Variants by HV2-3 Regions",
        "Mismatch by HV2-3 Regions",
        "Concordance_auto_2 (by HV2-3 Regions)",
        "Range by HV2-3 Regions Adjusted",
        "Variants by HV2-3 Regions Adjusted",
        "Mismatch by HV2-3 Regions Adjusted",
        "Concordance_auto_2 (by HV2-3 Regions Adjusted)",
        "Sample Profiles - Range",
        "Sample Profile - Flag",
        "HV1 Range",
        "HV1 Variants",
        "HV2-3 Range",
        "HV2-3 Variants",
        "Concordance_auto",
        "Concordance_auto_2",
        "HV1_PolyC_flag",
        "HV2-3_PolyC_flag",
        "Error_type",
    ]
    final_cols = [c for c in final_cols if c in df.columns]
    traces_df = df[final_cols].copy()

    lid_records = []
    for lid, g in traces_df.groupby("LID", dropna=False):
        row = g.iloc[0].to_dict()

        adj_ranges = parse_multi_ranges(row.get("Range by HV2-3 Regions Adjusted", ""))
        adj_vars = row.get("Variants by HV2-3 Regions Adjusted", "")
        adj_mismatch = row.get("Mismatch by HV2-3 Regions Adjusted", "")
        truth_hv23 = row.get("HV2-3 Variants", "")
        adj_conc = row.get("Concordance_auto_2 (by HV2-3 Regions Adjusted)", "")

        selected_pair = select_trace_pair(g)
        hv3r_php_count, hv3r_php_sites = hv3r_php_341_437_from_selected_trace(selected_pair.get("HV3R"))

        same_position_multicall_items = []
        same_position_multicall_items.extend(find_same_position_multicalls(selected_pair.get("HV2F")))
        same_position_multicall_items.extend(find_same_position_multicalls(selected_pair.get("HV3R")))
        same_position_multicall_text = format_same_position_multicall_text(same_position_multicall_items)

        flags = build_hv23_flags(
            row.get("HV2-3 Pair Status", ""),
            row.get("HV2-3 Truth Status", ""),
            adj_ranges,
            adj_vars,
            adj_mismatch,
            hv3r_php_count,
            same_position_multicall_text,
        )

        discordance = ""
        if safe_str(adj_conc) == "F":
            discordance = discordance_reason_hv23(adj_ranges, adj_vars, truth_hv23)

        severity = hv23_severity(flags, adj_conc, discordance, hv3r_php_count)
        review = hv23_review_class(flags)

        lid_records.append(
            {
                "Batch": row.get("Batch", ""),
                "Well": row.get("Well", ""),
                "LID": lid,
                "HV2-3 Pair Status": row.get("HV2-3 Pair Status", ""),
                "HV2-3 Truth Status": row.get("HV2-3 Truth Status", ""),
                "Range by HV2-3 Regions": row.get("Range by HV2-3 Regions", ""),
                "Variants by HV2-3 Regions": row.get("Variants by HV2-3 Regions", ""),
                "Mismatch by HV2-3 Regions": row.get("Mismatch by HV2-3 Regions", ""),
                "Concordance_auto_2 (by HV2-3 Regions)": row.get("Concordance_auto_2 (by HV2-3 Regions)", ""),
                "Range by HV2-3 Regions Adjusted": row.get("Range by HV2-3 Regions Adjusted", ""),
                "Variants by HV2-3 Regions Adjusted": adj_vars,
                "Mismatch by HV2-3 Regions Adjusted": adj_mismatch,
                "Concordance_auto_2 (by HV2-3 Regions Adjusted)": adj_conc,
                "HV3R PHP 341-437 Count": hv3r_php_count,
                "HV3R PHP 341-437 Sites": hv3r_php_sites,
                "HV2-3 Flag list": flags,
                "HV2-3 Review Class": review,
                "HV2-3 Severity": severity,
                "HV2-3 Discordance reason": discordance,
                "HV2-3 Range": row.get("HV2-3 Range", ""),
                "HV2-3 Variants": truth_hv23,
                "Sample Profiles - Range": row.get("Sample Profiles - Range", ""),
                "Sample Profile - Flag": row.get("Sample Profile - Flag", ""),
            }
        )

    lid_df = pd.DataFrame(lid_records)

    summary_rows = []
    summary_rows.append({"Metric": "Total LID", "Value": len(lid_df)})

    for k, v in lid_df["HV2-3 Review Class"].value_counts(dropna=False).items():
        summary_rows.append({"Metric": f"Review Class - {k}", "Value": v})

    for k, v in lid_df["HV2-3 Severity"].value_counts(dropna=False).items():
        summary_rows.append({"Metric": f"Severity - {k}", "Value": v})

    for k, v in lid_df["HV2-3 Pair Status"].value_counts(dropna=False).items():
        summary_rows.append({"Metric": f"Pair Status - {k}", "Value": v})

    for k, v in lid_df["Concordance_auto_2 (by HV2-3 Regions Adjusted)"].value_counts(dropna=False).items():
        label = k if safe_str(k) else "(blank)"
        summary_rows.append({"Metric": f"Adjusted Concordance - {label}", "Value": v})

    flag_counter = {}
    for flist in lid_df["HV2-3 Flag list"]:
        if safe_str(flist):
            for f in safe_str(flist).split("; "):
                flag_counter[f] = flag_counter.get(f, 0) + 1
    for k, v in sorted(flag_counter.items(), key=lambda x: (-x[1], x[0])):
        summary_rows.append({"Metric": f"Flag - {k}", "Value": v})

    summary_df = pd.DataFrame(summary_rows)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        traces_df.to_excel(writer, sheet_name="Traces", index=False)
        lid_df.to_excel(writer, sheet_name="LID", index=False)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)

        for sheet_name, worksheet in writer.sheets.items():
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for col in worksheet.columns:
                max_len = 0
                col_letter = col[0].column_letter
                for cell in col:
                    value = "" if cell.value is None else str(cell.value)
                    max_len = max(max_len, len(value))
                worksheet.column_dimensions[col_letter].width = min(max_len + 2, 50)

    print(f"Saved: {args.output}")
    print("Summary:")
    print(f"  Total LID: {len(lid_df)}")
    for k, v in lid_df["HV2-3 Review Class"].value_counts(dropna=False).items():
        print(f"  Review Class - {k}: {v}")
    for k, v in lid_df["HV2-3 Severity"].value_counts(dropna=False).items():
        print(f"  Severity - {k}: {v}")
    for k, v in lid_df["Concordance_auto_2 (by HV2-3 Regions Adjusted)"].value_counts(dropna=False).items():
        label = k if safe_str(k) else "(blank)"
        print(f"  Adjusted Concordance - {label}: {v}")


if __name__ == "__main__":
    main()
