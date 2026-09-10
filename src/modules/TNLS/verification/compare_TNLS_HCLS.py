#!/usr/bin/env python
# noqa: N999
"""
compare_TNLS_HCLS.py
-------------------
Compares HCLS base samples against TNLS target samples using a three-stage
matching pipeline, then computes summary statistics (1:1 / 1:N matching ratios,
exclusivity, variant totals) and writes annotated TSV outputs.

Input / output paths
--------------------
  TSV output  : results/output_compare_tnls_hcls.tsv
  Family summary: results/stats_TNLS_HCLS_family_summary.tsv
  Stats with meta: results/stats_TNLS_HCLS_<1toN>_with_metadata.tsv

CLI summary
-----------
  compare_TNLS_HCLS.py --base-samples HCLS.json --target-samples TNLS_merged.json \\
      --matching-rule rule.json [--filter-list lid.txt] [--tnls-metadata meta.xlsx] \\
      [--max-ratio N] [--results-dir DIR]
"""

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from src.core.variants import (
    format_variants_simplified,
    normalize_position,
    parse_variant_position_allele,
    pos_sort_key,
)
from src.modules.TNLS.verification.compare_1n_class import compare_1n
from src.modules.TNLS.verification.family_profiles import group_families, sample_coverage_bp
from src.modules.TNLS.verification.utils import (
    check_number_in_intervals,
    format_sample_variants,
    load_json,
)

# ── Statistics helpers ───────────────────────────────────────────────────────────


def _normalize_tnls_id(value: object | None) -> str:
    if pd.isna(value):  # type: ignore[reportGeneralTypeIssues]
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


def _build_tnls_metadata_lookup(tnls_metadata_path: Path | None) -> tuple[list[str] | None, dict[str, Any]]:
    """Return (columns, lookup_dict) for TNLS metadata, or (None, {})."""
    if not tnls_metadata_path or not tnls_metadata_path.exists():
        return None, {}

    df = pd.read_excel(tnls_metadata_path, sheet_name=0, dtype=str)
    df.columns = [str(col).strip() for col in df.columns]
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

    if "mtDNA LID" not in df.columns:
        logger.warning("Column 'mtDNA LID' not found in TNLS_metadata.xlsx")
        return None, {}

    df["mtDNA LID"] = df["mtDNA LID"].apply(_normalize_tnls_id)
    df = df[df["mtDNA LID"] != ""]
    df = df.drop_duplicates(subset=["mtDNA LID"], keep="first")

    columns = [col for col in df.columns if col != "mtDNA LID"]

    def _clean(v: object | None) -> str:
        if pd.isna(v):  # type: ignore[reportGeneralTypeIssues]
            return ""
        text = str(v).replace("\n", " ").replace("\r", " ").strip()
        return "" if text.lower() == "nan" else text

    lookup = df.set_index("mtDNA LID")[columns].map(_clean).to_dict(orient="index")
    logger.info(f"TNLS metadata rows: {len(lookup)}")
    return columns, lookup


def _write_clean_output(ratio_output_df: pd.DataFrame, results_dir: Path) -> None:
    """Write stats_TNLS_HCLS_clean.tsv (hardcoded column subset, bare names)."""
    clean_columns = [
        "HCLS",
        "status",
        "TNLS",
        "overlap_bp",
        "overlap_intervals",
        "matched_variants",
        "matched_bases",
        "variant_match",
        "HCLS_variants",
        "TNLS_variants",
        "CCCD",
        "họ và tên",
        "Chương trình thu mẫu",
        "Thông tin",
        "Liệt sĩ",
    ]
    available_clean_columns = [col for col in clean_columns if col in ratio_output_df.columns]
    clean_output_df = ratio_output_df[available_clean_columns].copy()
    output_clean = results_dir / "stats_TNLS_HCLS_clean.tsv"
    clean_output_df.to_csv(output_clean, index=False, sep="\t", encoding="utf-8-sig")
    logger.info(f"Clean output saved -> {output_clean}")


def _distinct_variant_totals(df: pd.DataFrame) -> pd.DataFrame:
    """Per-HCLS distinct matched variants, sites, and base pairs.
    Unions the matched variants across every matched TNLS pair for an HCLS so
    duplicate TNLS with identical variants are not double-counted (e.g. a 1:2
    match against two identical TNLS reports 5, not 10). The union string is
    parsed from the per-pair ``variant_match`` column (simplified display form)
    and rendered back (e.g. ``16093C 16189C 309.1C``). Falls back to per-pair
    sums when ``variant_match`` is unavailable (string set to ``"None"``).
    """
    if "variant_match" not in df.columns:
        fallback = pd.DataFrame(
            df[["base_sample", "matched_variants", "matched_bases"]].groupby("base_sample", as_index=False).sum()
        ).rename(columns={"matched_variants": "matched_variant_sites", "matched_bases": "matched_variant_bp"})
        fallback["matched_variants"] = "None"
        return fallback

    def _parse_simplified(raw: object) -> list[tuple[str, str]]:
        if not isinstance(raw, str) or not raw.strip() or raw.strip() == "None":
            return []
        pairs: list[tuple[str, str]] = []
        for token in raw.split():
            pos, allele = parse_variant_position_allele(token)
            if pos:
                pairs.append((pos, "-" if allele == "DEL" else allele))
        return pairs

    parsed = df["variant_match"].apply(_parse_simplified)
    rows: list[dict[str, Any]] = []
    for hcls, group in df.groupby("base_sample"):
        union: dict[str, str] = {}
        for pairs in parsed[group.index]:
            for pos, seq in pairs:
                union.setdefault(pos, seq)
        tokens = [
            pos + ("DEL" if seq == "-" else seq)
            for pos, seq in sorted(union.items(), key=lambda kv: pos_sort_key(kv[0]))
        ]
        rows.append(
            {
                "base_sample": hcls,
                "matched_variants": " ".join(tokens) if tokens else "None",
                "matched_variant_sites": len(union),
                "matched_variant_bp": sum(len(seq) for seq in union.values()),
            }
        )
    return pd.DataFrame(rows)


def _hcls_match_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Per-HCLS match summary: matched TNLS, ratio, exclusivity, variant totals.

    One row per HCLS base sample. ``matched_tnls_list`` is a Python list of TNLS
    ids (sorted); ``matched_variant_sites``/``matched_variant_bp`` are the
    distinct union across all matched TNLS pairs (so duplicate TNLS with
    identical variants are not double-counted). Shared by the per-pair stats and
    the family summary so the two views stay consistent.
    """
    hcls_to_tnls = (
        df.groupby("base_sample")["target_sample"]
        .apply(lambda x: sorted(x.unique().tolist()))
        .reset_index()
        .rename(columns={"target_sample": "matched_tnls_list"})
    )
    hcls_to_tnls["tnls_count"] = hcls_to_tnls["matched_tnls_list"].apply(len)
    hcls_to_tnls["ratio"] = hcls_to_tnls["tnls_count"].apply(lambda n: f"1:{n}")

    agg_df = _distinct_variant_totals(df)
    hcls_to_tnls = hcls_to_tnls.merge(agg_df, on="base_sample", how="left")  # type: ignore[reportArgumentType]
    hcls_to_tnls["matched_variant_sites"] = hcls_to_tnls["matched_variant_sites"].fillna(0).astype(int)  # type: ignore[reportAttributeAccessIssue,reportCallIssue]
    hcls_to_tnls["matched_variant_bp"] = hcls_to_tnls["matched_variant_bp"].fillna(0).astype(int)  # type: ignore[reportAttributeAccessIssue,reportCallIssue]
    hcls_to_tnls["matched_variants"] = hcls_to_tnls["matched_variants"].fillna("None")  # type: ignore[reportAttributeAccessIssue,reportCallIssue]
    return hcls_to_tnls


def compute_stats(  # noqa: C901,PLR0912,PLR0915
    results_dir: Path,
    input_tsv: Path,
    base_samples: dict[str, Any],
    target_samples: dict[str, Any],
    filter_list: Path | None,
    tnls_metadata: Path | None,
    max_ratio: int = 5,
) -> None:
    """Compute and write matching statistics from the comparison TSV."""

    # Load comparison TSV
    df = pd.read_csv(input_tsv, dtype=str, sep="\t")
    df["base_sample"] = df["base_sample"].str.strip()
    df["target_sample"] = df["target_sample"].str.strip()

    numeric_cols = ["overlap_bp", "no_mismatch_pos", "mismatch", "matched_variants", "matched_bases"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)  # type: ignore[reportAttributeAccessIssue,reportCallIssue]

    n_pairs = len(df)
    n_hcls_matched = int(df["base_sample"].nunique())  # type: ignore[reportCallIssue]
    n_tnls_matched = int(df["target_sample"].nunique())  # type: ignore[reportCallIssue]
    n_hcls_total = len(base_samples)
    n_tnls_raw = len(target_samples)

    # Filter-list analysis
    all_filters: set[str] = set()
    merged_keys: set[str] = set()
    if filter_list and filter_list.exists():
        with filter_list.open(encoding="utf-8") as f:
            all_filters = {line.strip() for line in f if line.strip()}
        merged_keys = set(target_samples.keys())

    n_filter_total = len(all_filters)
    n_filter_found = len(all_filters & merged_keys)

    tnls_filtered_data: dict = {}
    n_tnls_filtered = 0
    n_no_snps = 0
    if all_filters and merged_keys:
        tnls_filtered_data = {k: v for k, v in target_samples.items() if k in all_filters}
        for v in tnls_filtered_data.values():
            if v.get("variants", {}).get("snps") or v.get("variants", {}).get("insertions"):
                n_tnls_filtered += 1
            else:
                n_no_snps += 1

    logger.info(f"Stats — rows loaded: {n_pairs}")
    logger.info(f"HCLS: {n_hcls_total} total  (matched: {n_hcls_matched})")
    logger.info(f"TNLS raw: {n_tnls_raw}  (matched in comparison: {n_tnls_matched})")
    if n_filter_total:
        logger.info(
            f"Filter list: {n_filter_total}  (found in TNLS: {n_filter_found})"
            f"  (with variants: {n_tnls_filtered}, without: {n_no_snps})",
        )

    # Per-HCLS match summary (matched TNLS, ratio, exclusivity, distinct variant
    # totals). Shared with the family summary — see `_hcls_match_summary`.
    hcls_to_tnls = _hcls_match_summary(df)

    # Build expanded row-level output (1:1 to 1:max_ratio) with optional TNLS metadata
    metadata_columns, tnls_metadata_lookup = _build_tnls_metadata_lookup(tnls_metadata)

    ratio_df = hcls_to_tnls[hcls_to_tnls["tnls_count"].between(1, max_ratio)]
    expanded_rows = []
    for _, row in ratio_df.iterrows():
        hcls = row["base_sample"]
        status = row["ratio"]
        for tnls in row["matched_tnls_list"]:
            pair_row = df[(df["base_sample"] == hcls) & (df["target_sample"] == tnls)]
            if pair_row.empty:
                continue
            pr = pair_row.iloc[0]
            row_data = {
                "HCLS": hcls,
                "status": status,
                "TNLS": tnls,
                "overlap_bp": pr.get("overlap_bp", ""),
                "overlap_intervals": pr.get("list_overlap", ""),
                "matched_variants": pr.get("matched_variants", ""),
                "matched_bases": pr.get("matched_bases", ""),
                "variant_match": pr.get("variant_match", "None"),
                "HCLS_variants": format_sample_variants(base_samples[str(hcls)]),
                "TNLS_variants": format_sample_variants(target_samples[str(tnls)]),
            }
            if metadata_columns:
                tnls_meta = tnls_metadata_lookup.get(_normalize_tnls_id(tnls), {})
                for col in metadata_columns:
                    row_data[f"TNLS_{col}"] = tnls_meta.get(col)
            expanded_rows.append(row_data)

    output_with_meta = results_dir / "stats_TNLS_HCLS_with_metadata.tsv"
    ratio_output_df = pd.DataFrame(expanded_rows)
    ratio_output_df.to_csv(output_with_meta, index=False, sep="\t", encoding="utf-8-sig")
    logger.info(f"1:1 to 1:{max_ratio} output saved -> {output_with_meta}")

    # ── Clean subset (stats_TNLS_HCLS_clean.tsv) ──────────────────────────────────
    _write_clean_output(ratio_output_df, results_dir)


# ── Family grouping by mtDNA profile ──────────────────────────────────────────────


def _format_region_variants(sample: dict[str, Any], region: str) -> str:
    """Space-joined variant string for one HV region (SNPs + insertions + deletions).

    Reference-only display; not used for matching. Reuses
    ``format_variants_simplified`` (``src/core/variants.py``) which space-joins
    the per-record simplified tokens in position order (e.g.
    ``16093C 16189C 309.1C``). Returns ``""`` when the region has no variants.
    """
    intervals = sample.get("intervals", {}).get(region, [])
    records: list[dict[str, Any]] = []
    variants = sample.get("variants", {})
    for vtype in ("snps", "insertions", "deletions"):
        for rec in variants.get(vtype, []):
            pos = rec.get("pos")
            if pos is None or not check_number_in_intervals(intervals, normalize_position(pos)):
                continue
            records.append(rec)
    if not records:
        return ""
    return format_variants_simplified(records)


def compute_family_summary(
    results_dir: Path,
    input_tsv: Path,
    target_samples: dict[str, Any],
    matching_rule: dict[str, Any],
) -> None:
    """Group TNLS target samples into mtDNA-profile families and write a single
    family-level summary that also carries the per-HCLS match statistics.

    One row per family. Matching is aggregated from the existing per-pair
    ``CANNOT_EXCLUDE`` output (`input_tsv`): a family matches an HCLS remain if
    any member is ``CANNOT_EXCLUDE`` against it. The per-HCLS columns
    (``HCLS_matched_variants``, ``HCLS_matched_variant_sites``,
    ``HCLS_matched_variant_bp``, ``HCLS_ratios``) are comma-joined and aligned
    with
    ``Matched_HCLS_IDs``; they are empty for families with no matched HCLS. The
    matched TNLS per HCLS are not repeated — they are the family ``Member_IDs``
    (the exact HCLS->TNLS pairs live in the per-pair TSV), and the member / HCLS
    counts are dropped as derivable from ``Member_IDs`` / ``Matched_HCLS_IDs``.
    This merges the former ``stats_TNLS_HCLS_matching.tsv`` into the family
    summary. The mtDNA matching engine is unchanged.
    """
    minimum = matching_rule["minimum"]
    regions = matching_rule["region"]

    families = group_families(target_samples, regions, minimum)
    logger.info(f"Family grouping: {len(families)} families from {len(target_samples)} TNLS samples")

    df = pd.read_csv(input_tsv, dtype=str, sep="\t")
    if not df.empty:
        df["base_sample"] = df["base_sample"].str.strip()
        df["target_sample"] = df["target_sample"].str.strip()

    hcls_summary = (
        _hcls_match_summary(df)
        if not df.empty
        else pd.DataFrame(
            columns=[
                "base_sample",
                "matched_tnls_list",
                "tnls_count",
                "ratio",
                "matched_variants",
                "matched_variant_sites",
                "matched_variant_bp",
            ],
        )
    )
    summary_by_hcls: dict[str, dict[str, Any]] = {
        row["base_sample"]: row.to_dict()
        for _, row in hcls_summary.iterrows()  # type: ignore[reportGeneralTypeIssues]
    }

    rows: list[dict[str, Any]] = []
    for f_idx, members in enumerate(families):
        family_id = f"FAM-{f_idx + 1:03d}"
        profile_frequency = f"{len(members)}/{len(target_samples)}"
        if df.empty:
            matched_hcls: list[str] = []
        else:
            matched = df.loc[df["target_sample"].isin(members), "base_sample"]
            matched_hcls = sorted(matched.dropna().unique().tolist())

        # Representative = broadest total interval coverage; tie-break by sample id.
        representative = min(
            members,
            key=lambda member: (-sample_coverage_bp(target_samples[member], regions), member),
        )
        rep_sample = target_samples[representative]

        # Per-HCLS match stats, comma-joined and aligned with matched_hcls order.
        # The matched TNLS are the family members (Member_IDs); the exact HCLS
        # ->TNLS pairs live in the per-pair TSV, so they are not repeated here.
        hcls_variants: list[str] = []
        hcls_sites: list[str] = []
        hcls_bp: list[str] = []
        hcls_ratios: list[str] = []
        for hcls in matched_hcls:
            hs = summary_by_hcls.get(hcls)
            if hs is None:
                continue
            hcls_variants.append(str(hs["matched_variants"]))
            hcls_sites.append(str(int(hs["matched_variant_sites"])))
            hcls_bp.append(str(int(hs["matched_variant_bp"])))
            hcls_ratios.append(str(hs["ratio"]))

        hv1_variants = _format_region_variants(rep_sample, "HV1")
        hv2_variants = _format_region_variants(rep_sample, "HV2")
        hv3_variants = _format_region_variants(rep_sample, "HV3")
        profile_variants = " ".join(filter(None, (hv1_variants, hv2_variants, hv3_variants)))

        rows.append(
            {
                "Family_ID": family_id,
                "Member_IDs": ", ".join(members),
                "Profile_Frequency": profile_frequency,
                "Profile_Variants": profile_variants,
                "Matched_HCLS_IDs": ",".join(matched_hcls),
                "HCLS_matched_variants": ", ".join(hcls_variants),
                "HCLS_matched_variant_sites": ", ".join(hcls_sites),
                "HCLS_matched_variant_bp": ", ".join(hcls_bp),
                "HCLS_ratios": ", ".join(hcls_ratios),
                "HV1_variants": hv1_variants,
                "HV2_variants": hv2_variants,
                "HV3_variants": hv3_variants,
            }
        )

    output = results_dir / "stats_TNLS_HCLS_family_summary.tsv"
    pd.DataFrame(rows).to_csv(output, index=False, sep="\t", encoding="utf-8-sig")
    logger.info(f"Family summary saved -> {output}")


def main() -> None:
    script_dir = Path(__file__).parent.resolve()

    # ── CLI args ──────────────────────────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        description="Compare HCLS base samples vs TNLS target samples and compute matching statistics.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=script_dir / "results",
        help="Directory to write all output TSVs.",
    )
    parser.add_argument("--base-samples", required=True, type=Path, help="Path to HCLS base samples JSON")
    parser.add_argument("--target-samples", required=True, type=Path, help="Path to TNLS merged statistics JSON")
    parser.add_argument(
        "--filter-list",
        type=Path,
        default=None,
        help="Path to filter.txt — filters TNLS by LID (optional)",
    )
    parser.add_argument("--matching-rule", required=True, type=Path, help="Path to matching_rule.json")
    parser.add_argument(
        "--tnls-metadata",
        type=Path,
        default=None,
        help="Path to TNLS_metadata.xlsx. Optional; enables identity fields in output.",
    )
    parser.add_argument(
        "--max-ratio",
        type=int,
        default=5,
        help="Maximum 1:N ratio to include in the expanded output (default: 5).",
    )
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()

    # Initialize variables
    column_names = [
        "target_sample",
        "overlap_bp",
        "list_overlap",
        "no_mismatch_pos",
        "mismatch",
        "result",
        "base_sample",
    ]
    final_result = pd.DataFrame(columns=column_names)

    # Load required data
    matching_rule = load_json(args.matching_rule)
    base_samples = load_json(args.base_samples)
    target_raw = load_json(args.target_samples)

    # Apply optional LID filter
    if args.filter_list and args.filter_list.exists():
        with args.filter_list.open(encoding="utf-8") as f:
            lids = {line.strip() for line in f if line.strip()}
        n_raw = len(target_raw)
        target_samples = {k: v for k, v in target_raw.items() if k in lids}
        not_found = lids - set(target_raw.keys())
        logger.info(f"Filter: {len(lids)} LIDs loaded — kept {len(target_samples)}/{n_raw} (skipped: {len(not_found)})")
        if not_found:
            logger.warning(f"Filter LIDs not found in merged statistics: {sorted(not_found)}")
    else:
        target_samples = target_raw
        logger.info(f"No filter applied — using all {len(target_samples)} target samples.")

    n_base = len(base_samples)
    n_target = len(target_samples)
    logger.info(f"Loaded {n_base} HCLS base samples, {n_target} TNLS target samples.")

    # Loop through each base (HCLS) sample
    sample_ids = list(base_samples.keys())
    n_total = len(sample_ids)
    progress_interval = max(1, n_total // 10)  # log every ~10 % of samples
    for i, sample_id in enumerate(sample_ids):
        current_sample = base_samples[sample_id]
        result_json = compare_1n(current_sample, target_samples, matching_rule)
        if len(result_json) > 0:
            result_json["base_sample"] = [sample_id] * len(result_json)
            final_result = pd.concat([final_result, result_json], ignore_index=True)
        if (i + 1) % progress_interval == 0:
            pct = round((i + 1) / n_total * 100)
            logger.info(f"  Progress: {i + 1}/{n_total} HCLS samples ({pct}%)")

    output_path = results_dir / "output_compare_tnls_hcls.tsv"
    results_dir.mkdir(parents=True, exist_ok=True)
    final_result.to_csv(output_path, index=False, sep="\t")
    logger.info(f"Done. {len(final_result)} rows written to {output_path}")

    # ── Compute and write statistics ─────────────────────────────────────────────
    compute_stats(
        results_dir=results_dir,
        input_tsv=output_path,
        base_samples=base_samples,
        target_samples=target_samples,
        filter_list=args.filter_list,
        tnls_metadata=args.tnls_metadata,
        max_ratio=args.max_ratio,
    )

    # ── Family-level summary (mtDNA-profile grouping) ────────────────────────
    compute_family_summary(
        results_dir=results_dir,
        input_tsv=output_path,
        target_samples=target_samples,
        matching_rule=matching_rule,
    )


if __name__ == "__main__":
    main()
