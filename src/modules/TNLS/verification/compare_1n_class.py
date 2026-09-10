#!/usr/bin/env python
"""Three-stage HCLS-vs-TNLS matching engine.

One HCLS sample (base) is compared against all TNLS samples (targets) in
three sequential filters.  Stages 1 and 2 are fast pre-filters that eliminate
obviously non-matching pairs before the more expensive base-level comparison
in Stage 3.

Thresholds (minimum overlap, mismatch limits) are read from
``matching_rule.json`` at construction time.
"""

from typing import Any

import pandas as pd

from src.core.variants import get_overlap
from src.modules.TNLS.verification.utils import (
    combined_variants,
    compare_variants,
    count_position_mismatches,
    format_matched_variants,
    format_mismatch_variants,
    format_sample_variants,
    get_hcls_unknown_snp_positions,
    get_range_details,
    get_variant_positions,
)


def compare_1n(
    base_sample: dict[str, Any], target_samples: dict[str, Any], matching_rule: dict[str, Any]
) -> pd.DataFrame:
    """Compare one HCLS base sample against all TNLS target samples.

    This is the only public entry point.  Returns a DataFrame of
    ``CANNOT_EXCLUDE`` pairs only (all other classifications are filtered out
    by Stage 3).
    """
    return CompareManager(base_sample, target_samples, matching_rule).run()


class CompareManager:
    """Manages one base sample vs. all targets through the 3-stage pipeline."""

    def __init__(
        self, base_sample: dict[str, Any], target_samples: dict[str, Any], matching_rule: dict[str, Any]
    ) -> None:
        # Unwrap double-nested JSON wrapper if present (sample name is the top key).
        self.base = base_sample
        if "variants" not in self.base:
            self.base = self.base[next(iter(self.base))]

        self.target = target_samples  # dict: sample_name -> sample_data
        self.MR_limit = matching_rule["minimum"]
        self.MR_region = matching_rule["region"]
        self.MR_conclusion_threshold = matching_rule["conclusion_threshold"]
        self.MR_exclusion_threshold = matching_rule["exclusion_threshold"]

    def run(self) -> pd.DataFrame:
        """Convenience entry point wrapping the three stage methods."""
        return self._finish(self.s3_filter_variants_based(self.s2_filter_variants_based(self.s1_get_overlap_list())))

    # ── Internal helpers ─────────────────────────────────────────────────────

    def _get_target_combined_vars(self, sample_name: str) -> dict[str, Any]:
        """Build a typed-combined variant dict for one TNLS sample (used in Stage 3)."""
        return combined_variants(self.target[sample_name])

    def _get_target_variants_list(self) -> list[dict[str, Any]]:
        """Extract the raw variants dict from every TNLS sample, as a list.

        The list index corresponds to the DataFrame row index, so
        ``target_variants_list[i]`` is the variants dict for
        ``df.loc[i, "target_sample"]``.
        """
        return [data["variants"] for data in self.target.values() if isinstance(data, dict)]

    def _get_base_region(self) -> list[list[int]]:
        """Flatten all HV intervals from the base sample into a single list.

        ``get_overlap`` takes two interval lists and computes pairwise intersections.
        We pass the base intervals once and each target interval in turn.
        """
        intervals = []
        for region in self.MR_region:
            # A base sample may only have been sequenced for a subset of regions
            # (e.g. HCLS samples covering only HV1); missing region keys carry
            # no intervals. Use .get() to mirror the target-side handling.
            intervals.extend(self.base["intervals"].get(region, []))
        return intervals

    def _get_target_regions(self) -> list[list[list[int]]]:
        """Return per-target interval lists aligned with ``self.target.values()``.

        Same layout as ``_get_target_variants_list``: the i-th interval list
        belongs to the i-th target sample name.
        """
        return [
            [
                interval
                for region in self.MR_region
                for interval in (data["intervals"].get(region, []) if isinstance(data, dict) else [])
            ]
            for data in self.target.values()
        ]

    def _extract_hv_regions(self, overlap_details: list[str]) -> list[str]:
        """Strip the "HV1: 16024-16365" strings down to just "HV1"."""
        return [detail.split(":")[0] for detail in overlap_details if detail.startswith("HV")]

    def _finish(self, df: pd.DataFrame) -> pd.DataFrame:
        """Post-process Stage 3 output: add readable columns and filter."""
        # ``s3_filter_variants_based`` already calls ``_finish`` on the
        # CANNOT_EXCLUDE subset; ``run`` wraps that result in ``_finish`` again.
        # Guard against the double call so the simplified ``variant_match``
        # string is not reprocessed as a dict.
        if "matched_variants" in df.columns:
            return df
        # Add human-readable HV region labels for the overlap intervals.
        df["overlap_details"] = df["list_overlap"].apply(get_range_details)
        df["hv_region_overlap"] = df["overlap_details"].apply(self._extract_hv_regions)

        # Count matched variants and bases from the variant_match dict.
        def count_bases(d: dict[str, Any]) -> int:
            return sum(len(str(v.get("seq", ""))) for v in d.values())

        df["matched_variants"] = df["variant_match"].apply(len)
        df["matched_bases"] = df["variant_match"].apply(count_bases)

        # Human-readable simplified variant columns.
        df["HCLS_variants"] = format_sample_variants(self.base)
        df["TNLS_variants"] = df["target_sample"].apply(lambda name: format_sample_variants(self.target[name]))

        # Simplified display for the raw match/mismatch dicts as well. The matched
        # variant list lives in ``variant_match`` (the former separate
        # ``matched_variants`` string column was redundant with it).
        df["variant_match"] = df["variant_match"].apply(format_matched_variants)
        df["variant_mismatch"] = df["variant_mismatch"].apply(format_mismatch_variants)
        return df

    # ── Stage 1: Interval overlap ─────────────────────────────────────────────────

    def s1_get_overlap_list(self) -> pd.DataFrame:
        """Compute genomic interval overlap between the base sample and every target.

        For each TNLS sample, ``get_overlap`` finds the intersecting [start, end]
        intervals across all configured HV regions and sums their total base-pair
        length.  This is the first gate: if the overlapping region is shorter
        than ``MR_limit`` (150 bp), the pair is not worth comparing further.

        Output columns:
            target_sample   TNLS sample name
            overlap_bp    total overlapping base pairs
            list_overlap    list of [start, end] intervals in the overlap
        """
        base_intervals = self._get_base_region()
        target_intervals = self._get_target_regions()

        overlaps, overlap_bases = [], []
        for t_intervals in target_intervals:
            ov = get_overlap(base_intervals, t_intervals)
            overlaps.append(ov)
            # Inclusive 1-based interval length: [s, e] spans (e - s + 1) bp.
            overlap_bases.append(sum(e - s + 1 for s, e in ov))

        return pd.DataFrame(
            {
                "target_sample": list(self.target.keys()),
                "overlap_bp": overlap_bases,
                "list_overlap": overlaps,
            },
        )

    # ── Stage 2: Position-level fast filter ─────────────────────────────────────

    def s2_filter_variants_based(self, df: pd.DataFrame) -> pd.DataFrame:
        """Count variant-position mismatches (ignoring base values) within the overlap.

        This is a fast pre-filter before the expensive base-level comparison
        in Stage 3.  It looks only at *which positions* have variants, ignoring
        whether the actual bases match.

        If a sample has a variant at a position where the other sample does not,
        that counts as a position mismatch. HCLS SNP calls of ``N`` are removed
        from both position sets because that base is unresolved and must not reject
        the pair. Pairs with ``n_mismatch >= MR_exclusion_threshold`` (i.e. 2 or
        more position mismatches) are fast-tracked to ``EXCLUSION`` in Stage 3
        without doing any base-level comparison.

        Output column:
            no_mismatch_pos   count of positions that differ between samples
                             (-1 means overlap was below MR_limit — not evaluated)
        """
        base_variants = self.base["variants"]
        target_variants_list = self._get_target_variants_list()
        mismatch_counts = []

        for i, row in df.iterrows():
            if row["overlap_bp"] < self.MR_limit:
                # Stage 1 already disqualified this pair; mark as not evaluated.
                mismatch_counts.append(-1)
                continue

            target_variants = target_variants_list[i]  # type: ignore[reportCallIssue]
            cur_overlap = row["list_overlap"]

            # Ignore unresolved HCLS N calls on both sides of the position check.
            ignored_positions = get_hcls_unknown_snp_positions(
                base_variants,
                cur_overlap,  # type: ignore[reportArgumentType]
            )
            base_positions = get_variant_positions(base_variants, cur_overlap)  # type: ignore[reportArgumentType]
            target_positions = get_variant_positions(target_variants, cur_overlap)  # type: ignore[reportArgumentType]
            mismatch_counts.append(count_position_mismatches(base_positions, target_positions, ignored_positions))

        df["no_mismatch_pos"] = mismatch_counts
        return df

    # ── Stage 3: Base-level comparison ────────────────────────────────────────────

    def s3_filter_variants_based(self, df: pd.DataFrame) -> pd.DataFrame:
        """Classify each pair as CANNOT_EXCLUDE / INCONCLUSIVE / EXCLUSION.

        Only pairs that passed both Stage 1 (sufficient overlap) and Stage 2
        (fewer than MR_exclusion_threshold position mismatches) reach this stage.

        For those pairs, the full variant comparison runs:

        Classification thresholds (from ``matching_rule.json``):
            mismatches == 0 and matched variants > 0 -> CANNOT_EXCLUDE
            mismatches == 0 and matched variants == 0 -> INCONCLUSIVE
            mismatches == 1                          -> INCONCLUSIVE
            mismatches >= 2                          -> EXCLUSION

        Pairs with ``no_mismatch_pos >= MR_exclusion_threshold`` (i.e. 2+)
        are marked EXCLUSION without base-level comparison — they already failed
        the fast Stage 2 filter.

        Output columns:
            mismatch        count of base-level mismatches
            result          CANNOT_EXCLUDE | INCONCLUSIVE | EXCLUSION |
                           MR_MINIMUM_NOT_REACH
            variant_match   dict of matched variant keys -> {"ref", "seq"}
            variant_mismatch dict of mismatched keys -> {"base_sample", "target_sample"}

        Only ``CANNOT_EXCLUDE`` rows are returned.
        """
        base_combined = combined_variants(self.base)
        mismatch_counts, results, match_dicts, mismatch_dicts = [], [], [], []

        for _, row in df.iterrows():
            # ── Stage 1 gate ────────────────────────────────────────────────────
            if row["overlap_bp"] < self.MR_limit:
                results.append("MR_MINIMUM_NOT_REACH")
                mismatch_counts.append("-")
                match_dicts.append({})
                mismatch_dicts.append({})
                continue

            target_name = row["target_sample"]
            cur_overlap = row["list_overlap"]

            # ── Stage 2 gate: fast exclusion if too many position mismatches ───
            if row["no_mismatch_pos"] < self.MR_exclusion_threshold:
                target_combined = self._get_target_combined_vars(target_name)  # type: ignore[reportArgumentType]
                n_match, n_mismatch, mdict, mmdict = compare_variants(base_combined, target_combined, cur_overlap)  # type: ignore[reportArgumentType]
                mismatch_counts.append(n_mismatch)
                match_dicts.append(mdict)
                mismatch_dicts.append(mmdict)

                # A no-mismatch pair is only actionable when at least one variant
                # matches inside the overlap. Without matched variants, the pair
                # has no positive mtDNA evidence and must not be reported as
                # CANNOT_EXCLUDE.
                if n_mismatch <= self.MR_conclusion_threshold and n_match > 0:
                    results.append("CANNOT_EXCLUDE")
                elif n_mismatch < self.MR_exclusion_threshold:
                    results.append("INCONCLUSIVE")
                else:
                    results.append("EXCLUSION")
            else:
                # Stage 2 failed: ≥2 position mismatches — exclude without base check.
                mismatch_counts.append(-1)
                results.append("EXCLUSION")
                match_dicts.append({})
                mismatch_dicts.append({})

        df["mismatch"] = mismatch_counts
        df["result"] = results
        df["variant_match"] = match_dicts
        df["variant_mismatch"] = mismatch_dicts

        # Return only CANNOT_EXCLUDE pairs; post-process for readable columns.
        return self._finish(df[df["result"] == "CANNOT_EXCLUDE"])  # type: ignore[reportArgumentType]
