"""Mutation Surveyor pipeline module — ETL + auto-pass region filter + AB1 sort.

Runs the ETL stage: reads the combined final profiles xlsx, processes each
LID into a Sample, applies auto-pass region filtering, and writes batch JSON.
Then sorts AB1 trace files by autopass/review classification.

Steps 01 to 06 (trace processing → control QC → HV2-3 merge → HV1 merge →
final profiles → final profiles vs truth) are orchestrated by the shell
script ``scripts/tools/mutation_surveyor.sh``. This module handles step 08
(ETL + auto-pass filter) and step 09 (AB1 sort by classification) as a
single CLI entry point.

After ETL builds a Sample, auto-pass region filtering is applied inline
(via ``AUTOPASS_FLAG_TO_REGIONS`` and ``autopass_regions_from_flags()``):
samples are filtered to keep only auto-pass regions, and samples with no
auto-pass regions are excluded from the batch output.

Per-region intermediate JSONs (HV1_{LID}.json, HV2-3_{LID}.json) are written
directly under json_dir (defaults to output_dir). These are temporary
intermediate files used by the unify step for region-level provenance tracking.
The final unified JSON per LID follows the standard format defined in ARCHITECTURE.md.

Entry points:
- ``process_batch()`` — ETL + auto-pass filter + per-region JSON output from final_profiles xlsx.
- ``main()`` — CLI entry point for shell script step 08.
"""

import argparse
from pathlib import Path

import pandas as pd
from loguru import logger

from src.config import get_settings
from src.core.batch import Batch
from src.core.models import Sample
from src.core.region import region_groups_from_flags, validate_region_group_intervals, write_region_jsons
from src.core.sample import filter_sample_by_regions
from src.tools.mutation_surveyor.etl import process
from src.tools.mutation_surveyor.sort_review_data import sort_review_data
from src.tools.mutation_surveyor.utils import safe_str

# ---------------------------------------------------------------------------
# Auto-pass flag → region mapping (previously in filter_regions.py)
# ---------------------------------------------------------------------------

AUTOPASS_FLAG_TO_REGIONS: dict[str, list[str]] = {
    "Autopass HV1": ["HV1"],
    "Autopass HV2 and HV3": ["HV2", "HV3"],
}


def autopass_regions_from_flags(sample_flags: list[str]) -> set[str]:
    """Determine which regions are auto-pass from sample_flags.

    Returns:
        Set of region names (e.g. {"HV1"}, {"HV2", "HV3"}) that are auto-pass.
    """
    kept: set[str] = set()
    for flag in sample_flags:
        regions = AUTOPASS_FLAG_TO_REGIONS.get(flag)
        if regions:
            kept.update(regions)
    return kept


def _filter_sample_by_autopass(sample: Sample) -> Sample | None:
    """Filter a Sample to keep only auto-pass regions.

    Returns:
        Filtered Sample with only auto-pass region variants, or None if
        the sample has no auto-pass regions or no variants in those regions.
    """
    keep = autopass_regions_from_flags(sample.sample_flags)
    if not keep:
        return None

    return filter_sample_by_regions(sample, sorted(keep))


def _process_row(
    row: pd.Series,
    batch_id: str,
    final_profiles: Path,
    ref_path: str,
    region_json_dir: Path,
) -> tuple[str, Sample] | None:
    """Process a single row: ETL, auto-pass filter, and per-region JSON write.

    Returns:
        Tuple of (sample_id, filtered_sample) if the sample passed auto-pass
        filtering, or None if it should be skipped.
    """
    sample_id = safe_str(row.get("LID", ""))
    if not sample_id or sample_id.upper().startswith(("PC", "NTC")):
        return None

    try:
        sample = process(
            sample_id=sample_id,
            input_path=final_profiles,
            ref_path=ref_path,
            batch_id=batch_id,
        )
    except (ValueError, KeyError, OSError, RuntimeError) as e:
        logger.error("ETL failed for LID {}: {}", sample_id, e)
        return None

    # Auto-pass region filter
    filtered = _filter_sample_by_autopass(sample)
    if filtered is None or not filtered.variants:
        return None

    # Write per-region intermediate JSONs for auto-pass regions
    region_groups = region_groups_from_flags(sample.sample_flags)
    if region_groups:
        write_region_jsons(
            sample=filtered,
            output_dir=region_json_dir,
            ref_path=ref_path,
            region_groups=region_groups,
        )
        # Validate intervals for each region group (warning flags, not hard error)
        for group_name in region_groups:
            if filtered.intervals:
                warnings = validate_region_group_intervals(filtered.intervals, group_name)
                for w in warnings:
                    logger.warning("Range QC for {} {}: {}", sample_id, group_name, w)

    return sample_id, filtered


def process_batch(
    final_profiles: Path,
    output_dir: Path,
    ref_path: str = "ref/rCRS.fasta",
    batch_id: str | None = None,
    json_dir: Path | None = None,
) -> dict[str, Sample]:
    """Run ETL + auto-pass region filter from final_profiles xlsx.

    Reads the LID sheet from final_profiles, processes each LID into a
    Sample via ETL, applies auto-pass region filtering, and writes batch
    JSON (per-sample and statistic_fullbatch) via Batch.write().

    Samples without auto-pass regions are excluded from the batch output.

    Args:
        final_profiles: Path to final_profiles xlsx (from step 5).
        output_dir: Output directory for batch JSON.
        ref_path: Path to rCRS reference FASTA.
        batch_id: Optional batch identifier (derived from filename if not provided).
        json_dir: Optional directory for per-region JSON output. If None, uses output_dir.

    Returns:
        Dictionary mapping sample_id → Sample for all processed LIDs
        (only samples with auto-pass regions are included).
    """
    final_profiles = Path(final_profiles)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if batch_id is None:
        batch_id = final_profiles.stem
        for suffix in ("_final_profiles",):
            batch_id = batch_id.replace(suffix, "")

    region_json_dir = json_dir if json_dir is not None else output_dir

    # Read LID sheet from final_profiles
    logger.info("Reading final profiles from {}", final_profiles)
    df = pd.read_excel(final_profiles, sheet_name="LID")
    df.columns = [str(c).strip() for c in df.columns]

    results: dict[str, Sample] = {}
    skipped = 0
    total_rows = len(df)

    for _, row in df.iterrows():
        result = _process_row(row, batch_id, final_profiles, ref_path, region_json_dir)
        if result is None:
            skipped += 1
            continue
        sample_id, filtered = result
        results[sample_id] = filtered

    # Write batch JSON (only filtered samples)
    if results:
        Batch(list(results.values()), ref_path=ref_path).write(str(output_dir), batch_id)
        logger.info("Batch write: {} samples ({} skipped of {}) → {}", len(results), skipped, total_rows, output_dir)
    else:
        logger.warning("No samples to write from {}", final_profiles)

    logger.success("Pipeline complete — {} LIDs processed", len(results))

    return results


def main() -> None:
    """CLI entry point: ETL + auto-pass filter + AB1 sort by classification."""
    parser = argparse.ArgumentParser(
        description="ETL + auto-pass filter + AB1 sort by classification",
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        required=True,
        help="Path to final_profiles xlsx (from step 5)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for batch JSON",
    )
    parser.add_argument(
        "--batch-id",
        default=None,
        help="Batch identifier (derived from filename if not provided)",
    )
    parser.add_argument(
        "--ref-path",
        default="ref/rCRS.fasta",
        help="Path to rCRS reference FASTA (default: ref/rCRS.fasta)",
    )
    parser.add_argument(
        "--json-dir",
        type=Path,
        default=None,
        help="Output directory for per-region JSON files (default: same as --output-dir)",
    )

    args = parser.parse_args()

    # Step 1: ETL + auto-pass region filter
    process_batch(
        final_profiles=args.input_dir,
        output_dir=args.output_dir,
        ref_path=args.ref_path,
        batch_id=args.batch_id,
        json_dir=args.json_dir,
    )

    # Derive batch_id from input_dir if not provided
    batch_id = args.batch_id
    if batch_id is None:
        batch_id = Path(args.input_dir).stem
        for suffix in ("_final_profiles",):
            batch_id = batch_id.replace(suffix, "")

    # Step 2: AB1 sort by autopass/review classification
    excels_dir = Path(args.input_dir).parent

    # Derive raw_dir from settings (was --raw-dir, now config-driven)
    settings = get_settings()
    raw_dir = settings.directories.data / "raw" / batch_id
    data_dir = "/mnt/nas/bca/mtDNA/science/data"
    raw_dir = Path(data_dir, "raw", batch_id)

    sort_review_data(
        batch_id=batch_id,
        results_dir=excels_dir,
        raw_dir=raw_dir,
        excels_dir=excels_dir,
    )


if __name__ == "__main__":
    main()
