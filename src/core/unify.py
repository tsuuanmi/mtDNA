"""Batch unification: merge MS and Sequencher per-region JSONs into unified profiles.

Discovers per-region JSON files from both Mutation Surveyor and Sequencher
directories. Merges by region via ``merge_by_regions()`` with region_sources
provenance tracking.

Uses the ``regions/{group_name}/`` subdirectory layout under each tool
directory. Per-region JSONs are expected at
``{tool_dir}/regions/HV1/HV1_{LID}.json`` and
``{tool_dir}/regions/HV2-3/HV2-3_{LID}.json``.

For each sample with per-region JSONs:
  - Discovers HV1_{LID}.json and HV2-3_{LID}.json in each tool directory
  - Groups by region, merges by region with merge_by_regions()
  - Produces 1 unified {LID}.json per LID with region_sources in information

Final output: 1 unified JSON per LID in the standard format, plus
statistic_fullbatch.json.
"""

import json
from argparse import ArgumentParser
from pathlib import Path

from loguru import logger

from src.config import get_settings
from src.core.batch import Batch
from src.core.models import Sample
from src.core.mtdna_merger import MergeValidationError, merge_by_regions
from src.core.region import REGION_FILE_PATTERNS, read_region_json


def _resolve_region_dir(tool_dir: Path, group_name: str) -> Path | None:
    """Resolve the region directory for a tool/group.

    Returns ``{tool_dir}/regions/{group_name}`` when it exists, else ``None``.
    """
    region_dir = tool_dir / "regions" / group_name
    return region_dir if region_dir.is_dir() else None


def _extract_sample_id(filepath: Path, group_name: str) -> str | None:
    """Extract the sample_id from a region filename, or ``None`` if it doesn't match.

    Removes the group prefix and ``.json`` suffix, e.g.
    ``HV1_LN_26_AB3688.json`` -> ``LN_26_AB3688``.
    """
    stem = filepath.name.removesuffix(".json")
    prefix = group_name + "_"
    if not stem.startswith(prefix):
        return None
    return stem[len(prefix) :]


def _discover_region_samples(
    ms_dir: Path,
    seq_dir: Path,
) -> dict[str, dict[str, dict[str, Path]]]:
    """Discover per-region JSONs for all samples across MS and Sequencher dirs.

    Expects the ``regions/{group_name}/`` subdirectory layout under each
    tool directory:
        {tool_dir}/regions/HV1/HV1_{LID}.json
        {tool_dir}/regions/HV2-3/HV2-3_{LID}.json

    Returns:
        Dict mapping sample_id → group_name → tool_name → file_path.
        E.g. {"SAMPLE_001": {"HV1": {"mutation_surveyor": path, "sequencher": path}, ...}}
    """
    result: dict[str, dict[str, dict[str, Path]]] = {}

    for tool_name, tool_dir in [("mutation_surveyor", ms_dir), ("sequencher", seq_dir)]:
        if not tool_dir.is_dir():
            continue

        for group_name, pattern in REGION_FILE_PATTERNS.items():
            region_dir = _resolve_region_dir(tool_dir, group_name)
            if region_dir is None:
                continue

            for filepath in sorted(region_dir.glob(pattern.format(lid="*"))):
                if not filepath.is_file():
                    continue
                sample_id = _extract_sample_id(filepath, group_name)
                if sample_id is None:
                    continue  # skip files that don't match the expected pattern
                result.setdefault(sample_id, {}).setdefault(group_name, {})[tool_name] = filepath

    return result


def _process_per_region_samples(
    per_region_sample_ids: set[str],
    region_samples: dict[str, dict[str, dict[str, Path]]],
    ref_path: Path | None,
    unified_out: Path,
) -> list[Sample]:
    """Merge per-region samples and return the unified samples."""
    unified: list[Sample] = []

    for sample_id in sorted(per_region_sample_ids):
        region_data = region_samples[sample_id]

        # Collect region-group → list of Samples
        region_groups: dict[str, list[Sample]] = {}

        for group_name, tool_paths in region_data.items():
            for tool_name, filepath in tool_paths.items():
                try:
                    sample = read_region_json(filepath, sample_id)
                except (ValueError, KeyError, OSError) as e:
                    logger.error("Failed to load {} {} for {}: {}", group_name, tool_name, sample_id, e)
                    continue

                if group_name not in region_groups:
                    region_groups[group_name] = []
                region_groups[group_name].append(sample)

        if not region_groups:
            logger.warning("No region data loaded for sample {}, skipping", sample_id)
            continue

        try:
            result = merge_by_regions(
                region_samples=region_groups,
                sample_id=sample_id,
                ref_path=ref_path,
                output_dir=unified_out,
            )
        except MergeValidationError as e:
            logger.error("Merge failed for sample {}: {}", sample_id, e)
            raise

        unified.append(result.sample)
        if result.warnings:
            for w in result.warnings:
                logger.warning("Region QC warning for {}: {}", sample_id, w)

    return unified


def unify_batch(
    ms_dir: Path,
    seq_dir: Path,
    output_dir: Path,
    batch_id: str,
    ref_path: Path | None = None,
) -> int:
    """Merge MS and Sequencher per-region JSONs into unified profiles.

    Discovers per-region JSONs and merges them via merge_by_regions().

    Args:
        ms_dir: Directory with Mutation Surveyor per-region JSONs.
        seq_dir: Directory with Sequencher per-region JSONs.
        output_dir: Base output directory for unified results.
        batch_id: Batch identifier.
        ref_path: Path to rCRS reference FASTA.

    Returns:
        Number of unified samples produced.
    """
    settings = get_settings()
    if ref_path is None:
        ref_path = settings.directories.ref / "rCRS.fasta"

    # Step 1: Set up output directory
    unified_out = output_dir
    unified_out.mkdir(parents=True, exist_ok=True)

    # Step 2: Discover per-region samples
    region_samples = _discover_region_samples(ms_dir, seq_dir)

    # Step 3: Process per-region samples
    all_unified_samples = _process_per_region_samples(set(region_samples.keys()), region_samples, ref_path, unified_out)

    # Step 4: Write statistic_fullbatch.json
    if all_unified_samples:
        batch = Batch(all_unified_samples, ref_path=str(ref_path))
        batch_data = batch.to_json()
        batch_json_path = unified_out / "statistic_fullbatch.json"
        batch_json_path.parent.mkdir(parents=True, exist_ok=True)
        with batch_json_path.open("w", encoding="utf-8") as fh:
            json.dump(batch_data, fh, indent=2)
        logger.success(
            "Wrote batch JSON for {} ({} samples) to {}",
            batch_id,
            len(all_unified_samples),
            batch_json_path,
        )

    return len(all_unified_samples)


def main() -> None:
    """CLI entry point for batch unification."""
    parser = ArgumentParser(description="Merge MS and Sequencher per-region JSONs into unified profiles.")
    parser.add_argument(
        "--ms-dir",
        type=Path,
        required=True,
        help="Directory with Mutation Surveyor per-region JSONs",
    )
    parser.add_argument(
        "--seq-dir",
        type=Path,
        required=True,
        help="Directory with Sequencher per-region JSONs",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Base output directory for unified results",
    )
    parser.add_argument(
        "--batch-id",
        type=str,
        required=True,
        help="Batch identifier",
    )
    parser.add_argument(
        "--ref-path",
        type=Path,
        default=None,
        help="Path to rCRS.fasta (overrides default from settings)",
    )
    args = parser.parse_args()

    count = unify_batch(
        ms_dir=args.ms_dir,
        seq_dir=args.seq_dir,
        output_dir=args.output_dir,
        batch_id=args.batch_id,
        ref_path=args.ref_path,
    )
    logger.info("Unified {} samples", count)


if __name__ == "__main__":
    main()
