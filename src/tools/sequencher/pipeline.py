"""Sequencher pipeline module — batch orchestration and parallel processing.

Handles batch processing of Sequencher TXT files with parallel execution,
following the Tracy pipeline pattern (functional API, ProcessPoolExecutor).

The actual ETL (TXT → Sample) lives in ``src.tools.sequencher.etl``.
This module owns batch orchestration and parallelism only.

After batch processing, output is organized under ``OutputPaths``
subdirectories of the base output directory:

- ``json/`` — ``Batch.write()`` output (``statistic_fullbatch.json`` and
  per-sample ``{LID}/{LID}.json``). Because the base dir is batch-specific,
  batch JSONs are written flat (no extra ``<batch_id>`` segment).
- ``regions/`` — per-region intermediate JSONs (``HV1/{LID}``,
  ``HV2-3/{LID}``) used by the unify step for region-level provenance
  tracking. ``--json-dir`` overrides this directory (used by the unify
  workflow, which points ``--seq-dir`` at the same path).

``preprocess/`` holds the QC report (``qc_report.json``) written by the pre-ETL QC gate.

Supports two workflow modes:
  1. One full-region TXT per LID → split into per-region JSONs
  2. Two region-specific TXTs per LID → each produces its own region JSON

In both modes, region JSONs are written per-TXT immediately (no collision),
and same-LID Samples are combined for the batch-level statistic_fullbatch.json.
"""

import argparse
import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from loguru import logger

from src.config import get_settings
from src.core.batch import Batch
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Sample, Tool
from src.core.paths import OutputPaths
from src.core.region import REGION_TO_KEYS, validate_region_group_intervals, write_region_jsons
from src.core.sample import filter_sample_by_regions
from src.core.variants import pos_sort_key
from src.tools.sequencher.etl import process
from src.tools.sequencher.quality_control import reject_duplicate_lids, run_qc
from src.tools.sequencher.utils import detect_region_group_from_filename


@dataclass
class BatchOptions:
    """Configuration options for Sequencher batch processing."""

    ref_path: str = "ref/rCRS.fasta"
    batch_id: str | None = None
    json_dir: Path | None = None


def _merge_intervals(all_intervals: dict[str, list[list[int]]], sample: Sample) -> None:
    """Merge a sample's intervals into the accumulated intervals dict."""
    if not sample.intervals:
        return
    for region_key, spans in sample.intervals.items():
        if region_key in all_intervals:
            all_intervals[region_key].extend(spans)


def _merge_flags(
    all_sample_flags: list[str],
    all_variant_flags: dict[str, list[str]],
    sample: Sample,
) -> None:
    """Merge a sample's flags into the accumulated flag collections."""
    for flag in sample.sample_flags:
        if flag not in all_sample_flags:
            all_sample_flags.append(flag)
    for key, reasons in sample.variant_flags.items():
        if key not in all_variant_flags:
            all_variant_flags[key] = list(reasons)
        else:
            for r in reasons:
                if r not in all_variant_flags[key]:
                    all_variant_flags[key].append(r)


def _combine_same_tool_samples(samples: list[Sample], sample_id: str) -> Sample:
    """Combine region Samples from the same tool into one combined Sample.

    Concatenates variants, merges intervals, and unions flags from region
    Samples of the same tool (Sequencher). This is NOT the cross-tool merge
    that ``unify.py`` owns — it simply assembles a complete view from region
    parts. Since all inputs are same-tool, no conflict resolution is needed.

    Args:
        samples: List of Sample objects for the same LID and tool.
        sample_id: The unified sample ID for all region Samples.

    Returns:
        A combined Sample with all variants, intervals, and flags.
    """
    if len(samples) == 1:
        return samples[0]

    all_variants: list = []
    all_intervals: dict[str, list[list[int]]] = {"HV1": [], "HV2": [], "HV3": []}
    all_sample_flags: list[str] = []
    all_variant_flags: dict[str, list[str]] = {}
    batch_id: str | None = None

    for s in samples:
        all_variants.extend(s.variants)
        _merge_intervals(all_intervals, s)
        _merge_flags(all_sample_flags, all_variant_flags, s)
        if s.batch_id and not batch_id:
            batch_id = s.batch_id

    # Deduplicate variants by position
    pos_map: dict[Any, Any] = {}
    for v in all_variants:
        key = pos_sort_key(v.pos)
        if key not in pos_map:
            pos_map[key] = v
        else:
            # Same position, same tool — merge files and quality
            existing = pos_map[key]
            existing.files = sorted(set(existing.files + v.files))
            existing.quality = list(dict.fromkeys(list(existing.quality + v.quality)))

    merged_intervals = {k: v for k, v in all_intervals.items() if v}

    # Re-analyze the merged variant set so cross-region sample-level flags
    # (e.g. 459DEL + 16192T → "Complex: 459 deletion with 16192 issue") are
    # detected. Per-region ETL only sees its own variants. Recompute
    # variant_flags on the merged set so dedup sees all per-variant flags.
    merged_variant_list = sorted(pos_map.values(), key=lambda v: pos_sort_key(v.pos))
    if merged_variant_list:
        flagger = SampleFlagger(merged_variant_list, intervals=merged_intervals or None)
        _, merged_reasons, _ = flagger.analyze()
        for reason in merged_reasons:
            if reason not in all_sample_flags:
                all_sample_flags.append(reason)
        all_variant_flags = flag_variants(merged_variant_list)
        all_sample_flags = deduplicate_sample_flags(all_sample_flags, all_variant_flags)

    return Sample(
        sample_id=sample_id,
        variants=merged_variant_list,
        source_tool=Tool.SEQUENCHER,
        intervals=merged_intervals or None,
        sample_flags=all_sample_flags or [],
        variant_flags=all_variant_flags,
        information=None,
        batch_id=batch_id,
        hv1=None,
        hv2=None,
        hv3=None,
        no_snps=None,
        no_ins=None,
        no_dels=None,
        no_identicals=None,
        no_unread=None,
    )


def _write_single_region_jsons(
    sample: Sample,
    sample_id: str,
    region_groups: dict[str, list[str]],
    json_output_dir: Path,
    ref_path: str,
) -> None:
    """Write per-region JSONs for a single-region TXT file."""
    region_keys = next(iter(region_groups.values()))
    filtered = filter_sample_by_regions(sample, region_keys)
    if filtered is None:
        return
    write_region_jsons(
        sample=filtered,
        output_dir=json_output_dir,
        ref_path=ref_path,
        region_groups=region_groups,
    )
    for group_name in region_groups:
        if filtered.intervals:
            warnings = validate_region_group_intervals(filtered.intervals, group_name)
            for w in warnings:
                logger.warning("Range QC for {} {}: {}", sample_id, group_name, w)


def _write_full_region_jsons(
    sample: Sample,
    sample_id: str,
    json_output_dir: Path,
    ref_path: str,
) -> None:
    """Write per-region JSONs for a full-region (no prefix) TXT file."""
    region_groups: dict[str, list[str]] = {
        "HV1": REGION_TO_KEYS["HV1"],
        "HV2-3": REGION_TO_KEYS["HV2-3"],
    }
    write_region_jsons(
        sample=sample,
        output_dir=json_output_dir,
        ref_path=ref_path,
        region_groups=region_groups,
    )
    for group_name in region_groups:
        if sample.intervals:
            warnings = validate_region_group_intervals(sample.intervals, group_name)
            for w in warnings:
                logger.warning("Range QC for {} {}: {}", sample_id, group_name, w)


def _write_region_jsons_for_sample(
    sample: Sample,
    sample_id: str,
    filename: str,
    json_output_dir: Path,
    ref_path: str,
) -> None:
    """Dispatch per-region JSON writing based on filename region detection."""
    region_group = detect_region_group_from_filename(filename)

    if region_group is not None:
        region_keys = REGION_TO_KEYS.get(region_group, [])
        if region_keys:
            _write_single_region_jsons(sample, sample_id, {region_group: region_keys}, json_output_dir, ref_path)
    else:
        _write_full_region_jsons(sample, sample_id, json_output_dir, ref_path)


def _build_combined_samples(processed: list[dict[str, Any]]) -> dict[str, Sample]:
    """Group processed samples by LID and combine same-LID region Samples.

    Args:
        processed: List of dicts with 'sample_id' and 'sample' keys.

    Returns:
        Dict mapping sample_id to combined Sample.
    """
    per_lid_samples = defaultdict[str, list[Sample]](list)
    for item in processed:
        per_lid_samples[item["sample_id"]].append(item["sample"])

    combined: dict[str, Sample] = {}
    for sample_id, sample_list in per_lid_samples.items():
        if len(sample_list) == 1:
            combined[sample_id] = sample_list[0]
        else:
            combined[sample_id] = _combine_same_tool_samples(sample_list, sample_id)

    return combined


def process_batch(  # noqa: C901
    input_dir: Path,
    output_dir: Path,
    options: BatchOptions | None = None,
) -> dict[str, Sample]:
    """Process all Sequencher TXT files in a directory in parallel.

    Discovers TXT files, extracts sample IDs from filenames, processes
    each in parallel via ProcessPoolExecutor, writes per-region and
    batch-level JSON output.

    Args:
        input_dir: Directory containing Sequencher TXT files.
        output_dir: Directory for batch-level JSON output.
        options: Batch processing options (ref_path, batch_id, json_dir).
            If None, defaults are used. ``json_dir`` overrides the per-region
            JSON output directory (otherwise ``{output_dir}/regions`` is used).

    Returns:
        Dictionary mapping sample_id to processed Sample.
    """
    if options is None:
        options = BatchOptions()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    # Standard output directory structure (mirrors Tracy/BLASTn via OutputPaths):
    #   {output_dir}/json/      - Batch.write() output (statistic_fullbatch.json + {LID}/{LID}.json)
    #   {output_dir}/regions/   - per-region intermediate JSONs (HV1/{LID}, HV2-3/{LID})
    #   {output_dir}/preprocess/ - QC report (qc_report.json) from the pre-ETL QC gate
    paths = OutputPaths(output_dir)
    paths.json_dir.mkdir(parents=True, exist_ok=True)
    paths.preprocess_dir.mkdir(parents=True, exist_ok=True)
    # Only create the standard regions dir when it is actually used; when
    # ``--json-dir`` overrides it (unify workflow), leave it uncreated.
    if options.json_dir is None:
        paths.regions_dir.mkdir(parents=True, exist_ok=True)

    max_workers = get_settings().sequencher.max_workers or max(1, os.cpu_count() or 1)

    # Discover TXT files and run the pre-ETL QC gate. Samples that fail QC
    # (no LID, missing TXT variant-table header, or invalid analysis ranges) are
    # skipped before ETL — they never enter process(), never get per-region JSON,
    # and never enter the batch JSON.
    # Soft flags (range-missing) do NOT gate JSON; they
    # stay as flags attached by etl.process() for passing samples.
    txt_files = sorted(input_dir.glob("*.TXT"))
    tasks: list[dict[str, Any]] = []
    qc_results: list[dict[str, Any]] = []

    for txt_file in txt_files:
        qc = run_qc(txt_file.name, txt_file)
        qc_results.append(
            {
                "sample_id": qc.sample_id,
                "filename": txt_file.name,
                "status": "pass" if qc.passed else "fail",
                "reason": qc.reason,
            },
        )

    reject_duplicate_lids(qc_results)

    for txt_file, qc_result in zip(txt_files, qc_results, strict=True):
        if qc_result["status"] != "pass":
            logger.warning("QC failed for {}: {}", txt_file.name, qc_result["reason"])
            continue

        tasks.append(
            {
                "sample_id": qc_result["sample_id"],
                "input_path": txt_file,
            },
        )

    # Write the full QC report (pass + fail) immediately after discovery, before
    # processing, so it always reflects every discovered TXT regardless of any
    # later per-sample processing failure.
    passed_count = sum(1 for r in qc_results if r["status"] == "pass")
    qc_report = {
        "tool": "sequencher",
        "total": len(qc_results),
        "passed": passed_count,
        "failed": len(qc_results) - passed_count,
        "samples": qc_results,
    }
    qc_report_path = paths.preprocess_dir / "qc_report.json"
    qc_report_path.write_text(json.dumps(qc_report, indent=2), encoding="utf-8")
    logger.info(
        "QC report written to {} ({} passed, {} failed of {} discovered)",
        qc_report_path,
        passed_count,
        len(qc_results) - passed_count,
        len(qc_results),
    )

    logger.info("Found {} TXT files in {} ({} passed QC)", len(txt_files), input_dir, len(tasks))

    if not tasks:
        logger.warning("No valid samples to process in {}", input_dir)
        return {}

    # Process samples in parallel
    processed: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process,
                task["sample_id"],
                task["input_path"],
                options.batch_id,
            ): task
            for task in tasks
        }

        for future in as_completed(futures):
            task = futures[future]
            sample_id = task["sample_id"]
            try:
                sample = future.result()
            except (ValueError, KeyError, OSError, RuntimeError) as e:
                logger.error("Processing failed for {}: {}", sample_id, e)
                continue

            processed.append(
                {
                    "sample_id": sample_id,
                    "input_path": task["input_path"],
                    "sample": sample,
                },
            )

    if not processed:
        return {}

    # Write per-region intermediate JSONs per-TXT (immediate, no collision).
    # ``--json-dir`` overrides the regions directory (used by the unify workflow,
    # which points ``--seq-dir`` at the same path); otherwise region JSONs go to
    # the standard ``paths.regions_dir``.
    regions_output_dir = options.json_dir if options.json_dir is not None else paths.regions_dir
    for item in processed:
        _write_region_jsons_for_sample(
            sample=item["sample"],
            sample_id=item["sample_id"],
            filename=item["input_path"].name,
            json_output_dir=regions_output_dir,
            ref_path=options.ref_path,
        )

    # Group by sample_id, combine same-LID region Samples, write batch JSON.
    # The base output dir is already batch-specific (results/tools/sequencher/<batch_id>),
    # so Batch.write() writes flat into json_dir (nest_batch_id=False) to avoid a
    # duplicated <batch_id> path segment.
    combined = _build_combined_samples(processed)
    batch_results = list(combined.values())
    Batch(batch_results).write(
        str(paths.json_dir),
        options.batch_id or "unknown",
        nest_batch_id=False,
    )

    logger.info(
        "Sequencher batch processing complete: {} unique LIDs from {} TXT files",
        len(combined),
        len(tasks),
    )
    if len(processed) < len(tasks):
        failed = len(tasks) - len(processed)
        logger.warning("{} TXT file(s) failed processing", failed)
    return combined


def main() -> None:
    """CLI entry point for Sequencher batch processing."""

    parser = argparse.ArgumentParser(description="Sequencher batch processing")
    parser.add_argument("--input-dir", type=Path, required=True, help="Directory containing Sequencher TXT files")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for batch-level JSON output")
    parser.add_argument("--ref-path", default="ref/rCRS.fasta", help="Path to rCRS reference FASTA")
    parser.add_argument("--batch-id", default=None, help="Batch identifier")
    parser.add_argument("--json-dir", type=Path, default=None, help="Directory for per-region intermediate JSONs")

    args = parser.parse_args()

    options = BatchOptions(
        ref_path=args.ref_path,
        batch_id=args.batch_id,
        json_dir=args.json_dir,
    )

    process_batch(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        options=options,
    )


if __name__ == "__main__":
    main()
