"""Tracy pipeline module — batch orchestration and parallel processing.

Handles batch processing of Tracy AB1 trace files with parallel execution,
following the standard tool module pattern (ETO → pipeline → Batch.write).

The actual ETL (JSON → Sample) lives in ``src.tools.tracy.etl``.
This module owns batch orchestration and parallelism only.

After batch processing, output is organized into three subdirectories
under the base output directory:

- ``json/`` — ``Batch.write()`` output (``statistic_fullbatch.json`` and
  per-sample JSONs)
- ``regions/`` — per-region intermediate JSONs (``HV1_{LID}.json``,
  ``HV2-3_{LID}.json``) written by ``write_region_jsons()``
- ``preprocess/`` — Tracy decompose artifacts (``.align``, ``.decomp``,
  ``.bcf``, ``.json`` files)

The pipeline constructs these paths via ``OutputPaths`` and passes them to
the appropriate write functions. ``Batch.write()`` and
``write_region_jsons()`` signatures are unchanged.
"""

import argparse
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.core.batch import Batch
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Sample, Tool
from src.core.paths import OutputPaths
from src.core.region import REGION_TO_KEYS, validate_region_group_intervals, write_region_jsons
from src.core.sample import merge_intervals
from src.tools.tracy.etl import merge_variants, process
from src.tools.tracy.preprocessing import DecomposeConfig, decompose_sample


@dataclass
class BatchOptions:
    """Configuration options for Tracy batch processing.

    Tracy-specific parameters (trim, pratio, maxindel, quality_threshold,
    min_peak_value, heteroplasmy_threshold) are read from ``get_settings().tracy``
    at usage sites rather than duplicated here.
    """

    ref_path: str = "ref/rCRS.fasta"
    batch_id: str | None = None
    samples_path: Path | None = None


def _read_sample_ids(samples_path: Path) -> list[str]:
    """Read sample IDs from a TXT file (one per line)."""
    with samples_path.open() as f:
        return [line.strip() for line in f if line.strip()]


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

    Uses model-aware deduplication that preserves ``Variant`` objects and
    merges ``files`` and ``quality`` lists.
    """
    if len(samples) == 1:
        return samples[0]

    all_variants: list[Any] = []
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

    # Merge variants across reads — dedup by (pos, ref, seq) and collapse
    # conflicting same-position calls into IUPAC-consensus variants. The
    # variant-merge logic lives in the ETL layer (see
    # :func:`src.tools.tracy.etl.merge_variants`); this pipeline function only
    # aggregates intervals/flags and assembles the combined Sample.
    merged_variants = merge_variants(all_variants)

    # Re-merge overlapping interval spans per region.
    merged_intervals: dict[str, list[list[int]]] = {}
    for region, spans in all_intervals.items():
        if spans:
            merged_intervals[region] = merge_intervals(spans)

    # Re-analyze the merged variant set so cross-region sample-level flags
    # (e.g. 459DEL + 16192T → "Complex: 459 deletion with 16192 issue") are
    # detected. Per-region ETL only sees its own variants. Also recompute
    # variant_flags on merged variants — merge_variants may create IUPAC
    # consensus variants (e.g. 263G + 263A → 263R) not present in any
    # per-region sample, so all_variant_flags would miss their flags.
    if merged_variants:
        flagger = SampleFlagger(merged_variants, intervals=merged_intervals or None)
        _, merged_reasons, _ = flagger.analyze()
        for reason in merged_reasons:
            if reason not in all_sample_flags:
                all_sample_flags.append(reason)
        all_variant_flags = flag_variants(merged_variants)
        all_sample_flags = deduplicate_sample_flags(all_sample_flags, all_variant_flags)

    return Sample(
        sample_id=sample_id,
        variants=merged_variants,
        source_tool=Tool.TRACY,
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


def _write_full_region_jsons(
    sample: Sample,
    sample_id: str,
    regions_dir: Path,
    ref_path: str,
) -> None:
    """Write per-region JSONs for a combined (full-region) Tracy sample.

    Mirrors Sequencher's ``_write_full_region_jsons``: all region groups are
    emitted and ``write_region_jsons`` skips groups without data. Tracy does
    not use auto-pass flags, so region groups are always derived from
    ``REGION_TO_KEYS`` (data-driven) rather than from ``sample_flags``.

    Tracy combines forward/reverse primer traces per LID before this is
    called, so each LID produces one ``HV1_{LID}.json`` and one
    ``HV2-3_{LID}.json`` with merged primer data.

    Args:
        sample: Combined Sample (all primers merged) to write region JSONs for.
        sample_id: Sample identifier (LID) for range-QC logging.
        regions_dir: Directory for region JSON output (``paths.regions_dir``).
        ref_path: Path to rCRS reference FASTA.
    """

    region_groups: dict[str, list[str]] = {
        "HV1": REGION_TO_KEYS["HV1"],
        "HV2-3": REGION_TO_KEYS["HV2-3"],
    }
    write_region_jsons(
        sample=sample,
        output_dir=regions_dir,
        ref_path=ref_path,
        region_groups=region_groups,
    )
    for group_name in region_groups:
        if sample.intervals:
            warnings = validate_region_group_intervals(sample.intervals, group_name)
            for w in warnings:
                logger.warning("Range QC for {} {}: {}", sample_id, group_name, w)


def _process_single_sample(
    sample_id: str,
    input_dir: Path,
    preprocess_dir: Path,
    ref_seq: str,
    ref_path: str,
) -> list[dict[str, Any]]:
    """Process a single Tracy sample: decompose → ETL.

    Args:
        sample_id: Sample identifier (LID).
        input_dir: Directory containing AB1 trace files.
        preprocess_dir: Directory for decompose output (``paths.preprocess_dir``).
        ref_seq: Reference sequence string.
        ref_path: Path to rCRS reference FASTA file.

    Returns:
        List of dicts with sample_id and sample fields (one per region group).
    """
    settings = get_settings()
    tracy_cfg = settings.tracy
    decompose_config = DecomposeConfig(
        trim=tracy_cfg.trim,
        pratio=tracy_cfg.pratio,
        maxindel=tracy_cfg.maxindel,
    )

    # Step 1: Decompose AB1 files — output goes to preprocess_dir/sample_id/
    sample_dir = decompose_sample(
        sample=sample_id,
        input_dir=str(input_dir),
        output_dir=str(preprocess_dir),
        ref_path=ref_path,
        config=decompose_config,
    )

    if not sample_dir or not Path(sample_dir).is_dir():
        logger.error("Decompose failed for sample {}: no output directory", sample_id)
        return []

    # Step 2: Process JSON files from decompose output
    return _process_sample_json_files(sample_id, Path(sample_dir), ref_seq)


def _process_sample_json_files(
    sample_id: str,
    sample_dir: Path,
    ref_seq: str,
) -> list[dict[str, Any]]:
    """Process JSON files from Tracy decompose output for a single sample.

    Args:
        sample_id: Sample identifier (LID).
        sample_dir: Directory containing decompose JSON output.
        ref_seq: Reference sequence string.

    Returns:
        List of dicts with sample_id and sample fields (one per region group).
    """
    json_files = [p for p in sample_dir.glob("*.json") if "HV" in p.name]

    if not json_files:
        logger.warning("No JSON files found for sample {} in {}", sample_id, sample_dir)
        return []

    results: list[dict[str, Any]] = []

    for json_file in json_files:
        try:
            sample = process(
                sample_id=sample_id,
                input_path=json_file,
                ref_seq=ref_seq,
            )
        except (ValueError, KeyError, OSError, RuntimeError) as e:
            logger.error("ETL failed for sample {} file {}: {}", sample_id, json_file.name, e)
            continue

        results.append(
            {
                "sample_id": sample_id,
                "sample": sample,
            },
        )

    return results


def process_batch(
    input_dir: Path,
    output_dir: Path,
    options: BatchOptions | None = None,
) -> dict[str, Sample]:
    """Process Tracy samples in parallel.

    Reads sample IDs from a TXT file, runs Tracy decompose for each sample,
    processes each resulting JSON file through ETL in parallel, combines
    same-LID Samples, and writes output organized into json/, regions/,
    preprocess/ subdirectories.

    Args:
        input_dir: Directory containing AB1 trace files.
        output_dir: Base output directory (subdirectories created automatically).

    Returns:
        Dictionary mapping sample_id to processed Sample.
    """
    if options is None:
        options = BatchOptions()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    # Construct standard output directory structure
    paths = OutputPaths(output_dir)
    paths.create_dirs()

    # Load reference sequence (shared across all samples)
    ref_seq = str(SeqIO.read(options.ref_path, "fasta").seq)

    # Read sample IDs from TXT file (always required)
    if options.samples_path is None:
        logger.error("No samples_path provided — sample IDs must come from a TXT file")
        return {}

    sample_ids = _read_sample_ids(options.samples_path)

    if not sample_ids:
        logger.warning("No sample IDs found in {}", options.samples_path)
        return {}

    logger.info(
        "Loaded {} sample IDs from {}",
        len(sample_ids),
        options.samples_path,
    )

    max_workers = get_settings().tracy.max_workers or max(1, os.cpu_count() or 1)

    # Process each sample in parallel — decompose output goes to preprocess_dir
    processed: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_single_sample,
                sample_id,
                input_dir,
                paths.preprocess_dir,
                ref_seq,
                options.ref_path,
            ): sample_id
            for sample_id in sample_ids
        }
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                results = future.result()
                processed.extend(results)
            except (ValueError, KeyError, OSError, RuntimeError) as e:
                logger.error("Processing failed for sample {}: {}", sample_id, e)

    if not processed:
        logger.warning("No samples successfully processed")
        return {}

    # Group by sample_id, combine same-LID region Samples (merges forward/reverse
    # primer traces per region), then write per-region intermediate JSONs from the
    # combined samples so each LID produces one HV1/HV2-3 JSON with merged data.
    # The base output dir is already batch-specific (results/tools/tracy/<batch_id>),
    # so Batch.write() writes flat into json_dir (nest_batch_id=False) to avoid a
    # duplicated <batch_id> path segment. Output lands directly under json_dir:
    # statistic_fullbatch.json plus per-LID <LID>/<LID>.json.
    combined = _build_combined_samples(processed)
    for sample in combined.values():
        _write_full_region_jsons(
            sample=sample,
            sample_id=sample.sample_id,
            regions_dir=paths.regions_dir,
            ref_path=options.ref_path,
        )
    batch_results = list(combined.values())
    Batch(batch_results).write(
        str(paths.json_dir),
        options.batch_id or "unknown",
        nest_batch_id=False,
    )

    logger.info(
        "Tracy batch processing complete: {} unique LIDs from {} JSON files",
        len(combined),
        len(processed),
    )
    if len(processed) < len(sample_ids):
        failed = len(sample_ids) - len(processed)
        logger.warning("{} sample(s) failed processing", failed)

    return combined


def _build_combined_samples(processed: list[dict[str, Any]]) -> dict[str, Sample]:
    """Group processed results by sample_id and combine same-LID Samples.

    Args:
        processed: List of dicts with sample_id and sample fields.

    Returns:
        Dict mapping sample_id to combined Sample.
    """
    grouped: dict[str, list[Sample]] = defaultdict(list)
    for item in processed:
        grouped[item["sample_id"]].append(item["sample"])

    combined: dict[str, Sample] = {}
    for sample_id, samples in grouped.items():
        combined[sample_id] = _combine_same_tool_samples(samples, sample_id)

    return combined


def main() -> None:
    """CLI entry point for Tracy processing."""

    parser = argparse.ArgumentParser(description="Process Tracy AB1 trace results")
    parser.add_argument("-i", "--input-dir", type=str, required=True, help="Directory containing AB1 trace files")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="Base output dir (subdirs auto-created)")
    parser.add_argument("-s", "--samples", type=str, required=True, help="Path to TXT file with sample IDs")
    parser.add_argument("--batch-id", type=str, default=None, help="Batch identifier")
    parser.add_argument("--ref-path", type=str, default="ref/rCRS.fasta", help="Path to rCRS reference FASTA")

    args = parser.parse_args()

    options = BatchOptions(
        ref_path=args.ref_path,
        batch_id=args.batch_id,
        samples_path=Path(args.samples),
    )

    results = process_batch(
        input_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        options=options,
    )

    logger.info("Processed {} samples", len(results))


if __name__ == "__main__":
    main()
