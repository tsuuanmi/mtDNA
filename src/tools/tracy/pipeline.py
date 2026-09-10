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
- ``preprocess/<sample_id>/`` — Tracy artifacts and per-sample ``qc_report.json``

The pipeline constructs these paths via ``OutputPaths`` and passes them to
the appropriate write functions. ``Batch.write()`` and
``write_region_jsons()`` signatures are unchanged.
"""

import argparse
import json
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
from src.tools.tracy.noise_mask import apply_noise_mask
from src.tools.tracy.preprocessing import DecomposeConfig, decompose_sample
from src.tools.tracy.quality_control import (
    ExcludedVariant,
    NoiseConfig,
    TraceQCResult,
    analyze_trace_file,
    build_qc_report,
    error_result,
    unavailable_result,
    with_excluded_variants,
)


@dataclass(frozen=True)
class PreparedSample:
    """Decomposed Tracy inputs and their QC results for one sample."""

    sample_id: str
    json_paths: tuple[Path, ...]
    qc_results: tuple[TraceQCResult, ...]


@dataclass(frozen=True)
class ProcessedTrace:
    """Post-ETL trace sample and its noise-mask exclusions."""

    sample_id: str
    filename: str
    sample: Sample
    excluded_variants: tuple[ExcludedVariant, ...]


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


def _prepare_sample(
    sample_id: str,
    input_dir: Path,
    preprocess_dir: Path,
    ref_path: str,
    noise_config: NoiseConfig,
) -> PreparedSample:
    """Run Tracy decompose and QC for one sample without calling variants."""
    tracy_config = get_settings().tracy
    json_outputs = decompose_sample(
        sample=sample_id,
        input_dir=str(input_dir),
        output_dir=str(preprocess_dir),
        ref_path=ref_path,
        config=DecomposeConfig(
            trim=tracy_config.trim,
            pratio=tracy_config.pratio,
            maxindel=tracy_config.maxindel,
        ),
    )
    json_paths = tuple(path for path in json_outputs if "HV" in path.name)
    if not json_paths:
        logger.warning("No Tracy HV JSON files were produced for sample {}", sample_id)
        result = unavailable_result(sample_id, "No Tracy HV JSON files were produced")
        return PreparedSample(sample_id, (), (result,))

    qc_results = tuple(analyze_trace_file(path, sample_id=sample_id, config=noise_config) for path in json_paths)
    return PreparedSample(sample_id, json_paths, qc_results)


def _process_prepared_sample(
    prepared: PreparedSample,
    ref_seq: str,
    *,
    noise_mask_enabled: bool,
) -> list[ProcessedTrace]:
    """Run ETL and optional per-trace noise masking for one sample."""
    qc_by_filename = {result.filename: result for result in prepared.qc_results if result.filename is not None}
    results: list[ProcessedTrace] = []
    for json_path in prepared.json_paths:
        qc_result = qc_by_filename.get(json_path.name)
        if qc_result is None:
            msg = f"No QC result found for Tracy file {json_path.name}"
            raise ValueError(msg)
        try:
            sample = process(sample_id=prepared.sample_id, input_path=json_path, ref_seq=ref_seq)
        except (ValueError, KeyError, OSError, RuntimeError) as error:
            logger.error("ETL failed for sample {} file {}: {}", prepared.sample_id, json_path.name, error)
            continue

        if noise_mask_enabled:
            mask_result = apply_noise_mask(sample, qc_result.ranges)
            sample = mask_result.sample
            excluded_variants = mask_result.excluded_variants
        else:
            excluded_variants = ()
        results.append(
            ProcessedTrace(
                sample_id=prepared.sample_id,
                filename=json_path.name,
                sample=sample,
                excluded_variants=excluded_variants,
            ),
        )
    return results


def _write_sample_qc_reports(
    preprocess_dir: Path,
    prepared_samples: list[PreparedSample],
    processed_traces: list[ProcessedTrace],
    config: NoiseConfig,
) -> None:
    """Write one deterministic QC report with exclusions per sample."""
    exclusions = {(trace.sample_id, trace.filename): trace.excluded_variants for trace in processed_traces}
    for prepared in sorted(prepared_samples, key=lambda item: item.sample_id):
        sample_dir = preprocess_dir / prepared.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        report_path = sample_dir / "qc_report.json"
        qc_results = [
            with_excluded_variants(
                result,
                exclusions.get((prepared.sample_id, result.filename or ""), ()),
            )
            for result in prepared.qc_results
        ]
        report = build_qc_report(qc_results, config)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        logger.info("Wrote Tracy QC report to {}", report_path)


def _prepare_batch(
    sample_ids: list[str],
    input_dir: Path,
    preprocess_dir: Path,
    ref_path: str,
    noise_config: NoiseConfig,
    max_workers: int,
) -> list[PreparedSample]:
    """Decompose and analyze all requested samples in parallel."""
    prepared_samples: list[PreparedSample] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _prepare_sample,
                sample_id,
                input_dir,
                preprocess_dir,
                ref_path,
                noise_config,
            ): sample_id
            for sample_id in sample_ids
        }
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                prepared = future.result()
            except (ValueError, KeyError, OSError, RuntimeError) as error:
                logger.error("Tracy preparation failed for sample {}: {}", sample_id, error)
                prepared = PreparedSample(sample_id, (), (error_result(sample_id, str(error)),))
            prepared_samples.append(prepared)
    return prepared_samples


def _process_prepared_batch(
    prepared_samples: list[PreparedSample],
    ref_seq: str,
    max_workers: int,
    *,
    noise_mask_enabled: bool,
) -> list[ProcessedTrace]:
    """Run ETL and per-trace noise masking in parallel."""
    processed: list[ProcessedTrace] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_prepared_sample,
                prepared,
                ref_seq,
                noise_mask_enabled=noise_mask_enabled,
            ): prepared.sample_id
            for prepared in prepared_samples
            if prepared.json_paths
        }
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                processed.extend(future.result())
            except (ValueError, KeyError, OSError, RuntimeError) as error:
                logger.error("ETL processing failed for sample {}: {}", sample_id, error)
    return processed


def process_batch(
    input_dir: Path,
    output_dir: Path,
    options: BatchOptions | None = None,
) -> dict[str, Sample]:
    """Process Tracy samples in parallel.

    Reads sample IDs from a TXT file, runs decompose and QC, processes each
    resulting JSON through ETL and optional noise masking, combines same-LID
    samples, and writes json/, regions/, and preprocess/ outputs.

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

    # Remove the obsolete batch-level report if this output directory is reused.
    (paths.preprocess_dir / "qc_report.json").unlink(missing_ok=True)

    # Read sample IDs from TXT file (always required)
    noise_config = NoiseConfig.from_settings()
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

    tracy_settings = get_settings().tracy
    max_workers = tracy_settings.max_workers or max(1, os.cpu_count() or 1)

    # Phase 1: decompose and analyze traces.
    prepared_samples = _prepare_batch(
        sample_ids,
        input_dir,
        paths.preprocess_dir,
        options.ref_path,
        noise_config,
        max_workers,
    )
    # Phase 2: process and mask each trace before same-sample merging.
    processed = _process_prepared_batch(
        prepared_samples,
        ref_seq,
        max_workers,
        noise_mask_enabled=noise_config.mask_enabled,
    )
    _write_sample_qc_reports(paths.preprocess_dir, prepared_samples, processed, noise_config)

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
    failed = len(set(sample_ids).difference(combined))
    if failed:
        logger.warning("{} sample(s) failed processing", failed)

    return combined


def _build_combined_samples(processed: list[ProcessedTrace]) -> dict[str, Sample]:
    """Group processed traces by sample ID and combine same-sample results."""
    grouped: dict[str, list[Sample]] = defaultdict(list)
    for trace in processed:
        grouped[trace.sample_id].append(trace.sample)

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
