"""BLASTn pipeline module — batch orchestration and parallel processing.

Handles batch processing of BLASTN samples with parallel execution,
following the standard tool module pattern (ETO → pipeline → Batch.write).

The actual ETL (BLASTN TSV → Sample) lives in ``src.tools.blastn.etl``.
The preprocessing (AB1 → BLASTN TSV) lives in ``src.tools.blastn.preprocessing``.
This module owns batch orchestration and parallelism only.

After batch processing, output is organized into three subdirectories
under the base output directory:

- ``json/`` — ``Batch.write()`` output written flat under this dir
  (``statistic_fullbatch.json`` and per-sample ``<LID>/<LID>.json``)
- ``regions/`` — per-region intermediate JSONs (``HV1_{LID}.json``,
  ``HV2-3_{LID}.json``) written by ``write_region_jsons()``
- ``preprocess/`` — tool-specific preprocessing artifacts (concatenated
  ``.blastn`` plus ``.fasta``, ``.fastq``, ``.quality`` files)

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
from src.core.variants import pos_sort_key
from src.tools.blastn.etl import process
from src.tools.blastn.preprocessing import preprocess_sample
from src.tools.blastn.utils import get_files_by_id


@dataclass
class BatchOptions:
    """Configuration options for BLASTN batch processing.

    Tuning parameters (trim_threshold, word_size, max_workers) are loaded
    from ``get_settings().blastn`` (BlastnSettings) rather than passed
    through the CLI or this dataclass.

    Attributes:
        ref_path: Path to rCRS reference FASTA file.
        batch_id: Optional batch identifier.
        samples_path: Path to TXT file with sample IDs (one per line).
    """

    ref_path: str = "ref/rCRS.fasta"
    batch_id: str | None = None
    samples_path: Path | None = None


def _read_sample_ids(samples_path: Path) -> list[str]:
    """Read sample IDs from a TXT file (one per line).

    Args:
        samples_path: Path to TXT file with sample IDs.

    Returns:
        List of sample ID strings.
    """
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

    Concatenates variants, merges intervals, and unions flags from region
    Samples of the same tool (BLASTN). This is NOT the cross-tool merge
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
        source_tool=Tool.BLASTN,
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
    """Write per-region JSONs for a combined (full-region) BLASTN sample.

    Mirrors Tracy/Sequencher's ``_write_full_region_jsons``: all region groups
    are emitted and ``write_region_jsons`` skips groups without data. BLASTN
    does not set auto-pass flags, so region groups are always derived from
    ``REGION_TO_KEYS`` (data-driven) rather than from ``sample_flags``.

    BLASTN combines forward/reverse primer traces per LID before this is
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
    files: list[str],
    preprocess_dir: Path,
    options: BatchOptions,
) -> list[dict[str, Any]]:
    """Process a single BLASTN sample: preprocess → ETL.

    Args:
        sample_id: Sample identifier (LID).
        files: List of AB1 file paths for this sample.
        preprocess_dir: Directory for preprocessing output (``paths.preprocess_dir``).
        options: Batch processing options.

    Returns:
        List of dicts with sample_id and sample fields (one per region).
    """
    try:
        preprocess_result = preprocess_sample(
            sample_id=sample_id,
            list_ab1=files,
            output_dir=str(preprocess_dir),
            ref_path=options.ref_path,
        )
    except (ValueError, KeyError, OSError, RuntimeError) as e:
        logger.error("Preprocessing failed for sample {}: {}", sample_id, e)
        return []

    try:
        sample = process(
            sample_id=sample_id,
            input_path=preprocess_result.blastn_tsv_path,
            ref_path=options.ref_path,
            batch_id=options.batch_id,
        )
    except (ValueError, KeyError, OSError, RuntimeError) as e:
        logger.error("ETL failed for sample {}: {}", sample_id, e)
        return []

    return [{"sample_id": sample_id, "sample": sample}]


def process_batch(
    input_dir: Path,
    output_dir: Path,
    options: BatchOptions | None = None,
) -> dict[str, Sample]:
    """Process BLASTN samples in parallel.

    Reads sample IDs from a TXT file, preprocesses each sample (AB1→BLASTN),
    runs ETL (BLASTN TSV→Sample) in parallel, combines same-LID Samples,
    and writes output organized into json/, regions/, preprocess/ subdirectories.

    Args:
        input_dir: Directory containing AB1 files.
        output_dir: Base output directory (subdirectories created automatically).
        options: Batch processing options.

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

    # Load reference sequence
    str(SeqIO.read(options.ref_path, "fasta").seq)

    # Read sample IDs from TXT file
    if options.samples_path is None:
        logger.error("No samples_path provided — sample IDs must come from a TXT file")
        return {}

    sample_ids = _read_sample_ids(options.samples_path)
    if not sample_ids:
        logger.warning("No sample IDs found in {}", options.samples_path)
        return {}

    logger.info("Loaded {} sample IDs from {}", len(sample_ids), options.samples_path)

    # Get AB1 files for each sample
    mapping_id_file = get_files_by_id(str(options.samples_path), str(input_dir))

    max_workers = get_settings().blastn.max_workers or max(1, os.cpu_count() or 1)

    # Process each sample in parallel — preprocessing output goes to preprocess_dir
    processed: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_single_sample,
                sample_id,
                files,
                paths.preprocess_dir,
                options,
            ): sample_id
            for sample_id, files in mapping_id_file.items()
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
    # The base output dir is already batch-specific (e.g. results/tools/blastn/<batch_id>),
    # so Batch.write() writes flat into json_dir (nest_batch_id=False) to avoid a
    # duplicated <batch_id> path segment. Output lands directly under json_dir:
    # statistic_fullbatch.json plus per-LID <LID>/<LID>.json. (Matches the Tracy
    # reference pipeline.)
    combined = _build_combined_samples(processed)
    for sample in combined.values():
        _write_full_region_jsons(
            sample=sample,
            sample_id=sample.sample_id,
            regions_dir=paths.regions_dir,
            ref_path=options.ref_path,
        )

    # Write batch JSON to json_dir.
    batch_results = list(combined.values())
    Batch(batch_results, ref_path=options.ref_path).write(
        str(paths.json_dir),
        options.batch_id or "unknown",
        nest_batch_id=False,
    )

    logger.info(
        "BLASTN batch processing complete: {} unique LIDs from {} samples",
        len(combined),
        len(processed),
    )

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
    """CLI entry point for BLASTN batch processing."""

    parser = argparse.ArgumentParser(description="BLASTN batch processing")
    parser.add_argument("-i", "--input-dir", type=Path, required=True, help="Directory containing AB1 files")
    parser.add_argument("-o", "--output-dir", type=Path, required=True, help="Base output dir, subdirs auto-created")
    parser.add_argument("-s", "--samples", type=Path, required=True, help="Path to TXT file with sample IDs")
    parser.add_argument("--batch-id", type=str, default=None, help="Batch identifier")
    parser.add_argument("--ref-path", type=str, default="ref/rCRS.fasta", help="Path to rCRS reference FASTA")

    args = parser.parse_args()

    options = BatchOptions(
        ref_path=args.ref_path,
        batch_id=args.batch_id,
        samples_path=args.samples,
    )

    results = process_batch(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        options=options,
    )

    logger.info("Processed {} samples", len(results))


if __name__ == "__main__":
    main()
