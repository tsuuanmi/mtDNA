"""Regenerate gate — produce the filtered "final" JSON for downstream tasks.

The regenerate step is the boundary between raw tool output and the
authoritative final JSON consumed by downstream tasks (FASTA / report
generation). Raw tool JSON (``results/tools/<tool>/json``) keeps every
called variant, including ones that fall outside the sequenced regions
(coverage gaps) — these are intentionally retained so the comparison step
can surface them to the reviewer as "additional" variants.

This module loads each raw per-sample JSON, scopes variants to
``sample.intervals`` (the sequenced regions) via the existing
:func:`filter_sample_by_regions` helper, and writes the filtered result
as ``statistic_fullbatch.json`` plus per-sample ``{LID}/{LID}.json`` into
the regenerate (final) directory. Samples with no intervals (full-coverage
FIS/NGS path) pass through unchanged.

Pipeline wiring (``scripts/pipeline.sh``):

- ``results/tools/<tool>/json`` (raw)            -> comparison
- ``regenerate/<batch>/`` (this filtered output)  -> generate_fasta / generate_reports
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path

from loguru import logger

from src.core.batch import Batch
from src.core.models import Sample
from src.core.sample import filter_sample_by_regions, load_sample

ALL_REGIONS = ["HV1", "HV2", "HV3"]


def _load_raw_samples(input_dir: Path) -> list[Sample]:
    """Load every per-sample JSON under *input_dir*.

    Expects the standard tool output layout: ``{LID}/{LID}.json`` plus a
    top-level ``statistic_fullbatch.json`` (ignored here — samples are read
    from the per-sample files, which are the source of truth).
    """
    samples: list[Sample] = []
    for sample_dir in sorted(p for p in input_dir.iterdir() if p.is_dir()):
        json_path = sample_dir / f"{sample_dir.name}.json"
        if not json_path.exists():
            logger.warning("Skipping {}: no per-sample JSON", sample_dir.name)
            continue
        try:
            samples.append(load_sample(json_path))
        except (ValueError, KeyError, OSError) as e:
            logger.error("Failed to load {}: {}", json_path, e)
    return samples


def filter_to_intervals(sample: Sample) -> Sample:
    """Return a copy of *sample* with variants scoped to its intervals.

    Variants whose position falls outside any interval in ``sample.intervals``
    are dropped. Samples with no intervals (full-coverage path) are returned
    unchanged so all called variants are preserved.
    """
    if not sample.intervals:
        return sample
    region_keys = [k for k in ALL_REGIONS if k in sample.intervals]
    if not region_keys:
        return sample
    filtered = filter_sample_by_regions(sample, region_keys)
    return filtered if filtered is not None else sample


def regenerate(input_dir: str | Path, output_dir: str | Path, ref_path: str, batch_id: str) -> int:
    """Filter raw tool JSON into the regenerate (final) directory.

    Args:
        input_dir: Raw tool JSON directory (``{LID}/{LID}.json`` layout).
        output_dir: Destination regenerate directory (written flat).
        ref_path: Reference FASTA path (for HV sequence recomputation).
        batch_id: Batch identifier (informational; output is written flat).

    Returns:
        Number of samples written.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if not input_dir.is_dir():
        logger.error("Regenerate input dir not found: {}", input_dir)
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    samples = _load_raw_samples(input_dir)
    if not samples:
        logger.warning("No samples found in {}; regenerate output empty", input_dir)
        return 0

    filtered = [filter_to_intervals(s) for s in samples]
    Batch(filtered, ref_path=ref_path).write(str(output_dir), batch_id, nest_batch_id=False)
    logger.success("Regenerate wrote {} filtered samples to {}", len(filtered), output_dir)
    return len(filtered)


def parse_args() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Regenerate filtered final JSON from raw tool output")
    parser.add_argument("-i", "--input", required=True, help="Raw tool JSON input directory")
    parser.add_argument("-o", "--output", required=True, help="Regenerate (final) output directory")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA file")
    parser.add_argument("--batch-id", default="unknown", help="Batch identifier (informational)")
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    regenerate(args.input, args.output, args.reference, args.batch_id)


if __name__ == "__main__":
    main()
