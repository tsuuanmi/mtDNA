"""FASTA-input BLASTN pipeline — batch orchestration for pre-assembled FASTA.

This module is the FASTA-entry counterpart to :mod:`src.tools.blastn.pipeline`.
Where ``pipeline`` starts from AB1 trace files (AB1 -> ... -> BLASTN TSV via
:mod:`src.tools.blastn.preprocessing`), this module starts from already-assembled
FASTA files (one sample per FASTA, sample ID = filename stem) and runs a single
BLASTN alignment per sample to produce the outfmt-7 TSV the ETL consumes.

It reuses the shared ETL (:func:`src.tools.blastn.etl.process`), the batch
finalization helpers from :mod:`src.tools.blastn.pipeline`
(:func:`_build_combined_samples`, :func:`_write_full_region_jsons`),
:func:`src.core.batch.Batch.write`, and the standard
:class:`src.core.paths.OutputPaths` layout (``json/`` + ``regions/`` +
``preprocess/``). The only tool-specific step is FASTA -> BLASTN TSV; everything
downstream is identical to the AB1 pipeline, so output parity is preserved.

Output layout (same as the AB1 pipeline)::

    <output_dir>/
    |-- json/              # statistic_fullbatch.json + per-LID <LID>/<LID>.json
    |-- regions/           # HV1_{LID}.json, HV2-3_{LID}.json
    `-- preprocess/        # per-sample <LID>/<LID>.concatenate.blastn

FASTA files carry no per-base quality, so the ``.concatenate.quality`` file is
not produced; the ETL falls back to empty quality scores (consistent with the
legacy FASTA path, which also had no quality data).
"""

import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.core.batch import Batch
from src.core.models import Sample
from src.core.paths import OutputPaths
from src.tools.blastn.etl import process as etl_process
from src.tools.blastn.pipeline import _build_combined_samples, _write_full_region_jsons
from src.tools.blastn.utils import run_command

# BLASTN tabular output format used by the ETL. Must match the AB1 pipeline's
# `_blastn_multi_align` outfmt exactly so `parse_blastn_tsv` accepts the file.
_BLASTN_OUTFMT = (
    "7 qseqid sseqid sstrand nident pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore"
)

FASTA_EXTENSIONS: tuple[str, ...] = (".fasta", ".fa", ".fas", ".fna")


@dataclass
class FastaBatchOptions:
    """Configuration options for FASTA-input BLASTN batch processing.

    Tuning parameters (word_size, blastn_threads, max_workers) are loaded from
    ``get_settings().blastn`` (BlastnSettings) rather than passed through the
    CLI, matching the AB1 pipeline convention.

    Attributes:
        ref_path: Path to rCRS reference FASTA file.
        batch_id: Optional batch identifier.
    """

    ref_path: str = "ref/rCRS.fasta"
    batch_id: str | None = None


def _get_fasta_files(input_path: Path) -> list[Path]:
    """Collect FASTA files from a file or directory path.

    Args:
        input_path: A single FASTA file or a directory to search recursively.

    Returns:
        Sorted list of FASTA file paths.
    """
    if input_path.is_file():
        return [input_path]
    return sorted(f for f in input_path.rglob("*") if f.suffix.lower() in FASTA_EXTENSIONS)


def _blastn_fasta(
    sample_id: str,
    fasta_file: Path,
    preprocess_dir: Path,
    ref_path: str,
) -> str:
    """Run BLASTN on a single FASTA file and return the TSV path.

    Writes ``{preprocess_dir}/{sample_id}/{sample_id}.concatenate.blastn`` using
    outfmt 7 with the column set the ETL expects. The filename matches the AB1
    pipeline convention so :func:`src.tools.blastn.etl.process` can locate the
    (absent) quality sibling via its ``.concatenate.quality`` substitution.

    Args:
        sample_id: Sample identifier (FASTA filename stem).
        fasta_file: Path to the input FASTA file (BLASTN query).
        preprocess_dir: Base preprocess output directory.
        ref_path: Path to rCRS reference FASTA (BLASTN subject).

    Returns:
        Path to the produced BLASTN TSV file.
    """
    settings = get_settings()
    blastn_cfg = settings.blastn

    sample_dir = Path(preprocess_dir) / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    blastn_output = sample_dir / f"{sample_id}.concatenate.blastn"

    cmd = [
        str(settings.tools.blastn),
        "-query",
        str(fasta_file),
        "-subject",
        str(ref_path),
        "-word_size",
        str(blastn_cfg.word_size),
        "-outfmt",
        _BLASTN_OUTFMT,
        "-out",
        str(blastn_output),
        "-num_threads",
        str(blastn_cfg.blastn_threads),
    ]
    run_command(cmd, f"BLASTN analysis for {sample_id}")
    return str(blastn_output)


def _process_single_fasta(
    sample_id: str,
    fasta_file: Path,
    preprocess_dir: Path,
    options: FastaBatchOptions,
) -> list[dict[str, Any]]:
    """Process a single FASTA file: BLASTN alignment -> ETL -> Sample.

    Args:
        sample_id: Sample identifier (FASTA filename stem).
        fasta_file: Path to the input FASTA file.
        preprocess_dir: Directory for preprocessing output (``paths.preprocess_dir``).
        options: Batch processing options.

    Returns:
        List with a single dict ``{"sample_id", "sample"}``, or empty on failure.
    """
    try:
        blastn_tsv = _blastn_fasta(sample_id, fasta_file, preprocess_dir, options.ref_path)
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError, KeyError, OSError, RuntimeError) as e:
        logger.error("BLASTN failed for sample {}: {}", sample_id, e)
        return []

    try:
        sample = etl_process(
            sample_id=sample_id,
            input_path=blastn_tsv,
            ref_path=options.ref_path,
            batch_id=options.batch_id,
        )
    except (ValueError, KeyError, OSError, RuntimeError) as e:
        logger.error("ETL failed for sample {}: {}", sample_id, e)
        return []

    return [{"sample_id": sample_id, "sample": sample}]


def process_batch_fasta(
    input_path: Path,
    output_dir: Path,
    options: FastaBatchOptions | None = None,
) -> dict[str, Sample]:
    """Process pre-assembled FASTA files in parallel.

    Collects FASTA files (one sample per file, ID = filename stem), runs a
    BLASTN alignment per sample in parallel, runs ETL (BLASTN TSV -> Sample),
    combines same-LID Samples, and writes output organized into ``json/``,
    ``regions/``, and ``preprocess/`` subdirectories — identical to the AB1
    pipeline.

    Args:
        input_path: FASTA file or directory containing FASTA files.
        output_dir: Base output directory (subdirectories created automatically).
        options: Batch processing options.

    Returns:
        Dictionary mapping sample_id to processed Sample.
    """
    if options is None:
        options = FastaBatchOptions()

    input_path = Path(input_path)
    output_dir = Path(output_dir)

    # Construct standard output directory structure
    paths = OutputPaths(output_dir)
    paths.create_dirs()

    # Validate reference sequence is readable
    str(SeqIO.read(options.ref_path, "fasta").seq)

    fasta_files = _get_fasta_files(input_path)
    if not fasta_files:
        logger.warning("No FASTA files found in {}", input_path)
        return {}

    logger.info("Processing {} FASTA files from {}", len(fasta_files), input_path)

    max_workers = get_settings().blastn.max_workers or max(1, os.cpu_count() or 1)

    # Process each FASTA file in parallel — BLASTN output goes to preprocess_dir
    processed: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_single_fasta,
                fasta_file.stem,
                fasta_file,
                paths.preprocess_dir,
                options,
            ): fasta_file.stem
            for fasta_file in fasta_files
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

    # Group by sample_id, combine same-LID Samples, write per-region JSONs, then
    # batch JSON — mirrors the AB1 pipeline exactly. The base output dir is
    # already batch-specific, so Batch.write() uses nest_batch_id=False.
    combined = _build_combined_samples(processed)
    for sample in combined.values():
        _write_full_region_jsons(
            sample=sample,
            sample_id=sample.sample_id,
            regions_dir=paths.regions_dir,
            ref_path=options.ref_path,
        )

    batch_results = list(combined.values())
    Batch(batch_results, ref_path=options.ref_path).write(
        str(paths.json_dir),
        options.batch_id or "unknown",
        nest_batch_id=False,
    )

    logger.info(
        "FASTA BLASTN batch processing complete: {} unique LIDs from {} samples",
        len(combined),
        len(processed),
    )

    return combined


def main() -> None:
    """CLI entry point for FASTA-input BLASTN batch processing."""
    parser = argparse.ArgumentParser(
        description="BLASTN batch processing from pre-assembled FASTA files",
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="FASTA file or directory")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Base output dir, subdirs auto-created",
    )
    parser.add_argument(
        "-r",
        "--ref-path",
        type=str,
        default="ref/rCRS.fasta",
        help="Path to rCRS reference FASTA",
    )
    parser.add_argument("--batch-id", type=str, default=None, help="Batch identifier")

    args = parser.parse_args()

    options = FastaBatchOptions(ref_path=args.ref_path, batch_id=args.batch_id)

    results = process_batch_fasta(
        input_path=args.input,
        output_dir=args.output_dir,
        options=options,
    )

    logger.info("Processed {} samples", len(results))


if __name__ == "__main__":
    main()
