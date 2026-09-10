"""BLASTn preprocessing — AB1 → BLASTN TSV subprocess pipeline.

Runs the multi-step AB1-to-BLASTN pipeline: AB1 → FASTQ (UGENE) →
trimmed FASTQ (seqtk) → FASTA (seqtk) → concatenated BLASTN alignment.

Called by pipeline.py only; never by etl.py (per TOOLS_STANDARD.md §4.2).

This module is idempotent: running preprocess twice on the same input
with the same config produces the same output (assuming deterministic
external tools).
"""

import shutil
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.tools.blastn.utils import run_command


class PreprocessResult(NamedTuple):
    """Result of preprocessing a single sample.

    Attributes:
        blastn_tsv_path: Path to the concatenated BLASTN alignment TSV file.
            This is the primary input for the ETL module.
        quality_path: Path to the concatenated quality scores file, or None
            if quality extraction failed.
        concatenated_fasta_path: Path to the concatenated FASTA file, or None
            if concatenation failed.
    """

    blastn_tsv_path: Path
    quality_path: Path | None = None
    concatenated_fasta_path: Path | None = None


def _ab1_to_fastq(ab1_path: str | Path, fastq_path: str | Path, ugene_path: str) -> None:
    """Convert AB1 file to FASTQ format using UGENE.

    Args:
        ab1_path: Path to input AB1 file.
        fastq_path: Path to output FASTQ file.
        ugene_path: Path to UGENE executable.
    """
    cmd = [
        ugene_path,
        "convert-seq",
        f"--in={ab1_path}",
        f"--out={fastq_path}",
        "--format=fastq",
    ]
    run_command(cmd, f"AB1 to FASTQ conversion for {Path(ab1_path).name}")


def _fastq_trim(fastq_path: str | Path, out_path: str | Path, threshold: float, seqtk_path: str) -> None:
    """Trim FASTQ file based on quality threshold using seqtk.

    Args:
        fastq_path: Path to input FASTQ file.
        out_path: Path to output trimmed FASTQ file.
        threshold: Quality threshold for trimming.
        seqtk_path: Path to seqtk executable.
    """
    cmd = [
        seqtk_path,
        "trimfq",
        "-q",
        str(threshold),
        str(fastq_path),
    ]
    result = run_command(cmd, f"seqtk trimfq for {Path(fastq_path).name}")
    Path(out_path).write_text(result)


def _fastq_to_fasta(fastq_path: str | Path, out_path: str | Path, seqtk_path: str) -> None:
    """Convert FASTQ to FASTA format using seqtk.

    Args:
        fastq_path: Path to input FASTQ file.
        out_path: Path to output FASTA file.
        seqtk_path: Path to seqtk executable.
    """
    cmd = [
        seqtk_path,
        "seq",
        "-a",
        str(fastq_path),
    ]
    result = run_command(cmd, f"seqtk seq for {Path(fastq_path).name}")
    Path(out_path).write_text(result)


def _extract_phred_quality(fastq_path: str | Path, quality_out_path: str | Path) -> None:
    """Extract PHRED quality scores from a FASTQ file.

    Args:
        fastq_path: Path to input FASTQ file.
        quality_out_path: Path to output quality file.
    """
    fastq_path = Path(fastq_path)
    quality_out_path = Path(quality_out_path)

    with quality_out_path.open("w") as quality_file:
        for record in SeqIO.parse(str(fastq_path), "fastq"):
            quality_scores = record.letter_annotations["phred_quality"]
            quality_file.write(f">{record.id}\n")
            quality_file.write(",".join(map(str, quality_scores)) + "\n")


def _blastn_multi_align(
    query_path: str | Path,
    subject_path: str | Path,
    out_path: str | Path,
    word_size: int,
    num_threads: int,
    blastn_path: str,
) -> None:
    """Run BLASTN multiple sequence alignment with tabular output.

    Uses format 7 (with comment lines) for parsing compatibility with the
    ETL module. This produces the TSV file that the ETL reads.

    Args:
        query_path: Path to concatenated query FASTA file.
        subject_path: Path to subject (reference) FASTA file.
        out_path: Path to output BLASTN TSV file.
        params: BLASTN command parameters.
    """
    cmd = [
        blastn_path,
        "-query",
        str(query_path),
        "-subject",
        str(subject_path),
        "-word_size",
        str(word_size),
        "-outfmt",
        "7 qseqid sseqid sstrand nident pident length mismatch "
        "gapopen qstart qend sstart send qseq sseq evalue bitscore",
        "-out",
        str(out_path),
        "-num_threads",
        str(num_threads),
    ]
    run_command(cmd, f"BLASTN multi alignment for {Path(query_path).name}")


def _concatenate_fasta(output_dir: Path, out_path: str | Path) -> None:
    """Concatenate all trimmed FASTA files in output directory.

    Args:
        output_dir: Directory containing trimmed FASTA files.
        out_path: Path to output concatenated FASTA file.
    """
    with Path(out_path).open("w") as outfile:
        for fasta_file in sorted(output_dir.glob("*.trim.fasta")):
            with fasta_file.open() as infile:
                content = infile.read()
                if not content.endswith("\n"):
                    content += "\n"
                outfile.write(content)


def _concatenate_quality(output_dir: Path, out_path: str | Path) -> None:
    """Concatenate all quality files in output directory.

    Args:
        output_dir: Directory containing quality files.
        out_path: Path to output concatenated quality file.
    """
    with Path(out_path).open("w") as outfile:
        for quality_file in sorted(output_dir.glob("*.trim.quality")):
            with quality_file.open() as infile:
                outfile.write(infile.read())


def preprocess_sample(
    sample_id: str,
    list_ab1: list[str],
    output_dir: str | Path,
    ref_path: str | Path,
) -> PreprocessResult:
    """Preprocess a single sample: AB1 files → BLASTN TSV alignment.

    Runs the full multi-step pipeline:
    1. AB1 → FASTQ (UGENE convert-seq)
    2. FASTQ → trimmed FASTQ (seqtk trimfq)
    3. Trimmed FASTQ → FASTA (seqtk seq)
    4. Extract Phred quality scores
    5. Concatenate FASTA files
    6. Concatenate quality files
    7. BLASTN multi-sequence alignment (outfmt 7, primary ETL input)

    Each subprocess failure is logged and raises an exception. The
    caller (pipeline) should catch and skip failed samples.

    Args:
        sample_id: Sample identifier (LID).
        list_ab1: List of AB1 file paths for this sample.
        output_dir: Output directory for this sample.
        ref_path: Path to reference FASTA file (rCRS).

    Returns:
        PreprocessResult with the concatenated BLASTN TSV path as the
        primary ETL input, plus optional quality and FASTA paths.

    Raises:
        subprocess.CalledProcessError: If any subprocess fails.
        FileNotFoundError: If tool binaries are not found.
    """
    settings = get_settings()
    blastn_cfg = settings.blastn
    # Keep BLASTN single-threaded: samples already run in parallel via
    # ProcessPoolExecutor, so per-subprocess threads would oversubscribe the CPU.
    num_threads = blastn_cfg.blastn_threads

    output_dir = Path(output_dir) / sample_id
    output_dir.mkdir(parents=True, exist_ok=True)
    ref_path = Path(ref_path)

    # Validate tool binaries
    tools = settings.tools
    for tool_name, tool_attr in [("UGENE", "ugene"), ("BLASTN", "blastn"), ("seqtk", "seqtk")]:
        tool_path = getattr(tools, tool_attr)
        if shutil.which(str(tool_path)) is None:
            msg = f"{tool_name} binary not found: {tool_path}"
            raise FileNotFoundError(msg)

    # Step 1: Process individual AB1 files
    for ab1_file in list_ab1:
        name = Path(ab1_file).stem
        fastq_path = output_dir / f"{name}.fastq"
        trim_fastq_path = output_dir / f"{name}.trim.fastq"
        quality_path = output_dir / f"{name}.trim.quality"
        trim_fasta_path = output_dir / f"{name}.trim.fasta"

        # AB1 → FASTQ
        _ab1_to_fastq(ab1_file, str(fastq_path), str(tools.ugene))

        # FASTQ → trimmed FASTQ
        _fastq_trim(str(fastq_path), str(trim_fastq_path), blastn_cfg.trim_threshold, str(tools.seqtk))

        # Extract quality scores
        _extract_phred_quality(str(trim_fastq_path), str(quality_path))

        # FASTQ → FASTA
        _fastq_to_fasta(str(trim_fastq_path), str(trim_fasta_path), str(tools.seqtk))

    # Step 2: Concatenate FASTA and quality files
    concat_fasta_path = output_dir / f"{sample_id}.concatenate.fasta"
    concat_quality_path = output_dir / f"{sample_id}.concatenate.quality"

    _concatenate_fasta(output_dir, str(concat_fasta_path))
    _concatenate_quality(output_dir, str(concat_quality_path))

    # Step 3: BLASTN multi-sequence alignment (outfmt 7 — primary ETL input)
    concat_blastn_path = output_dir / f"{sample_id}.concatenate.blastn"
    _blastn_multi_align(
        str(concat_fasta_path),
        str(ref_path),
        str(concat_blastn_path),
        word_size=blastn_cfg.word_size,
        num_threads=num_threads,
        blastn_path=str(tools.blastn),
    )

    logger.info("Preprocessing complete for sample {}: {}", sample_id, concat_blastn_path)

    return PreprocessResult(
        blastn_tsv_path=concat_blastn_path,
        quality_path=concat_quality_path if concat_quality_path.exists() else None,
        concatenated_fasta_path=concat_fasta_path if concat_fasta_path.exists() else None,
    )
