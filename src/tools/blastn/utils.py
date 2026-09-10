"""BLASTn module utilities — TSV parsing, strand handling, quality extraction.

Pure utility functions with no Sample/Variant imports. BLASTN TSV parsing uses
pandas for format-7 comment-line handling, then converts to plain data structures
(lists of dicts) immediately. ETL never touches DataFrames directly.
"""

import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.Seq import Seq
from loguru import logger

# ---------------------------------------------------------------------------
# BLASTN TSV format constants
# ---------------------------------------------------------------------------

BLASTN_COLUMNS: list[str] = [
    "qseqid",
    "sseqid",
    "sstrand",
    "nident",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qseq",
    "sseq",
    "evalue",
    "bitscore",
]
"""Column names for BLASTN tabular output format 7 (with comment lines)."""

BLASTN_COMMENT_PREFIX: str = "#"
"""Comment line prefix for BLASTN format-7 output."""


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------


def parse_blastn_tsv(path: str | Path) -> list[dict[str, Any]]:
    """Parse a BLASTN format-7 TSV file into a list of row dicts.

    Uses pandas for initial parsing (handles comment lines and column
    assignment), then converts to plain data structures immediately.
    ETL works with list-of-dicts, not DataFrames.

    Args:
        path: Path to BLASTN TSV output file.

    Returns:
        List of dicts, one per alignment row. Each dict has keys from
        BLASTN_COLUMNS plus forward-facing fields added by
        ``normalize_strand()``.

    Raises:
        FileNotFoundError: If path does not exist.
        ValueError: If TSV has wrong number of columns.
    """
    path = Path(path)
    if not path.exists():
        msg = f"BLASTN TSV file not found: {path}"
        raise FileNotFoundError(msg)

    df = pd.read_csv(path, sep="\t", header=None, comment=BLASTN_COMMENT_PREFIX)

    if len(df.columns) != len(BLASTN_COLUMNS):
        msg = f"BLASTN TSV has {len(df.columns)} columns, expected {len(BLASTN_COLUMNS)}. Path: {path}"
        raise ValueError(msg)

    df.columns = BLASTN_COLUMNS

    # Convert to list of dicts immediately — ETL uses plain data structures
    return df.to_dict("records")


# ---------------------------------------------------------------------------
# Strand direction handling
# ---------------------------------------------------------------------------


def normalize_strand(row: dict[str, Any]) -> dict[str, Any]:
    """Compute forward-facing alignment fields from a BLASTN row.

    For plus-strand alignments, forward fields are identical to the
    original fields. For minus-strand alignments, reference and sample
    sequences are reverse-complemented, and start/end coordinates are
    swapped.

    This corresponds to ``Analysis._process_strand_directions()``
    applied row-by-row instead of via DataFrame operations.

    Args:
        row: A single BLASTN alignment row dict (from ``parse_blastn_tsv``).

    Returns:
        The input dict with four additional fields:
        ``forward_ref_seq``, ``forward_sample_seq``,
        ``forward_ref_start``, ``forward_ref_end``.
    """
    sstrand = row.get("sstrand", "plus")

    if sstrand == "plus":
        row["forward_ref_seq"] = row.get("sseq", "")
        row["forward_sample_seq"] = row.get("qseq", "")
        row["forward_ref_start"] = row.get("sstart", 0)
        row["forward_ref_end"] = row.get("send", 0)
    else:
        # Reverse-complement reference and sample sequences
        sseq = row.get("sseq", "")
        qseq = row.get("qseq", "")
        row["forward_ref_seq"] = str(Seq(sseq).reverse_complement()) if sseq else ""
        row["forward_sample_seq"] = str(Seq(qseq).reverse_complement()) if qseq else ""
        # Swap start/end for reverse strand
        row["forward_ref_start"] = row.get("send", 0)
        row["forward_ref_end"] = row.get("sstart", 0)

    return row


def normalize_strand_all(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Apply ``normalize_strand()`` to all rows in a BLASTN result set.

    Args:
        rows: List of BLASTN alignment row dicts.

    Returns:
        The same list with each row augmented with forward-facing fields.
    """
    for row in rows:
        normalize_strand(row)
    return rows


# ---------------------------------------------------------------------------
# Quality score extraction
# ---------------------------------------------------------------------------


def extract_quality_scores(quality_file_path: str | Path) -> dict[str, list[int]]:
    """Extract Phred quality scores from a quality file.

    Each line in the quality file starts with ``>sample_id`` followed by
    comma-separated Phred scores on the next line.

    Args:
        quality_file_path: Path to the quality file.

    Returns:
        Dict mapping sample IDs to lists of Phred quality scores.
    """
    quality_file_path = Path(quality_file_path)
    quality_scores: dict[str, list[int]] = {}

    if not quality_file_path.exists():
        logger.warning("Quality file not found: {}", quality_file_path)
        return quality_scores

    current_id: str | None = None
    with quality_file_path.open() as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith(">"):
                current_id = line[1:].strip()
                quality_scores[current_id] = []
            elif current_id is not None and line:
                try:
                    scores = [int(s) for s in line.split(",") if s.strip()]
                    quality_scores[current_id].extend(scores)
                except ValueError:
                    logger.warning("Invalid quality score in line: {}", line[:50])

    return quality_scores


def run_command(cmd: list[str], description: str) -> str:
    """Run a command and handle errors.

    Args:
        cmd: Command and arguments to run
        description: Description of the command for error logging

    Returns:
        Command stdout if successful.

    Raises:
        subprocess.CalledProcessError: If command fails
        Exception: For other errors during execution
    """
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"{description} failed: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error during {description}: {e}")
        raise
    return result.stdout


def get_files_by_id(list_id: str, ab1_dir: str) -> dict[str, list[str]]:
    """Get files by ID.

    Args:
        list_id: Path to list of IDs
        ab1_dir: AB1 directory

    Returns:
        Dictionary mapping IDs to file paths
    """
    lids: list[str] = []
    with Path(list_id).open() as f:
        lids.extend(line.replace("\n", "").split("\t")[0] for line in f)

    mapping_dict: dict[str, list[str]] = defaultdict(list)
    ab1_path = Path(ab1_dir)

    for lid in lids:
        for file in ab1_path.iterdir():
            if lid in file.name:
                mapping_dict[lid].append(str(file))

        if lid not in mapping_dict:
            logger.error(f"Cannot find any ab1 file for id: {lid}")

    return mapping_dict
