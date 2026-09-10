"""Core standardized data models for the mtDNA pipeline.

All tools (automate, sequencher, mutation_surveyor, tracy, blastn) must convert
their output to these models so the comparison and flagging layers are unified.
"""

from enum import StrEnum
from typing import Any

from pydantic import BaseModel, Field

# ---------------------------------------------------------------------------
# Position type — canonical source of truth
# ---------------------------------------------------------------------------

Position = int | str
"""
Canonical genomic position type for the mtDNA pipeline.

- Base positions are ``int``  (e.g. 73, 309, 16569).
- Insertion positions are ``str`` (e.g. "309.1", "217.10").

Every module must import ``Position`` from ``src.core.models`` rather than
re-defining the alias locally.  Use ``validate_position()`` at external
boundaries (Excel reads, API inputs) to normalise and enforce the type.
"""


def validate_position(pos: float | str) -> Position:
    """Normalise and validate a position value into canonical Position type.

    Enforces the Position contract at external boundaries (Excel reads, API
    inputs, etc.).  Rejects ``float`` — callers must convert explicitly
    (``int(x)`` for base positions, or ``str`` for insertion positions) so
    that precision loss is a conscious decision.

    Args:
        pos: Raw position value.  Accepted types:
            - ``int`` → returned as-is.
            - ``str`` like ``"309"`` → converted to ``int(309)``.
            - ``str`` like ``"309.1"`` → returned as ``"309.1"``.
            - ``float`` → **rejected** (raises TypeError).

    Returns:
        Canonical Position (int for bases, str for insertions).

    Raises:
        TypeError: If *pos* is a float (precision-loss risk) or an unsupported type.

    Examples:
        >>> validate_position(73)
        73
        >>> validate_position("309")
        309
        >>> validate_position("309.1")
        '309.1'
        >>> validate_position(309.1)
        Traceback (most recent call last):
            ...
        TypeError: ...

    """
    if isinstance(pos, int):
        return pos
    if isinstance(pos, float):
        msg = (
            f"Float position {pos!r} is ambiguous — "
            f"use int({int(pos)}) for base positions or "
            f'str("{pos}") for insertions to avoid precision loss.'
        )
        raise TypeError(
            msg,
        )
    if isinstance(pos, str):
        stripped = pos.strip().replace(",", "")
        # Try integer (base position)
        try:
            return int(stripped)
        except ValueError:
            pass
        # Has decimal → insertion position
        if "." in stripped:
            # Verify it's numeric
            try:
                float(stripped)
            except ValueError:
                msg = f"Cannot parse position string: {pos!r}"
                raise TypeError(msg) from None
            return stripped
        # Non-numeric string
        msg = f"Cannot parse position string: {pos!r}"
        raise TypeError(msg)
    msg = f"Unsupported position type {type(pos).__name__}: {pos!r}"
    raise TypeError(msg)


class Tool(StrEnum):
    """Enumeration of supported analysis tools.

    ``UNKNOWN`` is used as a sentinel when a batch/sample carries no tool
    attribution; it is not a real analysis tool and should display as
    "Unknown" in comparison output.
    """

    SEQUENCHER = "sequencher"
    MUTATION_SURVEYOR = "mutation_surveyor"
    TRACY = "tracy"
    BLASTN = "blastn"
    FASTA = "fasta"
    FIS = "fis"
    SANGER = "sanger"
    UNIFIED = "unified"
    UNKNOWN = "unknown"


class Variant(BaseModel):
    """Single variant entry — the atomic unit of comparison.

    All fields use snake_case for consistency with Pydantic conventions.
    pos must be pre-normalized (int for bases, str for insertions) before construction.
    """

    pos: int | str = Field(
        ...,
        description='1-based genomic pos (int for base positions; str for insertions, e.g. "309.1", "217.10")',
    )
    ref: str = Field(..., description="Reference base ( '-' for insertions)")
    seq: str = Field(..., description="Alternate sequence base ( '-' for deletions)")
    quality: list[int] = Field(
        default_factory=list, description="Per-trace Phred-like quality scores (one per source file)."
    )
    files: list[str] = Field(
        default_factory=list,
        description="Source trace filenames that contributed to this variant",
    )
    peaks: list[list[int | None]] | None = Field(
        default=None, description="Per-variant peak data: list of [A, C, G, T] peak lists, one per source file."
    )


class Sample(BaseModel):
    """Standardized result that every tool must produce.

    This is the single interchange format between ETL layers, comparison, and flagging.
    """

    sample_id: str = Field(..., description="Sample identifier (LID from lab)")
    variants: list[Variant] = Field(default_factory=list, description="List of detected variants")
    source_tool: Tool = Field(..., description="Primary analysis tool")
    sample_flags: list[str] = Field(
        default_factory=list,
        description="Sample-level flag reasons (e.g. range problems, QC issues).",
    )
    hv1: str | None = Field(default=None, description="HV1 consensus sequence (16024-16365)")
    hv2: str | None = Field(default=None, description="HV2 consensus sequence (73-340)")
    hv3: str | None = Field(default=None, description="HV3 consensus sequence (438-576)")
    intervals: dict[str, list[list[int]]] | None = Field(
        default=None,
        description="Per-HV-region intervals, e.g. {'HV1': [[16024, 16365]], 'HV2': [[73, 340]], 'HV3': [[438, 576]]}",
    )
    no_snps: int | None = Field(default=None, description="Number of SNP variants")
    no_ins: int | None = Field(default=None, description="Number of insertion variants")
    no_dels: int | None = Field(default=None, description="Number of deletion variants")
    no_identicals: int | None = Field(default=None, description="Number of positions identical to reference")
    no_unread: int | None = Field(default=None, description="Number of unread/N positions in HV regions")
    batch_id: str | None = Field(default=None, description="Batch identifier for the analysis run")
    variant_flags: dict[str, list[str]] = Field(
        default_factory=dict,
        description='Per-variant flag reasons keyed by "pos|ref|alt". '
        "Only includes variants that have at least one flag. "
        "Computed at serialization time via flag_variants().",
    )
    information: dict[str, Any] | None = Field(
        default=None,
        description="Sample metadata dict. Always includes source_tool and batch_id at "
        "serialization time. Tool-specific fields (e.g. type_sample, place_sample, "
        "test_time, test_date, well) are populated by ETL modules. "
        "Use batch_id (not batch) as the canonical key for the batch identifier.",
    )


class ToolResult(BaseModel):
    """Per-tool data extracted from ComparisonResult A/B pairs."""

    tool: Tool = Field(..., description="Source tool")
    analyzed_range: str = Field("", description="Analyzed range string")
    variants: list[dict[str, Any]] = Field(
        default_factory=list,
        description="Variant dicts in canonical {pos, ref, seq} format",
    )
    flagged: bool = Field(default=False, description="Whether this tool's result is flagged")
    flag_reasons: list[str] = Field(default_factory=list, description="Flag reasons for this tool")
    variant_flags: dict[str, list[str]] = Field(
        default_factory=dict,
        description='Per-variant flag reasons keyed by "pos|ref|alt". '
        "Only includes variants that have at least one flag.",
    )


class ComparisonResult(BaseModel):
    """Result of comparing two Sample objects."""

    sample_id: str = Field(..., description="Sample identifier")
    tool_a: ToolResult = Field(..., description="Tool A result data")
    tool_b: ToolResult = Field(..., description="Tool B result data")
    unique_variants_a: list[dict[str, Any]] = Field(
        default_factory=list,
        description="Variants unique to tool A (comparison-relative)",
    )
    unique_variants_b: list[dict[str, Any]] = Field(
        default_factory=list,
        description="Variants unique to tool B (comparison-relative)",
    )
    concordance: str = Field("N", description="Y/N concordance between the two results")


class Sequence(BaseModel):
    """Structured result from generate_sequence().

    Replaces the old 4/6-element tuple return, eliminating unpack-with-discard
    footguns and making each field self-documenting.
    """

    reference_seq: str = Field("", description="Full reference sequence string")
    consensus_seq: str = Field("", description="Full consensus sequence string")
    hv_seq_refs: dict[str, Any] = Field(
        default_factory=dict,
        description="Per-region reference sequences, e.g. {'HV1': {'start': ..., 'seq': [...]}}",
    )
    hv_seqs: dict[str, Any] = Field(
        default_factory=dict,
        description="Per-region consensus sequences, e.g. {'HV1': {'start': ..., 'seq': [...]}}",
    )
    ref_dict: dict[int | str, str] = Field(
        default_factory=dict,
        description="Position-keyed reference dict (int for bases, str for insertions like '309.1')",
    )
    con_dict: dict[int | str, str] = Field(
        default_factory=dict,
        description="Position-keyed consensus dict (int for bases, str for insertions like '309.1')",
    )
