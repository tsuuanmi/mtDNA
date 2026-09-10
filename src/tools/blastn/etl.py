"""BLASTn ETL module — BLASTN TSV → Sample.

Pure ETL: reads a single BLASTN alignment TSV file and produces a
Sample object with variants, intervals, and flags. No file I/O beyond
reading the input TSV. No subprocess calls. No side effects.

The entry point is ``process()`` which takes a BLASTN TSV file path and
returns a Sample. For batch processing, use
``src.tools.blastn.pipeline.process_batch``.
"""

import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Sample, Tool, Variant
from src.core.sample import merge_intervals
from src.core.variants import (
    get_hv_region_for_position,
    get_overlap,
    normalize_position,
    nucleotide_to_iupac,
    pos_base,
    pos_sort_key,
    split_variant_types,
)
from src.tools.blastn.transform import (
    _apply_conversion_rules,
    _apply_region_validation_rules,
    _handle_polyc_removal,
    _polyc_warning_flags,
    _right_norm_indel,
)
from src.tools.blastn.utils import (
    extract_quality_scores,
    normalize_strand_all,
    parse_blastn_tsv,
)

__all__ = ["process"]

# Minimum number of distinct strands for a variant to be considered
# homozygous (observed on both strands) during polyC consensus filtering.
_HOMOZYGOUS_MIN_STRANDS: int = 2

# ---------------------------------------------------------------------------
# Helper: load filter positions
# ---------------------------------------------------------------------------


def _load_filter_positions() -> set[str]:
    """Load positions to filter from rules/filter/filter.txt.

    Returns:
        Set of position strings to filter out entirely.
    """
    filter_positions: set[str] = set()
    filter_file = Path(str(get_settings().directories.filter_rules)) / "filter.txt"

    if filter_file.exists():
        with filter_file.open() as f:
            for raw_line in f:
                line = raw_line.strip()
                if line and not line.startswith("#") and not line.startswith("//"):
                    filter_positions.add(line)
    else:
        default_positions = ["16193", "455", "463", "573", "309", "309.1", "309.2"]
        logger.warning("Filter positions file not found at {}, using defaults", filter_file)
        filter_positions = set(default_positions)

    return filter_positions


def _load_no_filter_positions() -> set[str]:
    """Load positions that should never be filtered from rules/filter/no_filter.txt.

    Returns:
        Set of position strings that should not be filtered.
    """
    no_filter_positions: set[str] = set()
    no_filter_file = Path(str(get_settings().directories.filter_rules)) / "no_filter.txt"

    if no_filter_file.exists():
        with no_filter_file.open() as f:
            for raw_line in f:
                line = raw_line.strip()
                if line and not line.startswith("#") and not line.startswith("//"):
                    no_filter_positions.add(line)
    else:
        logger.warning("No-filter positions file not found at {}", no_filter_file)

    return no_filter_positions


# ---------------------------------------------------------------------------
# PolyC region helpers
# ---------------------------------------------------------------------------


def _is_in_polyc_region(pos: float | str) -> tuple[bool, str | None]:
    """Check if a position is within a polycytosine region."""
    pos_int = pos_base(pos)
    for region_name, region_range in get_settings().regions.POLYC_REGIONS.items():
        if region_range[0] <= pos_int <= region_range[1]:
            return True, region_name
    return False, None


def _is_after_polyc_region(pos: float | str) -> tuple[bool, str | None]:
    """Check if a position is after a polycytosine region (within 20 bases)."""
    pos_int = pos_base(pos)
    for region_name, region_range in get_settings().regions.POLYC_REGIONS.items():
        if pos_int > region_range[1] and pos_int <= region_range[1] + 20:
            return True, region_name
    return False, None


def _is_before_polyc_region(pos: float | str) -> tuple[bool, str | None]:
    """Check if a position is before a polycytosine region (within 20 bases)."""
    pos_int = pos_base(pos)
    for region_name, region_range in get_settings().regions.POLYC_REGIONS.items():
        if pos_int < region_range[0] and pos_int >= region_range[0] - 20:
            return True, region_name
    return False, None


def _is_polyc_stretch(region: str | None, ref_seq: str, context_seq: str | None = None) -> bool:
    """Check if there is a polycytosine stretch in the specified region."""
    if not ref_seq and not context_seq:
        logger.warning("No sequence available to check for polyC stretch")
        return False

    if not context_seq and ref_seq:
        if region is None:
            return False
        region_range = get_settings().regions.POLYC_REGIONS.get(region)
        if not region_range:
            return False
        start_idx = max(0, region_range[0] - 5 - 1)  # -1 for 0-based indexing
        end_idx = min(len(ref_seq), region_range[1] + 5)
        context_seq = ref_seq[start_idx:end_idx]

    if not context_seq:
        return False

    return bool(re.search(r"C{6,}", context_seq))


def _analyze_polyc_context(
    pos: float | str,
    ref_seq: str,
    sample_seq: str,
    local_index: int,
    strand: str = "plus",
) -> dict[str, Any]:
    """Analyze the context around a potential variant in a polycytosine region.

    Returns a dict of polyC context metadata including ``flagged_for_filtering``
    (True when a plus-strand variant is after a polyC region, or a minus-strand
    variant is before one).
    """
    result: dict[str, Any] = {}

    in_polyc, region_name = _is_in_polyc_region(pos)
    after_polyc, after_region = _is_after_polyc_region(pos)
    before_polyc, before_region = _is_before_polyc_region(pos)

    if not any([in_polyc, after_polyc, before_polyc]):
        return result

    context_start = max(0, local_index - 5)
    context_end = min(len(ref_seq), local_index + 6)

    ref_context = ref_seq[context_start:context_end]
    sample_context = sample_seq[context_start:context_end]

    should_filter = False
    if (strand == "plus" and after_polyc) or (strand == "minus" and before_polyc):
        should_filter = True

    has_polyc_stretch = False
    if in_polyc:
        region = region_name or after_region or before_region
        has_polyc_stretch = _is_polyc_stretch(region, ref_seq, ref_context)

        # Special case: 16189T>C creates polyC stretch
        if (
            region == "HV1"
            and pos_base(pos) == 16189  # noqa: PLR2004
            and ref_seq[local_index] == "T"
            and sample_seq[local_index] == "C"
        ):
            has_polyc_stretch = True

    result = {
        "polyc_region": region_name if in_polyc else None,
        "after_polyc_region": after_region if after_polyc else None,
        "before_polyc_region": before_region if before_polyc else None,
        "has_polyc_stretch": has_polyc_stretch,
        "ref_context": ref_context,
        "sample_context": sample_context,
        "strand": strand,
        "flagged_for_filtering": should_filter,
    }

    return {k: v for k, v in result.items() if v is not None}


# ---------------------------------------------------------------------------
# Per-row variant calling
# ---------------------------------------------------------------------------


def _call_row_variants(
    ref_seq: str,
    sample_seq: str,
    start_ref: int,
    quality_scores: list[int] | None = None,
    strand: str = "plus",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Call variants by comparing reference and sample alignment sequences.

    Produces SNPs, insertions (``ref == "-"``) and deletions (``seq == "-"``)
    with insertion-aware position arithmetic and per-variant polyC context
    metadata.
    """
    variants: list[dict[str, Any]] = []
    insertion_count = 0

    for i in range(len(ref_seq)):
        ref_base = ref_seq[i]
        seq_base = sample_seq[i]

        if ref_base == seq_base:
            continue

        pos = start_ref + i - insertion_count

        if seq_base == "-":  # Deletion
            pos = normalize_position(math.floor(pos))
        elif ref_base == "-":  # Insertion
            insertion_count += 1
            pos = normalize_position(f"{start_ref + i - insertion_count}.1")
        else:  # SNP
            pos = normalize_position(int(pos))

        variant: dict[str, Any] = {"pos": pos, "ref": ref_base, "seq": seq_base}

        if quality_scores and i < len(quality_scores):
            if ref_base == "-":  # Insertion - no ref quality
                variant["quality"] = None
            else:
                variant["quality"] = quality_scores[i] if i < len(quality_scores) else None

        polyc_context = _analyze_polyc_context(pos, ref_seq, sample_seq, i, strand)
        if polyc_context:
            variant.update(polyc_context)

        variants.append(variant)

    snps, insertions, deletions = split_variant_types(variants)
    return snps, insertions, deletions


# ---------------------------------------------------------------------------
# Collect variants by position
# ---------------------------------------------------------------------------


def _collect_variants_by_position(  # noqa: C901
    rows: list[dict[str, Any]],
    quality_data: dict[str, list[int]] | None,
    regions: dict[str, list[int]],
    ref_seq: str,
) -> dict[Any, dict[str, Any]]:
    """Collect variant data from alignment rows, grouped by position.

    For each alignment row, calls per-row variants (SNPs + insertions +
    deletions), right-normalizes indels, and stores into a position-grouped
    dict with strand tracking and polyC ``flagged_for_filtering`` metadata.
    """
    variants_by_pos: dict[Any, dict[str, Any]] = defaultdict(
        lambda: {
            "ref": None,
            "seq": [],
            "file": [],
            "strand_bases": defaultdict(list),
            "strands": [],
            "flagged_strands": [],
            "strand_files": defaultdict(list),
            "qualities": [],
        }
    )

    for row in rows:
        strand = row.get("sstrand", "plus")
        seq_id = row.get("qseqid", "")
        quality_scores = quality_data.get(seq_id) if quality_data else None

        ref_align = row.get("forward_ref_seq", "")
        sample_align = row.get("forward_sample_seq", "")
        ref_start = int(row.get("forward_ref_start", 0))

        snps, insertions, deletions = _call_row_variants(ref_align, sample_align, ref_start, quality_scores, strand)
        norm_variants = snps + _right_norm_indel(insertions + deletions, ref_seq)

        for variant in norm_variants:
            pos = variant["pos"]

            # Only keep variants within regions of interest
            if get_hv_region_for_position(pos, regions) is None:
                continue

            if variants_by_pos[pos]["ref"] is None:
                variants_by_pos[pos]["ref"] = variant["ref"]

            variants_by_pos[pos]["seq"].append(variant["seq"])
            variants_by_pos[pos]["file"].append(seq_id)

            if variant.get("quality") is not None:
                variants_by_pos[pos]["qualities"].append(variant["quality"])

            variants_by_pos[pos]["strand_bases"][strand].append(variant["seq"])
            variants_by_pos[pos]["strand_files"][strand].append(seq_id)

            if strand not in variants_by_pos[pos]["strands"]:
                variants_by_pos[pos]["strands"].append(strand)

            # Track polyC-flagged strands
            if variant.get("flagged_for_filtering", False) and strand not in variants_by_pos[pos]["flagged_strands"]:
                variants_by_pos[pos]["flagged_strands"].append(strand)

            # Preserve polyC region info
            if "polyc_region" in variant:
                variants_by_pos[pos]["polyc_region"] = variant["polyc_region"]
                variants_by_pos[pos]["has_polyc_stretch"] = variant.get("has_polyc_stretch", False)
            if "after_polyc_region" in variant:
                variants_by_pos[pos]["after_polyc_region"] = variant["after_polyc_region"]
            if "before_polyc_region" in variant:
                variants_by_pos[pos]["before_polyc_region"] = variant["before_polyc_region"]

    return variants_by_pos


# ---------------------------------------------------------------------------
# Compute consensus variants
# ---------------------------------------------------------------------------


def _compute_consensus_variants(  # noqa: C901
    variants_by_pos: dict[Any, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Compute consensus variants from position-grouped variant data.

    Determines the consensus base per position from strand-aware frequency
    data and applies strand-specific polyC filtering:

    - Heterozygous single-strand variants are removed.
    - Heterozygous variants flagged by polyC ``flagged_for_filtering`` on
      only one strand are removed; the same variant flagged on both strands
      (homozygous) is retained.
    - Non-polyC heterozygous variants are retained for flagging annotation.
    """
    filtered_variants: list[dict[str, Any]] = []

    for pos in sorted(variants_by_pos.keys(), key=pos_sort_key):
        var_info = variants_by_pos[pos]

        plus_files = set(var_info["strand_files"].get("plus", []))
        minus_files = set(var_info["strand_files"].get("minus", []))

        counter_dict = Counter(var_info["seq"])
        max_count = max(counter_dict.values()) if counter_dict else 0
        bases_with_max_count = [base for base, count in counter_dict.items() if count == max_count]

        if not bases_with_max_count:
            continue

        consensus_base = ""
        if len(bases_with_max_count) > 1:
            if "-" in bases_with_max_count:
                bases_with_max_count.remove("-")
            consensus_base = nucleotide_to_iupac(bases_with_max_count)
        else:
            consensus_base = bases_with_max_count[0] if bases_with_max_count else ""

        # Check heterozygosity
        is_heterozygous = False
        if len(consensus_base) == 1:
            iupac_codes = "RYSWKMBDHVN"
            is_heterozygous = consensus_base in iupac_codes
        else:
            all_bases = Counter(var_info["seq"])
            is_heterozygous = len([c for c in all_bases.values() if c > len(var_info["seq"]) / 4]) > 1

        # Check strand specificity
        strand_specific = (bool(plus_files) and not minus_files) or (bool(minus_files) and not plus_files)

        # Apply filtering
        should_keep = True
        is_flagged = len(var_info.get("flagged_strands", [])) > 0

        if is_heterozygous and strand_specific:
            should_keep = False
        elif is_flagged and is_heterozygous:
            # flagged_strands is a subset of all_strands.
            all_strands = set(var_info["strands"])
            flagged_strands = set(var_info.get("flagged_strands", []))

            # Filter if heterozygous and only flagged on one strand.
            if (flagged_strands != all_strands and len(flagged_strands) < len(all_strands)) or len(
                all_strands
            ) < _HOMOZYGOUS_MIN_STRANDS:
                should_keep = False

        if not should_keep:
            continue

        qualities = var_info.get("qualities", [])

        variant: dict[str, Any] = {
            "pos": pos,
            "ref": var_info["ref"] or "N",
            "seq": consensus_base,
            "quality": qualities,
            "file": list(set(var_info["file"])),
        }

        # Preserve polyC region info for flagging annotation
        if "polyc_region" in var_info:
            variant["polyc_region"] = var_info["polyc_region"]
            variant["has_polyc_stretch"] = var_info.get("has_polyc_stretch", False)

        filtered_variants.append(variant)

    return filtered_variants


def _call_variants(
    rows: list[dict[str, Any]],
    quality_data: dict[str, list[int]] | None = None,
    ref_seq: str = "",
) -> list[dict[str, Any]]:
    """Call variants from all alignment rows and build consensus.

    Processes all alignment rows, calls variants position-by-position
    (including indels and polyC context), right-normalizes indels, and
    produces a consolidated variant list with strand-aware consensus
    filtering.

    Args:
        rows: BLASTN alignment rows (already strand-normalized).
        quality_data: Quality scores keyed by read ID, or None.
        ref_seq: Reference sequence string (rCRS), used for indel
            right-normalization and polyC stretch checks.

    Returns:
        List of raw variant dicts (before position filtering and conversion).
    """
    settings = get_settings()
    regions = settings.regions.REGIONS
    variants_by_pos = _collect_variants_by_position(rows, quality_data, regions, ref_seq)

    return _compute_consensus_variants(variants_by_pos)


# ---------------------------------------------------------------------------
# Interval computation
# ---------------------------------------------------------------------------


def _compute_intervals(
    rows: list[dict[str, Any]],
    regions: dict[str, list[int]] | None = None,
) -> dict[str, list[list[int]]]:
    """Compute overlapping intervals between predefined regions and sequenced regions.

    Args:
        rows: BLASTN alignment rows (already strand-normalized).
        regions: Genomic regions dict. If None, uses settings.

    Returns:
        Dict mapping region names to list of [start, end] interval pairs.
    """
    if regions is None:
        regions = get_settings().regions.REGIONS

    # Extract sequenced regions from alignment rows
    seq_regions: list[list[int]] = []
    for row in rows:
        start = row.get("forward_ref_start", 0)
        end = row.get("forward_ref_end", 0)
        if start and end:
            seq_regions.append([int(start), int(end)])

    # Merge overlapping intervals
    merged = merge_intervals(seq_regions)

    # Find overlaps with predefined regions
    intervals: dict[str, list[list[int]]] = {}
    for region, bounds in regions.items():
        overlap = get_overlap([bounds], merged)
        if overlap:
            intervals[region] = overlap

    return intervals


# ---------------------------------------------------------------------------
# Position filtering
# ---------------------------------------------------------------------------


def _filter_positions(
    variants: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Remove variants at filter positions.

    This is the REMOVE step — distinct from flagging. Variants at positions
    in the filter set are excluded from the final list entirely, unless
    they are also in the no-filter set.

    Args:
        variants: List of variant dicts.
        filter_positions: Set of position strings to filter. Loaded from
            rules/filter/filter.txt if None.
        no_filter_positions: Set of position strings that should never be
            filtered. Loaded from rules/filter/no_filter.txt if None.

    Returns:
        Filtered list of variant dicts.
    """
    filter_positions = _load_filter_positions()
    no_filter_positions = _load_no_filter_positions()

    filtered: list[dict[str, Any]] = []
    for v in variants:
        pos = v.get("pos", "")
        # Normalize for consistent set membership (matches PostProcess.is_special_position).
        pos_str = str(normalize_position(pos))
        pos_base_val = pos_base(pos)

        # no_filter overrides everything — never remove these positions.
        if pos_str in no_filter_positions:
            filtered.append(v)
            continue
        # Explicit filter set membership.
        if pos_str in filter_positions:
            logger.debug("Filtering variant at position {}", pos_str)
            continue
        # Dynamic 309 rule: any 309.x position is filtered unless in no_filter.
        if pos_base_val == 309 and pos_str not in no_filter_positions:  # noqa: PLR2004
            logger.debug("Filtering variant at position {} (309 dynamic rule)", pos_str)
            continue
        filtered.append(v)

    return filtered


# ---------------------------------------------------------------------------
# Build Variant objects
# ---------------------------------------------------------------------------


def _build_variant_objects(variants: list[dict[str, Any]]) -> list[Variant]:
    """Convert raw variant dicts to Pydantic Variant objects.

    Args:
        variants: List of variant dicts with pos, ref, seq, quality, files.

    Returns:
        List of Variant model instances.
    """
    variant_objects: list[Variant] = []
    for v in variants:
        pos = v.get("pos", 0)
        # Normalize position
        if isinstance(pos, str) and "." in pos:
            pos_norm = pos  # Keep insertion positions as strings
        else:
            pos_norm = int(pos) if isinstance(pos, (int, float)) else int(pos) if str(pos).isdigit() else pos

        # Build quality list with explicit int conversion (Phred scores are ints;
        # Variant.quality is list[int]).
        raw_quality = v.get("quality", [])
        if isinstance(raw_quality, list):
            quality_list: list[int] = [int(q) for q in raw_quality if q is not None]
        elif raw_quality is not None:
            quality_list = [int(raw_quality)]
        else:
            quality_list = []

        # Build files list with explicit str conversion
        raw_files = v.get("file", [])
        if isinstance(raw_files, list):
            files_list: list[str] = [str(f) for f in raw_files if f]
        elif raw_files:
            files_list = [str(raw_files)]
        else:
            files_list = []

        variant_objects.append(
            Variant(
                pos=pos_norm,
                ref=str(v.get("ref", "")),
                seq=str(v.get("seq", "")),
                quality=quality_list,
                files=files_list,
                peaks=None,
            )
        )
    return variant_objects


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def process(
    sample_id: str,
    input_path: str | Path,
    ref_path: str,
    batch_id: str | None = None,
) -> Sample:
    """Parse a BLASTN TSV file and produce a Sample object.

    This is the main ETL entry point. It reads BLASTN alignment output,
    calls variants, applies filtering and conversion rules, and produces
    a standardized Sample object.

    Args:
        sample_id: Sample identifier (LID from lab).
        input_path: Path to BLASTN TSV output file.
        ref_path: Path to reference FASTA file (rCRS).
        batch_id: Optional batch identifier.

    Returns:
        Sample with source_tool=Tool.BLASTN, variants, intervals,
        and flags. hv1/hv2/hv3 and statistics are None (filled by Batch).
    """
    settings = get_settings()
    regions = settings.regions.REGIONS

    # Load reference sequence
    ref_seq = str(SeqIO.read(ref_path, "fasta").seq)

    # Step 1: Parse BLASTN TSV and normalize strands
    rows = parse_blastn_tsv(input_path)
    rows = normalize_strand_all(rows)

    # Step 2: Load quality data
    quality_path = Path(str(input_path).replace(".concatenate.blastn", ".concatenate.quality"))
    quality_data = extract_quality_scores(quality_path) if quality_path.exists() else {}

    # Step 3: Call variants from alignment
    raw_variants = _call_variants(rows, quality_data, ref_seq)

    # Step 4: Compute intervals from alignment ranges
    intervals = _compute_intervals(rows, regions)

    # Step 5: Remove polyC-region artifacts; set polyC metadata on kept variants.
    filtered_variants = _handle_polyc_removal(raw_variants, ref_seq)

    # Step 6: Filter special positions (incl. the 309.x dynamic rule)
    filtered_variants = _filter_positions(filtered_variants)

    # Step 7: Apply conversion rules (YAML Add/Remove)
    filtered_variants = _apply_conversion_rules(filtered_variants)

    # Step 8: Apply region validation rules (requires no_filter positions)
    filtered_variants = _apply_region_validation_rules(filtered_variants)

    # Step 9: Build Variant objects
    variant_list = _build_variant_objects(filtered_variants)

    # Step 10: Flag variants (annotate, distinct from filter/remove)
    variant_flags_dict = flag_variants(variant_list)

    # Step 11: Merge polyC warning flags into variant_flags
    polyc_flags = _polyc_warning_flags(filtered_variants)
    for key, reasons in polyc_flags.items():
        if key in variant_flags_dict:
            for r in reasons:
                if r not in variant_flags_dict[key]:
                    variant_flags_dict[key].append(r)
        else:
            variant_flags_dict[key] = list(reasons)

    flagger = SampleFlagger(variant_list, intervals=intervals or None)
    _, analyze_reasons, _ = flagger.analyze()
    sample_flags_list = list(analyze_reasons)

    # Per-variant flags belong only in variant_flags; drop any sample-level
    # flag that duplicates a per-variant flag string.
    sample_flags_list = deduplicate_sample_flags(sample_flags_list, variant_flags_dict)

    # Step 12: Build and return Sample
    return Sample(
        sample_id=sample_id,
        variants=variant_list,
        source_tool=Tool.BLASTN,
        intervals=intervals or None,
        sample_flags=sample_flags_list or [],
        variant_flags=variant_flags_dict,
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
