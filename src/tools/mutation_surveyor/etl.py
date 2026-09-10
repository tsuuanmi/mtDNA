"""Mutation Surveyor ETL module — lid_merge xlsx → Sample.

Adapter layer: reads lid_merge output xlsx, parses per-LID variant tokens
into Variant objects and builds Sample objects with flags and intervals.

ETL is pure: it populates only variants, intervals, and flags. HV sequences
and statistics are computed by Batch at aggregation time.

Structured MSFlag objects (via combined_flags) route flags into
sample_flags vs variant_flags based on their FlagLevel. When combined_flags
are not available for a LID, flags default to empty lists.
"""

import re
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Position, Sample, Tool, Variant
from src.tools.mutation_surveyor.flagging import (
    MSFlag,
    build_region_status_flags,
    sample_flags_from,
    variant_flags_from,
)
from src.tools.mutation_surveyor.utils import (
    UNORDERED_MAP,
    parse_variant_tokens,
    safe_str,
)

# ---------------------------------------------------------------------------
# rCRS reference loader
# ---------------------------------------------------------------------------


def _load_rcrs(ref_path: str) -> dict[int, str]:
    """Load rCRS reference into {position: base} dict (1-based)."""
    for record in SeqIO.parse(ref_path, "fasta"):
        return {i + 1: str(base) for i, base in enumerate(record.seq)}
    logger.error("Failed to load rCRS reference from {}", ref_path)
    return {}


# ---------------------------------------------------------------------------
# Token → Variant mapping
# ---------------------------------------------------------------------------


def _token_to_variant(
    token: str,
    pos: Position,
    allele: str,
    ref_dict: dict[int, str],
    quality: list[float] | None = None,
    files: list[str] | None = None,
) -> Variant | None:
    """Convert a parsed variant token into a Variant object.

    Token types:
        SNP:        "73G"           → Variant(pos=73, ref=rCRS[73], seq="G")
        Insertion:  "309.1C"        → Variant(pos="309.1", ref="-", seq="C")
        Deletion:   "16189DEL"      → Variant(pos=16189, ref=rCRS[16189], seq="-")
        PHP/IUPAC:  "73R"           → Variant(pos=73, ref=rCRS[73], seq="R")
        Het SNP:    "309C_het"      → Variant(pos=309, ref=rCRS[309], seq="Y")
                    where seq is IUPAC({ref, allele}) so core flag_variants()
                    detects "Heteroplasmy at 309" automatically.
        Het del:    "16189het_DEL"  → Variant(pos=16189, ref=rCRS[16189], seq="DEL")
                    Het deletions keep seq="DEL" (not "-") so the het signal
                    is preserved; core flag_variants() detects
                    "Heteroplasmy at {pos}" via the seq="DEL" check.

    Args:
        token: Raw token string (e.g. "73G", "309.1C", "16189DEL").
        pos: Position from parse_variant_tokens (int for bases, str for insertions).
        allele: Allele string from parse_variant_tokens.
        ref_dict: rCRS reference dict for ref base lookup.
        quality: Optional list of per-trace quality scores.
        files: Optional list of source trace filenames.

    Returns:
        Variant or None if token cannot be mapped.
    """
    is_insertion = "." in token and not allele.upper().startswith("DEL")
    is_deletion = allele.upper() == "DEL"
    is_het = "_het" in allele.lower() or allele.lower().startswith("het_")

    # Parse position
    if is_insertion:
        # e.g. "309.1C" → pos_str = "309.1"
        m = re.match(r"^(\d+\.\d+)", token)
        pos_str = m.group(1) if m else str(pos)
        seq_val = allele.replace("_het", "").replace("_HET", "").removeprefix("het_").removeprefix("HET_")
        if is_het:
            # Het insertions: lowercase seq signals low-confidence (ref="-",
            # no IUPAC); core catches this via has_lowercase_variants().
            seq_val = seq_val.lower()
        return Variant(
            pos=pos_str,
            ref="-",
            seq=seq_val,
            quality=quality or [],
            files=files or [],
            peaks=None,
        )

    base_pos = int(pos)
    ref_base = ref_dict.get(base_pos, "?")

    if is_deletion:
        # Het deletions keep seq="DEL" (not "-") so the het signal is
        # preserved for downstream consumers; core flag_variants still
        # detects "Heteroplasmy at {pos}" via the seq="DEL" check.
        seq_val = "DEL" if is_het else "-"
        return Variant(
            pos=base_pos,
            ref=ref_base,
            seq=seq_val,
            quality=quality or [],
            files=files or [],
            peaks=None,
        )

    # SNP or PHP/IUPAC
    seq_val = allele.replace("_het", "").replace("_HET", "").removeprefix("het_").removeprefix("HET_")

    if is_het and len(seq_val) == 1 and ref_base != "?":
        # Het SNP: compute IUPAC ambiguity code from {ref_base, seq_val}
        # so core flag_variants() detects "Heteroplasmy at {pos}".
        iupac_key = "".join(sorted(ref_base.upper() + seq_val.upper()))
        iupac_code = UNORDERED_MAP.get(iupac_key)
        if iupac_code:
            seq_val = iupac_code
        # If IUPAC lookup fails (unknown combination), seq_val stays as-is;
        # core will still see it as a variant at this position.

    return Variant(
        pos=base_pos,
        ref=ref_base,
        seq=seq_val,
        quality=quality or [],
        files=files or [],
        peaks=None,
    )


def _to_interval(bounds: list[int]) -> list[list[int]]:
    """Wrap a flat [start, end] bounds list into [[start, end]]."""
    return [[bounds[0], bounds[1]]]


def _parse_ranges_from_text(range_text: str) -> dict[str, list[list[int]]]:  # noqa: C901
    """Parse 'Sample profile - Range' text into intervals dict.

    E.g. "73-340 438-576 16024-16365" →
        {"HV2": [[73, 340]], "HV3": [[438, 576]], "HV1": [[16024, 16365]]}

    Falls back to default full-region intervals if parsing fails.
    """
    regions = get_settings().regions.REGIONS
    intervals: dict[str, list[list[int]]] = {}

    if not range_text or range_text.strip() in ("", "nan", "None"):
        for name, bounds in regions.items():
            intervals[name] = _to_interval(bounds)
        return intervals

    # Parse space-separated "start-end" pairs
    pairs: list[list[int]] = []
    for m in re.finditer(r"(\d+)\s*-\s*(\d+)", range_text):
        start, end = int(m.group(1)), int(m.group(2))
        if start > end:
            start, end = end, start
        pairs.append([start, end])

    if not pairs:
        for name, bounds in regions.items():
            intervals[name] = _to_interval(bounds)
        return intervals

    # Assign pairs to HV regions by overlap
    for name, bounds in regions.items():
        for region_bound in _to_interval(bounds):
            lo, hi = region_bound[0], region_bound[1]
            matched = [p for p in pairs if p[0] <= hi and p[1] >= lo]
            if matched:
                intervals[name] = matched
                break

    # Fill missing regions with defaults
    for name, bounds in regions.items():
        if name not in intervals:
            intervals[name] = _to_interval(bounds)

    return intervals


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def process(  # noqa: C901, PLR0912, PLR0915
    sample_id: str,
    input_path: Path,
    ref_path: str = "ref/rCRS.fasta",
    batch_id: str | None = None,
    combined_flags: dict[str, list[MSFlag]] | None = None,
) -> Sample:
    """Process a single LID from a lid_merge xlsx into a Sample object.

    When combined_flags (structured MSFlag objects per LID) are provided,
    flags are routed into sample_flags and variant_flags based on their
    FlagLevel. When combined_flags is None or the LID is not found,
    flags default to empty lists.

    Args:
        sample_id: LID to extract.
        input_path: Path to lid_merge xlsx file.
        ref_path: Path to rCRS reference FASTA.
        batch_id: Optional batch identifier.
        combined_flags: Optional dict mapping LID → list[MSFlag] from pipeline.

    Returns:
        Sample with parsed variants, intervals, source provenance, and flags.
    """
    input_path = Path(input_path)
    ref_dict = _load_rcrs(ref_path)

    # Read LID sheet
    df = pd.read_excel(input_path, sheet_name="LID")
    df.columns = [str(c).strip() for c in df.columns]

    # Find matching row
    lid_rows = df[df["LID"].astype(str).str.strip() == sample_id]
    if lid_rows.empty:
        logger.warning("LID {} not found in {}", sample_id, input_path)
        return Sample(
            sample_id=sample_id,
            variants=[],
            source_tool=Tool.MUTATION_SURVEYOR,
            intervals=None,
            batch_id=batch_id,
            sample_flags=[],
            variant_flags={},
            information=None,
            hv1=None,
            hv2=None,
            hv3=None,
            no_snps=None,
            no_ins=None,
            no_dels=None,
            no_identicals=None,
            no_unread=None,
        )

    row = lid_rows.iloc[0]

    # Parse variant tokens
    variants_col = "Sample profile - Variants"
    variant_text = safe_str(row.get(variants_col, ""))
    parsed = parse_variant_tokens(variant_text)

    # Parse per-variant quality scores and source files
    quality_text = safe_str(row.get("Sample profile - Quality", ""))
    quality_tokens = quality_text.split() if quality_text else []
    files_text = safe_str(row.get("Sample profile - Files", ""))
    source_files = files_text.split() if files_text else []

    variant_list: list[Variant] = []

    for i, (pos, allele, token) in enumerate(parsed):
        q_val: list[float] = []
        if i < len(quality_tokens):
            q_tok = quality_tokens[i]
            if q_tok not in ("", "."):
                try:
                    q_val = [float(q_tok)]
                except ValueError:
                    q_val = []
        v = _token_to_variant(token, pos, allele, ref_dict, quality=q_val, files=source_files)
        if v is not None:
            variant_list.append(v)

    # Parse ranges → intervals
    range_text = safe_str(row.get("Sample profile - Range", ""))
    intervals = _parse_ranges_from_text(range_text)

    # Route flags into sample_flags vs variant_flags
    if combined_flags and sample_id in combined_flags:
        flags_for_lid = combined_flags[sample_id]
        sample_flags = sample_flags_from(flags_for_lid)
        variant_flags_dict = variant_flags_from(flags_for_lid)
    else:
        sample_flags = []
        variant_flags_dict = {}

    # Per-variant flags via core flagging (polyC, het, etc.)
    core_variant_flags = flag_variants(variant_list)
    # Merge core variant flags into variant_flags_dict
    for key, vals in core_variant_flags.items():
        variant_flags_dict.setdefault(key, []).extend(vals)
        # Deduplicate while preserving order
        seen = set()
        unique = []
        for v in vals:
            if v not in seen:
                seen.add(v)
                unique.append(v)
        variant_flags_dict[key] = unique

    # Information metadata
    information: dict[str, Any] = {}
    # Region status flags from per-region review classes
    hv1_review_class = safe_str(row.get("HV1 Review Class", ""))
    hv23_review_class = safe_str(row.get("HV2-3 Review Class", ""))
    region_status_flags = build_region_status_flags(hv1_review_class, hv23_review_class)
    sample_flags.extend(f.name for f in region_status_flags)

    # Run core sample-level flagging on the full variant set. Only
    # per-variant flag_variants() was called above, so sample-level flags
    # (e.g. "Complex: 459 deletion with 16192 issue") were never computed.
    if variant_list:
        flagger = SampleFlagger(variant_list, intervals=intervals or None)
        _, core_sample_reasons, _ = flagger.analyze()
        for reason in core_sample_reasons:
            if reason not in sample_flags:
                sample_flags.append(reason)

    # Per-variant flags belong only in variant_flags; drop any sample-level
    # flag that duplicates a per-variant flag string.
    sample_flags = deduplicate_sample_flags(sample_flags, variant_flags_dict)

    well = safe_str(row.get("Well", ""))
    if well:
        information["well"] = well

    return Sample(
        sample_id=sample_id,
        variants=variant_list,
        source_tool=Tool.MUTATION_SURVEYOR,
        intervals=intervals or None,
        batch_id=batch_id,
        sample_flags=sample_flags,
        variant_flags=variant_flags_dict,
        information=information or None,
        hv1=None,
        hv2=None,
        hv3=None,
        no_snps=None,
        no_ins=None,
        no_dels=None,
        no_identicals=None,
        no_unread=None,
    )
