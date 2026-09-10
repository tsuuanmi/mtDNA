"""Sample I/O: serialization and deserialization of Sample objects.

This module is the canonical source for converting between Sample objects and
their JSON representation. All Sample-to-JSON and JSON-to-Sample conversion
logic lives here.

- sample_to_dict()  — Sample → canonical JSON dict (grouped variants + flags + information)
- load_sample_batch() — JSON file → List[Sample]
- load_sample() — single per-sample JSON file → Sample
- filter_sample_by_regions() — tool-agnostic post-ETL region filter
- _normalize_keyed_sample() — grouped/flat variant normalization
- _parse_variant_dict() — single variant dict normalization
"""

import json
from pathlib import Path
from typing import Any

from Bio import SeqIO
from loguru import logger

from src.config import get_settings
from src.core.flagging import flag_variants, recompute_sample_flags
from src.core.models import Sample, Sequence, Tool, Variant
from src.core.variants import (
    normalize_position,
    pos_base,
    pos_sort_key,
    split_variant_types,
    statistic_variants,
    variant_to_dict,
)

# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------


def merge_intervals(intervals: list[list[int]]) -> list[list[int]]:
    """
    Merge overlapping intervals.

    Args:
        intervals: List of intervals

    Returns:
        List of merged intervals
    """
    # Handle empty input
    if not intervals:
        return []

    # Sort the array on the basis of start values of intervals.
    intervals.sort()
    stack = []
    # insert first interval into stack
    stack.append(intervals[0])
    for i in intervals[1:]:
        # Check for overlapping interval,
        # if interval overlap
        if stack[-1][0] <= i[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], i[-1])
        else:
            stack.append(i)

    return stack


# ---------------------------------------------------------------------------
# Reference/consensus sequence building helpers
# ---------------------------------------------------------------------------


def _load_reference_sequence(ref_path: str) -> dict:
    """Load reference sequence from FASTA file into position->base dict.

    Args:
        ref_path: Path to reference FASTA file

    Returns:
        Dictionary with positions (1-based int) as keys and bases as values.
        Returns empty dict if file cannot be read.
    """
    ref_dict = {}
    for record in SeqIO.parse(ref_path, "fasta"):
        for i, base in enumerate(record.seq):
            ref_dict[i + 1] = base  # 1-based positions
        return ref_dict
    logger.error(f"Failed to load reference sequence from {ref_path}")
    return {}


def _apply_variants_to_dict(ref_dict: dict, variants: dict | list) -> dict:
    """Apply variants to the reference dictionary.

    Args:
        ref_dict: Reference sequence dictionary (position -> base)
        variants: Variant data as Dict (position->{{ref,seq}}) or List ({{pos,ref,seq}}).
            Dict format is normalized to List internally.

    Returns:
        Modified dictionary with variants applied
    """
    # Normalize Dict input to List format
    # Dict format: {position: {"ref": ..., "seq": ...}} -> [{"pos": pos, "ref": ..., "seq": ...}]
    variants_list = [{"pos": pos, **val} for pos, val in variants.items()] if isinstance(variants, dict) else variants

    result = dict(ref_dict)

    for variant in variants_list:
        pos = variant["pos"]
        ref = variant.get("ref", "")
        seq = variant.get("seq", "")

        # Validate reference base matches the reference genome (Bug 2 fix)
        if ref not in {"-", ""}:
            ref_base = ref_dict.get(pos)
            if ref_base is not None and ref_base != ref:
                logger.warning(f"Reference base mismatch at position {pos}: expected {ref}, found {ref_base}")

        if ref == "-":
            # This is an insertion - normalize position for consistent key types
            # (int for bases, str for insertions like "217.1")
            result[normalize_position(pos)] = seq
        elif seq == "-":
            # This is a deletion - mark with '-'
            result[normalize_position(pos)] = "-"
        else:
            # This is a SNP - replace the base
            result[normalize_position(pos)] = seq

    return result


def _mark_non_region_as_n(seq_dict: dict, seq_regions: list) -> dict:
    """Mark positions outside specified regions as 'N'.

    Args:
        seq_dict: Sequence dictionary (position -> base)
        seq_regions: List of [start, end] intervals

    Returns:
        Modified dictionary with non-region positions marked as 'N'
    """
    if not seq_regions:
        return seq_dict

    result = dict(seq_dict)

    # Get all base positions (not insertions) in the dictionary.
    # With the always-str Position contract, keys can be int (from ref_dict) or str.
    # Only mark non-insertion positions as N; insertion positions stay as-is.
    base_positions = [k for k in result if not (isinstance(k, str) and "." in k)]

    for pos in base_positions:
        # Check if position is within any region
        pos_int = pos if isinstance(pos, int) else pos_base(pos)
        in_region = any(start <= pos_int <= end for start, end in seq_regions)
        if not in_region:
            result[pos] = "N"

    return result


def _extract_region(seq_dict: dict, start: int, end: int) -> str:
    """Extract a region from the sequence dictionary.

    Gets all keys in range, sorts naturally, removes deletions.
    Handles multi-base insertions as individual bases (Bug 5 fix).

    Args:
        seq_dict: Sequence dictionary (position -> base)
        start: Start position (inclusive, 1-based)
        end: End position (inclusive, 1-based)

    Returns:
        Extracted sequence as string
    """
    # Get all keys in the range [start, end]
    # This includes integer positions and fractional insertion positions
    keys_in_range = [k for k in seq_dict if start <= (pos_base(k) if isinstance(k, str) else k) <= end]

    # Sort naturally using pos_sort_key
    keys_in_range = sorted(keys_in_range, key=pos_sort_key)

    # Build sequence, skipping deletions
    # Bug 5 fix: Handle multi-base insertions as individual bases
    bases = []
    for key in keys_in_range:
        base = seq_dict[key]
        if base == "-":  # Skip deletion markers
            continue
        # Handle multi-character insertions (e.g., "ABCD" as 4 separate bases)
        bases.extend(base)

    return "".join(bases)


def generate_sequence(
    sample: Sample,
    ref_path: str | None = None,
) -> Sequence:
    """Generate consensus and reference sequences from a Sample model.

    Accepts a Sample model and returns a structured Sequence result object.
    Derives seq_regions from sample.intervals, regions from config defaults,
    and ref_path from config when not provided. ``intervals=None`` preserves
    legacy unspecified coverage; an explicit intervals mapping with no spans
    represents no trusted coverage and produces all-N sequences.

    Args:
        sample: A Sample model with variants, intervals, and sample_id.
        ref_path: Path to reference FASTA. Defaults to config ref directory.

    Returns:
        Sequence with reference_seq, consensus_seq, hv_seq_refs, hv_seqs,
        ref_dict, con_dict.
    """
    if ref_path is None:
        ref_path = str(get_settings().directories.ref / "rCRS.fasta")

    # Derive regions from config
    regions = get_settings().regions.REGIONS

    # Derive seq_regions from sample.intervals
    if sample.intervals:
        seq_regions = [[start, end] for intervals in sample.intervals.values() for start, end in intervals]
    else:
        seq_regions = []

    # Load reference sequence
    ref_dict = _load_reference_sequence(ref_path)
    if not ref_dict:
        return Sequence(reference_seq="", consensus_seq="")

    # Apply variants — convert List[Variant] to canonical dicts
    variants_dicts = [v for v in (variant_to_dict(var) for var in sample.variants) if v is not None]
    con_dict = _apply_variants_to_dict(ref_dict, variants_dicts)

    # Mark non-region positions as N
    if seq_regions:
        con_dict = _mark_non_region_as_n(con_dict, seq_regions)
        ref_dict = _mark_non_region_as_n(ref_dict, seq_regions)
    elif sample.intervals is not None:
        # Explicit empty coverage is NO_CALL, not an implicit reference call.
        # Use reference positions only: no untrusted insertion or deletion may
        # change the length of a sequence with no trusted evidence.
        ref_dict = dict.fromkeys(ref_dict, "N")
        con_dict = dict(ref_dict)

    # Extract HV region sequences
    hv_seqs = {}
    hv_seq_refs = {}
    for region, region_ranges in regions.items():
        # Handle both flat [start, end] and nested [[start, end], ...] formats
        if isinstance(region_ranges[0], list):
            region_start = region_ranges[0][0]
            region_end = region_ranges[0][-1]
        else:
            region_start = region_ranges[0]
            region_end = region_ranges[1]

        ref_bases = list(_extract_region(ref_dict, region_start, region_end))
        con_bases = list(_extract_region(con_dict, region_start, region_end))

        hv_seq_refs[region] = {
            "start": region_start,
            "end": region_end,
            "seq": ref_bases,
            "ranges": [[region_start, region_end]],
        }
        hv_seqs[region] = {
            "start": region_start,
            "end": region_end,
            "seq": con_bases,
            "ranges": [[region_start, region_end]],
        }

    # Convert dicts to strings
    ref_max_pos = max((k for k in ref_dict if isinstance(k, int)), default=16569)
    ref_str = "".join(ref_dict.get(i, "N") for i in range(1, ref_max_pos + 1))

    con_int_keys = [k for k in con_dict if isinstance(k, int)]
    con_max_pos = max(con_int_keys, default=ref_max_pos)
    con_str = "".join(con_dict.get(i, "N") for i in range(1, con_max_pos + 1) if con_dict.get(i) != "-")

    return Sequence(
        reference_seq=ref_str,
        consensus_seq=con_str,
        hv_seq_refs=hv_seq_refs,
        hv_seqs=hv_seqs,
        ref_dict=ref_dict,
        con_dict=con_dict,
    )


def sample_to_dict(sample: Sample, ref_path: str = "ref/rCRS.fasta") -> dict[str, Any]:
    """Serialize a single Sample into the canonical per-entry dict.

    Computes HV sequences via generate_sequence, writes them back to
    Sample.hv1/hv2/hv3, then serializes. Downstream consumers read from Sample
    instead of recomputing. Uses ``sample.intervals`` for region boundaries.

    The returned dict uses canonical batch field names (``HV1``/``HV2``/``HV3``,
    grouped variants) and omits ``sample_id`` (the caller provides it as the
    top-level key).
    """
    # Build per-sample intervals; Sample always has intervals set.
    intervals = dict(sample.intervals) if sample.intervals else {}

    hv1: str = ""
    hv2: str = ""
    hv3: str = ""
    no_snps = 0
    no_ins = 0
    no_dels = 0
    no_identicals = 0
    no_unread = 0

    computation_error = False
    try:
        result = generate_sequence(sample, ref_path)
        hv1 = "".join(result.hv_seqs.get("HV1", {}).get("seq", []))
        hv2 = "".join(result.hv_seqs.get("HV2", {}).get("seq", []))
        hv3 = "".join(result.hv_seqs.get("HV3", {}).get("seq", []))

        # Write computed sequences back to Sample
        sample.hv1 = hv1
        sample.hv2 = hv2
        sample.hv3 = hv3
    except (ValueError, KeyError, OSError) as e:
        logger.warning("generate_sequence failed for {}: {}", sample.sample_id, e)
        result = Sequence(reference_seq="", consensus_seq="")
        computation_error = True

    # statistic_variants expects Dict {pos: {ref, seq}}, not List[Variant]
    variants_dict = {str(v.pos): {"ref": v.ref, "seq": v.seq} for v in sample.variants}
    try:
        no_snps, no_ins, no_dels, no_identicals, no_unread = statistic_variants(
            variants_dict,
            result.hv_seq_refs,
            result.hv_seqs,
        )
    except (ValueError, KeyError, TypeError) as e:
        logger.warning("statistic_variants failed for {}: {}", sample.sample_id, e)
        computation_error = True

    # Inline grouping: normalize variants then split by type
    variant_list = sample.variants
    simple = [v for v in (variant_to_dict(v) for v in variant_list) if v is not None]
    snps, insertions, deletions = split_variant_types(simple)
    grouped_variants = {
        "snps": snps,
        "insertions": insertions,
        "deletions": deletions,
    }

    # Compute variant-level flags (general, applies to all tools).
    # Sample-level flags are tool-specific and set by ETL directly.
    variant_flags_dict = flag_variants(sample.variants)

    # Recompute core sample-level flags against THIS JSON's own variant set
    # so per-region JSONs don't carry stale cross-region flags (e.g.
    # "Complex: 459 deletion with 16192 issue" needs 459DEL + 16192T; an HV1-only
    # JSON lacks 459DEL and must not inherit the flag from the combined sample).
    # Tool-specific flags (not produced by SampleFlagger.analyze) are preserved.
    # This also handles "No 315.1 variant" (region-level gating: fires only when
    # the region containing position 315 — HV2 — was analyzed).
    sample_flags = recompute_sample_flags(
        sample.sample_flags,
        sample.variants,
        sample.intervals,
        variant_flags_dict,
    )

    sample_dict: dict[str, Any] = {
        "no_snps": no_snps,
        "no_ins": no_ins,
        "no_dels": no_dels,
        "no_identicals": no_identicals,
        "no_unread": no_unread,
        "HV1": hv1,
        "HV2": hv2,
        "HV3": hv3,
        "variants": grouped_variants,
        "intervals": intervals,
        "sample_flags": sample_flags,
        "variant_flags": variant_flags_dict,
        "information": {
            **({"batch_id": sample.batch_id} if sample.batch_id else {}),
            **({"source_tool": sample.source_tool.value} if sample.source_tool else {}),
            **({k: v for k, v in sample.information.items() if k != "batch"} if sample.information else {}),
            **({"computation_error": computation_error} if computation_error else {}),
        },
    }

    return sample_dict


def _normalize_keyed_sample(raw: dict[str, Any]) -> dict[str, Any]:
    """Normalize a keyed sample entry from grouped variant format to flat list.

    Converts {snps: [...], insertions: [...], deletions: [...]}
    into a single flat list of {pos, ref, seq, ...} dicts with a variant_type tag.
    """
    variants_raw = raw.get("variants", {})

    # Type-grouped format: {snps: [...], insertions: [...], deletions: [...]}
    # Region-grouped format: {HV1: [...], HV2: [...], HV3: [...]}
    if isinstance(variants_raw, dict):
        flat: list[dict[str, Any]] = []
        if "snps" in variants_raw:
            for group_name in ("snps", "insertions", "deletions"):
                for entry in variants_raw.get(group_name, []):
                    norm_entry = dict(entry)
                    norm_entry["variant_type"] = group_name
                    flat.append(norm_entry)
        elif any(k in variants_raw for k in ("HV1", "HV2", "HV3")):
            for region in ("HV1", "HV2", "HV3"):
                flat.extend(dict(entry) for entry in variants_raw.get(region, []))
        else:
            flat = []
        if flat:
            variants_raw = flat

    # Already a flat list: [{pos, ref, seq, ...}]
    if isinstance(variants_raw, list):
        return {**raw, "variants": variants_raw}

    return {**raw, "variants": []}


def _parse_variant_dict(v: dict[str, Any], sample_id: str) -> dict[str, Any]:
    """Normalize a single variant dict for Sample model validation."""
    result = dict(v)
    pos = result.get("pos", result.get("position"))
    if pos is None:
        logger.warning(f"Variant missing position in {sample_id}: {v}")
        return result

    result["position"] = pos
    result.setdefault("ref", result.get("ref", "-"))
    result.setdefault("alt", result.get("seq", result.get("alt", "-")))

    # Preserve both 'seq' and 'alt' for backward compat
    if "seq" in result and "alt" not in result:
        result["alt"] = result["seq"]
    if "alt" in result and "seq" not in result:
        result["seq"] = result["alt"]

    result.setdefault("files", [])

    # Normalize quality: scalar → [scalar], list stays list, None/missing → []
    quality = result.get("quality")
    if quality is None:
        result["quality"] = []
    elif isinstance(quality, (int, float)):
        result["quality"] = [float(quality)]
    # list stays as-is (already a list of floats)

    return result


def load_sample_entries(data: dict) -> list[Sample]:
    """Load a dict of {sample_id: raw_entry} into a list of Sample objects."""
    samples: list[Sample] = []
    for sample_id, raw in data.items():
        entry = _normalize_keyed_sample(raw)
        entry["sample_id"] = sample_id

        # Map batch JSON keys (HV1/HV2/HV3) to Sample field names (hv1/hv2/hv3)
        for json_key, field_name in [("HV1", "hv1"), ("HV2", "hv2"), ("HV3", "hv3")]:
            if json_key in entry and field_name not in entry:
                entry[field_name] = entry.pop(json_key)

        variant_dicts = entry.get("variants", [])
        entry["variants"] = [_parse_variant_dict(v, sample_id) for v in variant_dicts]

        # source_tool and batch_id may be top-level (legacy) or inside information
        information = entry.get("information") or {}
        source_tool_val = entry.get("source_tool") or information.get("source_tool")
        if isinstance(source_tool_val, str) and source_tool_val:
            try:
                entry["source_tool"] = Tool(source_tool_val)
            except ValueError:
                logger.warning(f"Unknown source_tool '{source_tool_val}' for {sample_id}, marking as unknown")
                entry["source_tool"] = Tool.UNKNOWN
        else:
            # No tool attribution in the batch: surface as "Unknown" rather
            # than silently defaulting to a concrete tool.
            entry["source_tool"] = Tool.UNKNOWN
        batch_id_val = entry.get("batch_id") or information.get("batch_id")
        if batch_id_val:
            entry["batch_id"] = batch_id_val

        samples.append(Sample.model_validate(entry))

    return samples


def load_sample_batch(path: str | Path) -> list[Sample]:
    """Load a batch JSON file into a list of Sample objects."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return load_sample_entries(data)


def load_sample(path: str | Path) -> Sample:
    """Load a single per-sample JSON file into a Sample object.

    Tool outputs write one standard Sample entry per file at
    ``<json_dir>/<sample_id>/<sample_id>.json``; the file stem is used as the
    sample_id. Variants may be grouped by type (``snps``/``insertions``/
    ``deletions``) or by region (``HV1``/``HV2``/``HV3``); both are normalized.
    """
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    return load_sample_entries({path.stem: data})[0]


# ---------------------------------------------------------------------------
# Region-filtered sample helper
# ---------------------------------------------------------------------------

# Canonical region name set — all filtering uses these keys.
ALL_REGIONS = {"HV1", "HV2", "HV3"}


def _position_in_any_interval(pos: int | str, intervals: list[list[int]]) -> bool:
    """Return True if *pos* falls inside any [start, end] pair in *intervals*.

    For insertion positions (e.g. "309.1"), the integer base (309) is used
    for the range check, matching how the merger validates variant coverage.
    """
    base = pos_base(pos)
    return any(start <= base <= end for start, end in intervals)


def _validate_regions(keep_regions: list[str]) -> None:
    """Raise ValueError if keep_regions contains invalid region names."""
    invalid = set(keep_regions) - ALL_REGIONS
    if invalid:
        msg = f"Unknown region name(s): {invalid}. Must be subset of {ALL_REGIONS}"
        raise ValueError(msg)


def _kept_intervals(
    sample: Sample,
    keep_regions: list[str],
) -> dict[str, list[list[int]]] | None:
    """Return intervals dict for the kept regions, or None if none have data."""
    src_intervals = sample.intervals or {}
    kept: dict[str, list[list[int]]] = {}
    for region in keep_regions:
        region_ivs = src_intervals.get(region)
        if region_ivs:
            kept[region] = region_ivs
    return kept or None


def _filter_variants_by_intervals(
    variants: list[Variant],
    intervals: dict[str, list[list[int]]],
) -> list[Variant]:
    """Keep only variants whose position falls within any interval."""
    result: list[Variant] = []
    for v in variants:
        for ivs in intervals.values():
            if _position_in_any_interval(v.pos, ivs):
                result.append(v)
                break
    return result


def _filter_variant_flags(
    variant_flags: dict[str, list[str]],
    kept_positions: set[str],
) -> dict[str, list[str]]:
    """Keep only variant_flags entries whose position is in kept_positions."""
    kept: dict[str, list[str]] = {}
    for key, reasons in variant_flags.items():
        pos_part = key.split("|")[0] if "|" in key else key
        if pos_part in kept_positions:
            kept[key] = reasons
    return kept


def filter_sample_by_regions(
    sample: Sample,
    keep_regions: list[str],
) -> Sample | None:
    """Return a new Sample containing only variants and intervals for *keep_regions*.

    This is a **tool-agnostic** post-ETL filter.  ETL modules stay pure —
    region scoping is applied after they have produced a full Sample.

    Args:
        sample: Source Sample (typically produced by an ETL module).
        keep_regions: List of region names to keep (e.g. ``["HV1"]``,
            ``["HV2", "HV3"]``).  Must be a subset of ``{"HV1", "HV2", "HV3"}``.

    Returns:
        A new Sample with only the kept region variants and intervals, or
        ``None`` when none of the requested regions have data (empty intervals
        or no matching variants).

    Behaviour:
        - Variants are kept when their position falls within **any** interval
          of a kept region.
        - The returned ``Sample.intervals`` dict contains **only** the kept
          region keys (Shape A — empty regions are omitted entirely).
        - HV sequence strings (``hv1``, ``hv2``, ``hv3``) for dropped regions
          are set to ``None``; kept regions retain their value.
        - ``sample_flags``, ``variant_flags``, ``information``, and
          ``batch_id`` are carried forward unchanged.
        - ``source_tool`` is preserved (filtering does **not** change tool
          attribution — only merging does).

    """
    _validate_regions(keep_regions)

    if not keep_regions:
        return None

    kept_intervals = _kept_intervals(sample, keep_regions)
    if not kept_intervals:
        return None

    filtered_variants = _filter_variants_by_intervals(sample.variants, kept_intervals)

    kept_hv = {
        "hv1": sample.hv1 if "HV1" in kept_intervals else None,
        "hv2": sample.hv2 if "HV2" in kept_intervals else None,
        "hv3": sample.hv3 if "HV3" in kept_intervals else None,
    }

    kept_positions = {str(v.pos) for v in filtered_variants}
    kept_variant_flags = _filter_variant_flags(sample.variant_flags, kept_positions)

    return Sample(
        sample_id=sample.sample_id,
        variants=filtered_variants,
        source_tool=sample.source_tool,
        intervals=kept_intervals,
        sample_flags=sample.sample_flags,
        variant_flags=kept_variant_flags,
        information=sample.information,
        batch_id=sample.batch_id,
        hv1=kept_hv["hv1"],
        hv2=kept_hv["hv2"],
        hv3=kept_hv["hv3"],
        no_snps=sample.no_snps,
        no_ins=sample.no_ins,
        no_dels=sample.no_dels,
        no_identicals=sample.no_identicals,
        no_unread=sample.no_unread,
    )
