"""Tracy ETL module — Tracy decompose JSON → Sample.

Pure ETL: reads a single Tracy decompose JSON output file and produces a
Sample object with variants, intervals, and flags. No file I/O beyond
reading the input JSON. No subprocess calls. No side effects.

The entry point is ``process()`` which takes a single JSON file and returns
a Sample. For batch processing, use ``src.tools.tracy.pipeline.process_batch``.

Variant transformation logic (position-specific, polyC, and strand transforms)
lives in ``src.tools.tracy.transforms``.
"""

import json
from pathlib import Path
from typing import Any, NamedTuple

from loguru import logger

from src.config import get_settings
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Sample, Tool, Variant
from src.core.sample import merge_intervals
from src.core.variants import (
    is_position_in_intervals,
    nucleotide_to_iupac,
    pos_base,
    pos_sort_key,
)
from src.tools.tracy.transforms import (
    _create_variant,
    _CreateParams,
    _PositionParams,
    apply_all_transformations,
    update_variant_positions,
)
from src.tools.tracy.utils import (
    POS_73,
    POS_16000,
    POS_16189,
    POS_16193,
    PrimerTypeInfo,
    alignment_bounds,
    detect_heteroplasmy,
    detect_primer_type,
    detect_variant_conditions,
    extract_peak_data,
    normalize_ref_positions,
    validate_peak_quality,
)

# ---------------------------------------------------------------------------
# Helper: load Tracy decompose JSON
# ---------------------------------------------------------------------------


def _load_tracy_json(input_path: Path) -> dict[str, Any]:
    """Load a Tracy decompose JSON file.

    Args:
        input_path: Path to the Tracy decompose JSON file.

    Returns:
        Parsed JSON data dict.
    """
    with input_path.open() as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Parameter dataclasses to reduce function argument counts
# ---------------------------------------------------------------------------


class _CallParams(NamedTuple):
    """Parameters for calling a variant at a single alignment position."""

    data: dict[str, Any]
    deletions: int
    filename: str
    heteroplasmy_threshold: float
    sample: str | None


class _FilterConfig(NamedTuple):
    """Configuration for variant filtering."""

    quality_threshold: int
    min_peak_value: int
    pratio: float
    sample: str | None
    variant_regions: dict[str, list[int]]


class _ProcessConfig(NamedTuple):
    """Configuration parameters for the Tracy ETL process function."""

    quality_threshold: int
    min_peak_value: int
    pratio: float
    heteroplasmy_threshold: float
    variant_regions: dict[str, list[int]] | None = None
    batch_id: str | None = None
    sample: str | None = None


# ---------------------------------------------------------------------------
# Helper: call a variant at a single alignment index
# ---------------------------------------------------------------------------


def _call_variant_at_index(index: int, params: _CallParams) -> dict[str, Any] | None:
    """Call a variant at a single alignment index.

    Peak data is extracted at ``peak_index = index - deletions`` (the alignment
    index minus deletions seen so far). ``index`` keeps the raw alignment index;
    ``peak_index`` keeps the corrected basecall index.
    """
    data = params.data
    alt1 = data["alt1align"][index]
    ref1 = data["ref1align"][index]

    # Skip if alt matches ref (reference position)
    if alt1 == ref1:
        return None

    pos_raw = index  # alignment index; transformed to absolute pos later
    peak_index = index - params.deletions  # corrected basecall index

    # Extract peak data at the corrected basecall index
    peak_a, peak_c, peak_g, peak_t, quality = extract_peak_data(data, peak_index)
    peaks: list[float | None] = [peak_a, peak_c, peak_g, peak_t]

    variant: dict[str, Any] = {
        "pos": pos_raw,
        "ref": ref1,
        "seq": alt1,
        "peaks": peaks,
        "index": index,
        "peak_index": peak_index,
        "quality": quality,
    }

    # Detect heteroplasmy whenever an alternate base is present, so insertions
    # (ref="-", alt=base) also get IUPAC heteroplasmy calls.
    if alt1 != "-":
        detected_base = detect_heteroplasmy(peaks, pos_raw, params.heteroplasmy_threshold, params.sample)
        if detected_base != "N":
            variant["seq"] = detected_base

    return variant


# ---------------------------------------------------------------------------
# Helper: filter raw variants
# ---------------------------------------------------------------------------


def _should_skip_variant(
    variant: dict[str, Any],
    config: _FilterConfig,
) -> bool:
    """Check if a variant should be skipped during filtering.

    Filter conditions (in order):
    1. Skip variants marked for removal (``variant.get("remove")``)
    2. Skip 'N' bases
    3. Skip variants outside regions (``is_position_in_intervals``)
    4. Skip HV1 deletions except 16189 and 16193
    5. Skip low quality variants (except position 73)
    6. Validate peak quality

    Index field removal happens later in ``process()`` after filtering.
    """
    # Skip variants marked for removal
    if variant.get("remove", False):
        return True

    # Skip 'N' bases
    if variant["seq"] == "N":
        return True

    # Skip variants outside regions
    if not is_position_in_intervals(variant["pos"], [(int(v[0]), int(v[1])) for v in config.variant_regions.values()]):
        return True

    # Skip HV1 deletions except 16189 and 16193
    if (
        pos_base(variant["pos"]) > POS_16000
        and pos_base(variant["pos"]) not in [POS_16189, POS_16193]
        and variant["seq"] == "-"
    ):
        return True

    # Skip low quality variants (except position 73)
    if (
        variant.get("quality") is not None
        and variant["quality"] < config.quality_threshold
        and pos_base(variant["pos"]) != POS_73
    ):
        return True

    # Validate peak quality
    return not validate_peak_quality(variant, config.min_peak_value, config.pratio, config.sample)


# ---------------------------------------------------------------------------
# Helper: build Variant objects from filtered dicts
# ---------------------------------------------------------------------------


def _build_variant_objects(filtered_variants: list[dict[str, Any]], filename: str) -> list[Variant]:
    """Build Variant model objects from filtered variant dicts.

    Wraps raw dicts in ``Variant`` model objects for type safety and
    validation. Peaks and quality are wrapped in outer lists to match the
    ``Variant`` model's ``list[list[int | None]]`` and ``list[int]`` type
    annotations. Peak heights and Phred quality scores are integer-valued;
    the ``int`` model typing (Pydantic v2) coerces integer-valued floats and
    rejects fractional floats at construction.
    """
    variant_list: list[Variant] = []
    for variant in filtered_variants:
        try:
            # ETL produces flat peaks [A, C, G, T] and scalar quality;
            # Variant model expects list[list] and list[int].
            raw_peaks = variant.get("peaks", [])
            raw_quality = variant.get("quality", [])
            wrapped_peaks = [raw_peaks] if raw_peaks and not isinstance(raw_peaks[0], list) else raw_peaks
            wrapped_quality: list[int] = (
                [int(raw_quality)] if isinstance(raw_quality, (int, float)) else [int(q) for q in raw_quality]
            )

            v = Variant(
                pos=variant["pos"],
                ref=variant["ref"],
                seq=variant["seq"],
                files=[filename],
                peaks=wrapped_peaks or None,
                quality=wrapped_quality,
            )
            variant_list.append(v)
        except (ValueError, KeyError) as e:
            logger.warning("Skipping invalid variant: {}", e)
            continue

    return variant_list


# ---------------------------------------------------------------------------
# Helper: compute intervals for each file
# ---------------------------------------------------------------------------


def _compute_intervals(
    file_intervals: dict[str, list[list[int]]],
    variant_regions: dict[str, list[int]],
) -> dict[str, list[list[int]]]:
    """Compute coverage intervals per file.

    Takes ``file_intervals`` and ``variant_regions`` as parameters and uses
    ``merge_intervals()`` from ``src.core.sample`` for interval tracking.
    """
    intervals: dict[str, list[list[int]]] = {}

    for region_name, region_bounds in variant_regions.items():
        region_start, region_end = int(region_bounds[0]), int(region_bounds[1])
        intervals[region_name] = []

        for file_interval in file_intervals.values():
            for start, end in file_interval:
                overlap_start = max(start, region_start)
                overlap_end = min(end, region_end)

                if overlap_start <= overlap_end:
                    intervals[region_name].append([overlap_start, overlap_end])

        # Merge overlapping intervals
        if intervals[region_name]:
            intervals[region_name] = merge_intervals(intervals[region_name])
        else:
            intervals[region_name] = []

    return intervals


def merge_variants(variants: list[Variant]) -> list[Variant]:
    """Merge variants from multiple Tracy reads into one IUPAC consensus per (pos, ref).

    Variants sharing the same ``(pos, ref)`` are collapsed into a single
    variant whose ``seq`` is the IUPAC ambiguity code for the observed bases
    (e.g. C>T and C>A collapse to C>W). IUPAC encoding also unifies identical
    calls (T + T -> "T"), so there is no separate deduplication pass: every
    variant at a position is merged via IUPAC, not deduplicated by exact
    ``(pos, ref, seq)`` first. ``files``/``peaks``/``quality`` are combined
    from every source in first-appearance order.

    Identical deletion calls at the same ``(pos, ref)`` (e.g. A>- from HV2F
    and HV3R) collapse to a single deletion variant, combining
    ``files``/``peaks``/``quality``. Only mixed deletion/substitution
    conflicts at the same position are kept separate, since IUPAC does not
    cover deletions.
    """
    # Group by (pos, ref) — one consensus per position. No dedup pass:
    # IUPAC encoding subsumes exact-match collapsing (T + T -> "T").
    by_pos_ref: dict[Any, list[Variant]] = {}
    for v in variants:
        by_pos_ref.setdefault((pos_sort_key(v.pos), v.ref), []).append(v)

    merged: list[Variant] = []
    for group in by_pos_ref.values():
        if len(group) == 1:
            merged.append(group[0])
            continue
        seqs = [v.seq for v in group]
        seq_set = set(seqs)
        has_del = "-" in seq_set
        # Mixed deletion/substitution at the same (pos, ref): keep separate,
        # since IUPAC does not cover deletions. Identical calls (including
        # all-deletion groups, e.g. A>- from HV2F + HV3R) collapse to one.
        if has_del and len(seq_set) > 1:
            merged.extend(group)
            continue
        consensus = group[0].model_copy()
        # All-deletion group keeps seq "-"; otherwise IUPAC-encode the bases.
        consensus.seq = "-" if has_del else nucleotide_to_iupac(seqs)
        consensus.files = []
        consensus.quality = []
        consensus.peaks = None
        for v in group:
            consensus.files = consensus.files + v.files
            consensus.quality = consensus.quality + v.quality
            consensus.peaks = (consensus.peaks or []) + (v.peaks or [])
        merged.append(consensus)

    return sorted(merged, key=lambda v: pos_sort_key(v.pos))


# ---------------------------------------------------------------------------
# Main ETL: process()
# ---------------------------------------------------------------------------


def _maybe_add_hv2f_pos73(
    variants: list[dict[str, Any]],
    primers: PrimerTypeInfo,
    params: _CreateParams,
    file_intervals: dict[str, list[list[int]]],
    filename: str,
    sample: str | None,
) -> None:
    """Add an SNP variant at position 73 for HV2F traces when needed.

    When the primer is HV2F and the first variant genomic position is beyond
    73, create an SNP variant at position 73 and append the [73, 75] interval.
    """
    if not primers.is_hv2f or not variants or pos_base(variants[0]["pos"]) <= POS_73:
        return
    variant_73 = _create_variant(POS_73, "snp", params)
    if variant_73 is None or variant_73["ref"] == variant_73["seq"]:
        return
    logger.info("Sample {}: Add variant at position {}", sample, variant_73["pos"])
    variants.append(variant_73)
    if filename in file_intervals:
        file_intervals[filename].append([POS_73, 75])


def process(  # noqa: PLR0915
    sample_id: str,
    input_path: Path,
    *,
    ref_seq: str,
    config: _ProcessConfig | None = None,
) -> Sample:
    """Process a single Tracy decompose JSON file into a Sample.

    Flow:
    1. Load JSON data
    2. Normalize ref1pos (Tracy >= 0.7.8 compat)
    3. Compute trim bounds
    4. Parse raw variants from the alignment
    5. Detect variant conditions
    6. Apply transforms
    7. Update positions
    8. Filter variants
    9. Build Sample with flags

    Args:
        sample_id: Sample identifier.
        input_path: Path to Tracy decompose JSON output file.
        ref_seq: Reference genome sequence string.
        config: Processing configuration (uses defaults if None).

    Returns:
        Sample with parsed variants, intervals, source provenance, and flags.
    """
    if config is None:
        config = _ProcessConfig(
            quality_threshold=get_settings().tracy.quality_threshold,
            min_peak_value=get_settings().tracy.min_peak_value,
            pratio=get_settings().tracy.pratio,
            heteroplasmy_threshold=get_settings().tracy.heteroplasmy_threshold,
        )

    variant_regions = config.variant_regions
    if variant_regions is None:
        variant_regions = get_settings().regions.REGIONS

    quality_threshold = config.quality_threshold
    min_peak_value = config.min_peak_value
    pratio = config.pratio
    heteroplasmy_threshold = config.heteroplasmy_threshold
    batch_id = config.batch_id
    # ``config.sample`` is an optional free-text label for logging; fall back to
    # the real sample identifier so log messages never show "Sample None".
    sample = config.sample if config.sample is not None else sample_id

    # Load JSON data
    data = _load_tracy_json(input_path)

    # Normalize ref1pos from integer to position array (Tracy >= 0.7.8 compat)
    ref1pos_start = normalize_ref_positions(data)

    # Detect primer type from filename
    filename = input_path.stem
    primers = detect_primer_type(filename)

    # Compute trim bounds
    _alt_start, _alt_end, _ref_start, _ref_end, start_index, end_index = alignment_bounds(data)

    # Parse raw variants from the alignment
    call_params = _CallParams(
        data=data,
        deletions=0,
        filename=filename,
        heteroplasmy_threshold=heteroplasmy_threshold,
        sample=sample,
    )

    raw_variants: list[dict[str, Any]] = []
    deletions = 0
    i = start_index if start_index is not None else 0
    # Iterate the trimmed range inclusively: end_index is the last analysed
    # index, so use ``<=`` (fallback subtracts one for an exclusive-end shape).
    end = end_index if end_index is not None else len(data["alt1align"]) - 1

    while i <= end:
        call_params = call_params._replace(deletions=deletions)
        variant = _call_variant_at_index(i, call_params)
        if variant is not None:
            raw_variants.append(variant)
            if variant["ref"] != "-" and variant["seq"] == "-":
                deletions += 1
        i += 1

    # Update variant positions before transforms so condition detection and
    # transforms operate on genomic/decimal positions rather than raw indices.
    transform_params = _CreateParams(
        ref_seq=ref_seq,
        data=data,
        filename=filename,
        sample=sample,
        heteroplasmy_threshold=heteroplasmy_threshold,
    )
    position_params = _PositionParams(
        primers=primers,
        ref1pos=ref1pos_start,
        ref_align=data["ref1align"],
        filename=filename,
        sample=sample,
        heteroplasmy_threshold=heteroplasmy_threshold,
    )
    file_intervals = update_variant_positions(
        raw_variants,
        position_params,
        start_index,
        end_index,
    )

    # HV2F position-73 special handling.
    _maybe_add_hv2f_pos73(raw_variants, primers, transform_params, file_intervals, filename, sample)

    # Detect variant conditions after position update, so insertion positions
    # are decimal strings (e.g. "513.1") that conditions and remap rely on.
    conditions = detect_variant_conditions(raw_variants)

    # Apply transformations
    apply_all_transformations(raw_variants, conditions, primers, transform_params)

    # Filter variants
    filter_config = _FilterConfig(
        quality_threshold=quality_threshold,
        min_peak_value=min_peak_value,
        pratio=pratio,
        sample=sample,
        variant_regions=variant_regions,  # type: ignore[arg-type]
    )
    filtered_variants: list[dict[str, Any]] = []
    for variant in raw_variants:
        if _should_skip_variant(variant, filter_config):
            if variant.get("remove", False):
                reason = variant.get("reason", "marked for removal during conversion")
                logger.warning(
                    "Sample {} trace {}: Remove variant {} {}>{}, quality={} ({})",
                    sample,
                    filename,
                    variant["pos"],
                    variant["ref"],
                    variant["seq"],
                    variant.get("quality"),
                    reason,
                )
            continue

        # Remove index fields before building Variant objects
        variant.pop("index", None)
        variant.pop("peak_index", None)
        filtered_variants.append(variant)

    # Build Variant objects
    variant_list = _build_variant_objects(filtered_variants, filename)

    # Compute intervals
    intervals = _compute_intervals(file_intervals, variant_regions)  # type: ignore[arg-type]

    # Compute flags: per-variant and sample-level
    variant_flags_dict = flag_variants(variant_list)
    flagger = SampleFlagger(variant_list, intervals=intervals or None)
    _, analyze_reasons, _ = flagger.analyze()
    sample_flags_list = list(analyze_reasons)

    # Per-variant flags belong only in variant_flags; drop any sample-level
    # flag that duplicates a per-variant flag string.
    sample_flags_list = deduplicate_sample_flags(sample_flags_list, variant_flags_dict)

    # Build Sample — ETL is pure: only variants, intervals, and flags
    return Sample(
        sample_id=sample_id,
        variants=variant_list,
        source_tool=Tool.TRACY,
        intervals=intervals or None,
        batch_id=batch_id,
        sample_flags=sample_flags_list or [],
        variant_flags=variant_flags_dict,
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
