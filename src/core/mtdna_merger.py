"""Single-pass mtDNA profile merger.

Load multiple per-sample JSON files, auto-merge variants (with strict overlap consistency),
and write per-sample JSON to {sample_id}/{sample_id}.json.

Merged results use ``source_tool = "unified"`` (Tool.UNIFIED) and include
``region_sources`` in the ``information`` dict to track which tool produced each region's data.

``merged_from`` (sample-level tool list) has been replaced by ``region_sources`` (per-region provenance).

Conflicting alleles at the same position raise MergeValidationError (hard error).
Quality-based conflict resolution is not yet implemented; see ARCHITECTURE.md for
the deferred Regenerate Gate.

HV strings are always rebuilt from the merged Sample via
sample_to_dict() -> generate_sequence() so no single input file dominates the output.

Refactored to use Sample objects internally for consistency with the rest of
the core layer.  Output format matches sample_to_dict() (grouped variants + flags + information).
"""

import json
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from src.config import get_settings
from src.core.models import Sample, Tool, Variant
from src.core.region import REGION_TO_KEYS, validate_region_group_intervals
from src.core.sample import merge_intervals, sample_to_dict
from src.core.variants import (
    Position,
    normalize_position,
    pos_base,
    pos_sort_key,
)

MIN_SPAN_LENGTH = 2
MAX_WARNINGS_PREVIEW = 10


class MergeValidationError(ValueError):
    """Raised when merged variants fail validation checks or contain overlapping conflicts."""


class MergeResult:
    """Result of the auto-merge phase.

    Attributes:
        sample: Merged Sample object.
        warnings: List of warning strings from validation.

    """

    def __init__(self, sample: Sample, warnings: list[str] | None = None) -> None:
        self.sample = sample
        self.warnings = warnings or []


@dataclass
class MergedData:
    """Container for merged variant, interval, and flag data passed to build_merged_sample."""

    merged_variants: list[Variant]
    intervals: dict[str, list[list[int]]] | None
    all_sample_flags: list[str]
    all_variant_flags: dict[str, list[str]]
    merged_info: dict
    merged_from_tools: list[str]


def _variant_from_flat(flat: dict[Position, dict]) -> list[Variant]:
    """Convert flat variant dict {pos: {ref, seq, file, quality}} to Variant objects."""
    variants: list[Variant] = []
    for pos in sorted(flat.keys(), key=pos_sort_key):
        entry = flat[pos]
        raw_file = entry.get("file", [])
        files = [raw_file] if isinstance(raw_file, str) else list(raw_file)
        raw_quality = entry.get("quality", [])
        if isinstance(raw_quality, (int, float)):
            quality_list: list[int] = [int(raw_quality)]
        elif raw_quality:
            quality_list = [int(q) for q in raw_quality]
        else:
            quality_list = []
        variants.append(
            Variant(
                pos=pos,
                ref=entry["ref"],
                seq=entry["seq"],
                files=files,
                quality=quality_list,
                peaks=None,
            ),
        )
    return variants


def _flat_from_variant_list(variants: list[Variant]) -> dict[Position, dict]:
    """Convert a list of Variant objects to flat variant dict {pos: {ref, seq, file, quality}}."""
    flat: dict[Position, dict] = {}
    for v in variants:
        pos = normalize_position(v.pos)
        flat[pos] = {
            "ref": v.ref,
            "seq": v.seq,
            "file": sorted(set(v.files)) if v.files else [],
            "quality": list(v.quality) if v.quality else [],
        }
    return flat


class MtDnaMerger:
    """Single-pass merger for per-sample mtDNA JSON profiles.

    Enforces strict consistency:
    - All variants must lie inside their own profile's declared intervals.
    - Ref allele must exactly match rCRS (hard error).
    - Overlap regions must have identical alleles (explicit + implicit conflicts).

    Uses Sample objects internally and sample_to_dict() for output,
    ensuring the merged JSON matches the standard format used by the rest
    of the pipeline.
    """

    def __init__(self, ref_path: str | Path | None = None) -> None:
        self.ref_path: Path = Path(ref_path) if ref_path else get_settings().directories.ref / "rCRS.fasta"
        self._rcrs_ref: dict[int, str] | None = None

    # ----------------------------------------------------------------------
    # Private helpers
    # ----------------------------------------------------------------------

    def _load_profile(self, path: Path, _sample_id: str) -> dict:
        """Load one JSON file in unwrapped per-sample format.

        Per-sample JSONs use unwrapped format: top-level keys are profile
        fields (e.g. "no_snps", "HV1", "variants"), not a sample_id wrapper.
        """
        with path.open(encoding="utf-8") as fh:
            return json.load(fh)

    def _flat_variant_dict(self, profile: dict) -> dict[Position, dict]:
        """Flatten grouped variants {snps, insertions, deletions} into {pos: {ref, seq, file, quality}}."""
        flat: dict[Position, dict] = {}
        variants_raw = profile.get("variants", {})

        # Grouped format: {snps: [...], insertions: [...], deletions: [...]}
        if isinstance(variants_raw, dict) and "snps" in variants_raw:
            for typ in ("snps", "insertions", "deletions"):
                for var in variants_raw.get(typ, []):
                    pos = normalize_position(var["pos"])
                    raw_file = var.get("file", [])
                    file_list: list[str] = [raw_file] if isinstance(raw_file, str) else list(raw_file)
                    file_list = sorted(set(file_list))
                    raw_quality = var.get("quality", [])
                    quality_list: list[float] = (
                        [float(raw_quality)]
                        if isinstance(raw_quality, (int, float))
                        else list(raw_quality)
                        if isinstance(raw_quality, list)
                        else []
                    )
                    flat[pos] = {
                        "ref": var["ref"],
                        "seq": var["seq"],
                        "file": file_list,
                        "quality": quality_list,
                    }
            return flat

        return flat

    def merge_profiles(self, flat_variants_list: list[dict[Position, dict]]) -> dict[Position, dict]:
        """Merge flat variant dicts.
        - Same pos + same seq → union file[] and quality[]
        - Same pos + diff seq → raise MergeValidationError immediately
        """
        merged: dict[Position, dict] = {}

        for flat in flat_variants_list:
            for pos, var in flat.items():
                if pos not in merged:
                    merged[pos] = dict(var)
                elif merged[pos]["seq"] == var["seq"]:
                    merged[pos]["file"] = sorted(set(merged[pos]["file"] + var["file"]))
                    merged[pos]["quality"] = merged[pos].get("quality", []) + var.get("quality", [])
                else:
                    all_sources = sorted(set(merged[pos].get("file", []) + var.get("file", [])))
                    msg = (
                        f"Conflicting alleles at position {pos}: "
                        f"{merged[pos]['seq']} vs {var['seq']} "
                        f"(sources: {all_sources})"
                    )
                    raise MergeValidationError(
                        msg,
                    )
        return merged

    def _union_intervals(self, profiles: list[dict]) -> dict[str, list[list[int]]]:
        """Union-merge intervals per region using merge_intervals()."""
        regions = ["HV1", "HV2", "HV3"]
        intervals_by_region: dict[str, list[list[int]]] = {r: [] for r in regions}

        for profile in profiles:
            for region in regions:
                intervals_by_region[region].extend(profile.get("intervals", {}).get(region, []))

        return {region: merge_intervals(intervals_by_region[region]) for region in regions}

    def _load_rcrs_reference(self) -> dict[int, str]:
        """Load rCRS reference as {pos: base} dict using pure Python. Cached."""
        if self._rcrs_ref is not None:
            return self._rcrs_ref

        with self.ref_path.open(encoding="utf-8") as fh:
            seq = "".join(line.strip().upper() for line in fh if not line.startswith(">"))

        ref_dict: dict[int, str] = {}
        for i, base in enumerate(seq):
            if base in "ACGT":
                ref_dict[i + 1] = base

        self._rcrs_ref = ref_dict
        if not ref_dict:
            msg = f"rCRS reference file is empty or malformed: {self.ref_path}"
            raise MergeValidationError(msg)
        return self._rcrs_ref

    # ----------------------------------------------------------------------
    # Validation helpers
    # ----------------------------------------------------------------------

    def _validate_individual_profile(self, profile: dict, filename: str) -> None:
        """Validate a single profile: variants inside intervals, ref alleles match rCRS.

        Raises MergeValidationError on any hard error.
        """
        intervals = profile.get("intervals", {})
        bounds: list[tuple] = []
        for spans in intervals.values():
            bounds.extend(
                (span[0], span[1]) for span in spans if isinstance(span, (list, tuple)) and len(span) >= MIN_SPAN_LENGTH
            )

        rcrs = self._load_rcrs_reference()
        flat = self._flat_variant_dict(profile)

        for pos, var in flat.items():
            anchor = pos_base(pos)

            # Check ref allele matches rCRS
            if anchor in rcrs and var["ref"] != rcrs[anchor] and var["ref"] != "-":
                msg = (
                    f"ref allele '{var['ref']}' does not match rCRS "
                    f"base '{rcrs[anchor]}' at position {pos} in '{filename}'"
                )
                raise MergeValidationError(
                    msg,
                )

            # Check variant lies within declared intervals
            if bounds:
                in_any = any(start <= anchor <= end for start, end in bounds)
                if not in_any:
                    msg = (
                        f"Variant pos={pos} (anchor={anchor}) in '{filename}' "
                        f"lies outside the profile's declared intervals."
                    )
                    raise MergeValidationError(
                        msg,
                    )

    def _check_overlap_consistency(
        self,
        profiles: list[dict],
        merged_flat: dict[Position, dict],
        file_paths: list[str],
    ) -> None:
        """Check for implicit overlap conflicts: where one profile calls a variant
        and another covers that position but calls the rCRS reference base.

        This is a hard error because the merged result would be ambiguous.
        """
        rcrs = self._load_rcrs_reference()

        for _i, (profile, fp) in enumerate(zip(profiles, file_paths, strict=False)):
            intervals = profile.get("intervals", {})
            bounds: list[tuple] = []
            for spans in intervals.values():
                bounds.extend(
                    (span[0], span[1])
                    for span in spans
                    if isinstance(span, (list, tuple)) and len(span) >= MIN_SPAN_LENGTH
                )

            if not bounds:
                continue

            # Check positions covered by this profile's intervals
            for pos_str, var_data in merged_flat.items():
                anchor = pos_base(pos_str)

                # Is this position covered by profile i's intervals?
                in_interval = any(start <= anchor <= end for start, end in bounds)
                if not in_interval:
                    continue

                # Does profile i have this variant?
                flat_i = self._flat_variant_dict(profile)
                if pos_str in flat_i:
                    continue  # Explicit variant — already handled by merge

                # Profile covers this position but doesn't call a variant →
                # it implicitly calls rCRS reference base
                rcrs_base = rcrs.get(anchor)
                if rcrs_base and var_data["seq"] != rcrs_base:
                    msg = (
                        f"Overlap conflict at position {pos_str}: "
                        f"merged allele is '{var_data['seq']}', but "
                        f"{fp} covers this position and would call rCRS base '{rcrs_base}'"
                    )
                    raise MergeValidationError(
                        msg,
                    )

    def _validate_merged_profile(
        self,
        merged_flat: dict[Position, dict],
        merged_profile: dict,
    ) -> list[str]:
        """Light validation on merged result. Returns list of warning strings."""
        warnings: list[str] = []

        intervals = merged_profile.get("intervals", {})
        bounds: list[tuple] = []
        for spans in intervals.values():
            bounds.extend(
                (span[0], span[1]) for span in spans if isinstance(span, (list, tuple)) and len(span) >= MIN_SPAN_LENGTH
            )

        for pos, var in merged_flat.items():
            anchor = pos_base(pos)
            seq = var["seq"]

            # Check for non-standard bases (including N)
            is_insertion = isinstance(pos, str) and "." in str(pos)
            non_standard_base = (seq.upper() not in ("A", "C", "G", "T", "-") and not is_insertion) or (
                is_insertion and len(seq) == 1 and seq.upper() not in ("A", "C", "G", "T")
            )
            if non_standard_base:
                warnings.append(f"Variant at pos={pos} has non-standard base '{seq}'")

            # Check variant lies within intervals
            if bounds:
                in_any = any(start <= anchor <= end for start, end in bounds)
                if not in_any:
                    warnings.append(f"Variant pos={pos} (anchor={anchor}) not in any declared interval.")

        return warnings

    # ----------------------------------------------------------------------
    # Public API
    # ----------------------------------------------------------------------

    def _load_and_flatten_profiles(
        self,
        json_files: list[str | Path],
        sample_id: str,
    ) -> tuple[list[dict], list[dict[Position, dict]], list[str]]:
        """Load profiles, validate individually, and flatten variant dicts."""
        profiles: list[dict] = []
        flat_variants_list: list[dict[Position, dict]] = []
        file_paths = [str(Path(fp)) for fp in json_files]

        for fp in json_files:
            path = Path(fp)
            profile = self._load_profile(path, sample_id)
            profiles.append(profile)
            self._validate_individual_profile(profile, path.name)
            flat_variants_list.append(self._flat_variant_dict(profile))

        return profiles, flat_variants_list, file_paths

    @staticmethod
    def _collect_flags_and_info(
        profiles: list[dict],
    ) -> tuple[list[str], dict[str, list[str]], dict, list[str]]:
        """Collect deduplicated flags, information, and source tools from profiles."""
        all_sample_flags: list[str] = []
        all_variant_flags: dict[str, list[str]] = {}
        merged_info: dict = {}
        merged_from_tools: list[str] = []

        for profile in profiles:
            info = profile.get("information", {})
            if isinstance(info, dict):
                merged_info.update(info)
                st = info.get("source_tool")
                if st and st not in merged_from_tools:
                    merged_from_tools.append(st)

            for flag in profile.get("sample_flags", []):
                if flag not in all_sample_flags:
                    all_sample_flags.append(flag)
            for key, reasons in profile.get("variant_flags", {}).items():
                if key not in all_variant_flags:
                    all_variant_flags[key] = reasons

        return all_sample_flags, all_variant_flags, merged_info, merged_from_tools

    @staticmethod
    def build_merged_sample(
        sample_id: str,
        merged_data: MergedData,
        batch_id: str | None,
    ) -> Sample:
        """Construct the merged Sample object with unified metadata."""
        merged_data.merged_info["source_tool"] = Tool.UNIFIED.value
        merged_data.merged_info["region_sources"] = dict.fromkeys(("HV1", "HV2", "HV3"), merged_data.merged_from_tools)

        return Sample(
            sample_id=sample_id,
            variants=merged_data.merged_variants,
            source_tool=Tool.UNIFIED,
            intervals=merged_data.intervals,
            sample_flags=merged_data.all_sample_flags,
            variant_flags=merged_data.all_variant_flags,
            information=merged_data.merged_info or None,
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

    def merge(
        self,
        json_files: list[str | Path],
        sample_id: str,
        output_dir: str | Path = ".",
    ) -> MergeResult:
        """Merge multiple per-sample mtDNA JSON profiles into one Sample.

        Reads per-sample JSON files (grouped variant format),
        validates consistency, merges variants, and returns a MergeResult
        containing the merged Sample and any warnings.

        Writes per-sample JSON to {output_dir}/{sample_id}/{sample_id}.json
        (unwrapped format, same as other pipeline tools).

        Raises MergeValidationError for any inconsistency.

        Returns:
            MergeResult with merged Sample and warnings list.

        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        profiles, flat_variants_list, file_paths = self._load_and_flatten_profiles(json_files, sample_id)

        merged_flat = self.merge_profiles(flat_variants_list)
        self._check_overlap_consistency(profiles, merged_flat, file_paths)
        intervals = self._union_intervals(profiles)
        merged_variants = _variant_from_flat(merged_flat)

        all_sample_flags, all_variant_flags, merged_info, merged_from_tools = self._collect_flags_and_info(profiles)

        batch_id = merged_info.pop("batch_id", None) or merged_info.pop("batch", None)

        merged_data = MergedData(
            merged_variants=merged_variants,
            intervals=intervals or None,
            all_sample_flags=all_sample_flags,
            all_variant_flags=all_variant_flags,
            merged_info=merged_info,
            merged_from_tools=merged_from_tools,
        )
        sample = self.build_merged_sample(
            sample_id=sample_id,
            merged_data=merged_data,
            batch_id=batch_id,
        )

        merged_profile_dict = sample_to_dict(sample, str(self.ref_path))
        all_warnings = self._validate_merged_profile(merged_flat, merged_profile_dict)
        for w in all_warnings:
            logger.warning(f"[merge] {w}")

        if all_warnings:
            warning_preview = "\n".join(f"- {w}" for w in all_warnings[:MAX_WARNINGS_PREVIEW])
            if len(all_warnings) > MAX_WARNINGS_PREVIEW:
                warning_preview += f"\n- ... and {len(all_warnings) - MAX_WARNINGS_PREVIEW} more warning(s)"
            msg = f"Merge failed validation with {len(all_warnings)} warning(s):\n{warning_preview}"
            raise MergeValidationError(msg)

        sample_dir = output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_path = sample_dir / f"{sample_id}.json"
        with sample_path.open("w", encoding="utf-8") as fh:
            json.dump(merged_profile_dict, fh, indent=2)

        return MergeResult(sample=sample, warnings=all_warnings)


# ----------------------------------------------------------------------
# Convenience: merge from Sample objects directly
# ----------------------------------------------------------------------


def _collect_sample_flags_and_info(
    samples: list[Sample],
) -> tuple[list[str], dict[str, list[str]], dict, str | None]:
    """Collect deduplicated flags, information, and batch_id from Sample objects."""
    all_sample_flags: list[str] = []
    all_variant_flags: dict[str, list[str]] = {}
    merged_info: dict = {}
    batch_id = None

    for s in samples:
        for flag in s.sample_flags:
            if flag not in all_sample_flags:
                all_sample_flags.append(flag)
        for key, reasons in s.variant_flags.items():
            if key not in all_variant_flags:
                all_variant_flags[key] = reasons
        if s.information:
            merged_info.update(s.information)
        if s.batch_id and not batch_id:
            batch_id = s.batch_id

    return all_sample_flags, all_variant_flags, merged_info, batch_id


def _collect_source_tools(samples: list[Sample]) -> list[str]:
    """Collect deduplicated source tool names from Sample objects."""
    merged_from_tools: list[str] = []
    for s in samples:
        tool_name = s.source_tool.value
        if tool_name not in merged_from_tools:
            merged_from_tools.append(tool_name)
    return merged_from_tools


def _write_sample_json(
    sample: Sample,
    sample_id: str,
    output_dir: Path,
    ref_path: str | Path,
) -> Path:
    """Write a Sample to per-sample JSON and return the path."""
    sample_dict = sample_to_dict(sample, str(ref_path))
    sample_dir = output_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    sample_path = sample_dir / f"{sample_id}.json"
    with sample_path.open("w", encoding="utf-8") as fh:
        json.dump(sample_dict, fh, indent=2)
    return sample_path


def merge_samples(
    samples: list[Sample],
    ref_path: str | Path | None = None,
    output_dir: str | Path = ".",
) -> MergeResult:
    """Merge a list of Sample objects into a single merged Sample.

    Higher-level API that accepts Sample objects directly instead of JSON files.
    Validates consistency and merges variants from all input samples.

    Args:
        samples: List of Sample objects to merge (must all have the same sample_id).
        ref_path: Path to rCRS reference FASTA.
        output_dir: Directory for output files.

    Returns:
        MergeResult with merged Sample and warnings list.

    Raises:
        MergeValidationError: If variants conflict across samples.
        ValueError: If samples have different sample_ids.

    """
    if not samples:
        msg = "Cannot merge empty sample list."
        raise ValueError(msg)

    sample_ids = {s.sample_id for s in samples}
    if len(sample_ids) > 1:
        msg = f"All samples must have the same sample_id, got: {sample_ids}"
        raise ValueError(msg)

    sample_id = samples[0].sample_id
    merger = MtDnaMerger(ref_path=ref_path)

    flat_variants_list = [_flat_from_variant_list(s.variants) for s in samples]
    merged_flat = merger.merge_profiles(flat_variants_list)

    all_intervals: dict[str, list[list[int]]] = {"HV1": [], "HV2": [], "HV3": []}
    for s in samples:
        if s.intervals:
            for region in ("HV1", "HV2", "HV3"):
                all_intervals[region].extend(s.intervals.get(region, []))
    merged_intervals = {region: merge_intervals(all_intervals[region]) for region in ("HV1", "HV2", "HV3")}

    all_sample_flags, all_variant_flags, merged_info, batch_id = _collect_sample_flags_and_info(samples)
    merged_from_tools = _collect_source_tools(samples)
    merged_variants = _variant_from_flat(merged_flat)

    merged_data = MergedData(
        merged_variants=merged_variants,
        intervals=merged_intervals,
        all_sample_flags=all_sample_flags,
        all_variant_flags=all_variant_flags,
        merged_info=merged_info,
        merged_from_tools=merged_from_tools,
    )
    result = MtDnaMerger.build_merged_sample(
        sample_id=sample_id,
        merged_data=merged_data,
        batch_id=batch_id,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_sample_json(result, sample_id, output_dir, merger.ref_path)

    return MergeResult(sample=result, warnings=[])


TOOL_PRIORITY: dict[str, int] = {
    "sequencher": 1,
    "mutation_surveyor": 2,
    "tracy": 3,
    "blastn": 4,
}


@dataclass
class RegionGroupState:
    """Mutable accumulator for merging region group data across samples."""

    all_intervals: dict[str, list[list[int]]]
    all_sample_flags: list[str]
    all_variant_flags: dict[str, list[str]]
    merged_info: dict
    batch_id: str | None = None

    def merge_sample(self, sample: Sample, region_keys: list[str]) -> None:
        """Merge a single sample's data into this state."""
        if sample.intervals:
            for rk in region_keys:
                if rk in sample.intervals:
                    self.all_intervals[rk].extend(sample.intervals[rk])

        for flag in sample.sample_flags:
            if flag not in self.all_sample_flags:
                self.all_sample_flags.append(flag)
        for key, reasons in sample.variant_flags.items():
            if key not in self.all_variant_flags:
                self.all_variant_flags[key] = reasons

        if sample.information:
            self.merged_info.update(sample.information)
        if sample.batch_id and not self.batch_id:
            self.batch_id = sample.batch_id


def _process_region_group(
    group_name: str,
    samples: list[Sample],
    state: RegionGroupState,
) -> tuple[list[Variant], str, str | None, dict[str, str], list[str]]:
    """Process one region group: collect variants, intervals, flags, info.

    Returns (winning_variants, winning_tool, updated_batch_id, region_sources, qc_warnings).
    """
    region_keys = REGION_TO_KEYS.get(group_name, [group_name])

    tools_in_group = sorted(
        {s.source_tool.value for s in samples},
        key=lambda t: TOOL_PRIORITY.get(t, 99),
    )
    winning_tool = tools_in_group[0]
    winning_samples = [s for s in samples if s.source_tool.value == winning_tool]

    winning_variants: list[Variant] = []
    for s in winning_samples:
        winning_variants.extend(s.variants)

    for s in samples:
        state.merge_sample(s, region_keys)

    region_sources = dict.fromkeys(region_keys, winning_tool)

    qc_warnings: list[str] = []
    for s in samples:
        if s.intervals:
            qc_warnings.extend(validate_region_group_intervals(s.intervals, group_name))

    return winning_variants, winning_tool, state.batch_id, region_sources, qc_warnings


def _check_variant_conflicts(all_variants: list[Variant]) -> dict[Position, Variant]:
    """Build pos_map from variants, checking for allele conflicts.

    Raises MergeValidationError if conflicting alleles exist.
    """
    pos_map: dict[Position, Variant] = {}
    for v in all_variants:
        key = normalize_position(v.pos)
        if key in pos_map:
            existing = pos_map[key]
            if existing.ref != v.ref or existing.seq != v.seq:
                msg = (
                    f"Conflicting alleles at position {key}: "
                    f"existing={existing.ref}>{existing.seq} vs new={v.ref}>{v.seq}"
                )
                raise MergeValidationError(msg)
            existing.files = sorted(set(existing.files + v.files))
            existing.quality = list(dict.fromkeys(list(existing.quality + v.quality)))
        else:
            pos_map[key] = v
    return pos_map


def merge_by_regions(
    region_samples: dict[str, list[Sample]],
    sample_id: str,
    ref_path: str | Path | None = None,
    output_dir: str | Path = ".",
) -> MergeResult:
    """Merge per-region Samples into a unified Sample with region-level provenance.

    This is the region-aware merge that replaces sample-level merged_from with
    per-region region_sources tracking. Each region group maps to the tool(s)
    that produced data for that region.

    Per-region intermediate JSONs are temporary — the output is a single
    unified Sample JSON per LID in the standard format.

    Args:
        region_samples: Dict mapping region group name → list of Sample objects
            for that region. E.g. {"HV1": [ms_hv1_sample], "HV2-3": [seq_hv23_sample]}.
        sample_id: The unified sample ID for all region Samples.
        ref_path: Path to rCRS reference FASTA.
        output_dir: Directory for output files.

    Returns:
        MergeResult with unified Sample and warnings list.

    Raises:
        MergeValidationError: If conflicting alleles exist in overlapping regions.
        ValueError: If no region samples are provided.

    """
    if not region_samples:
        msg = "Cannot merge empty region_samples dict."
        raise ValueError(msg)

    if ref_path is None:
        ref_path = get_settings().directories.ref / "rCRS.fasta"

    merger = MtDnaMerger(ref_path=ref_path)

    all_variants: list[Variant] = []
    all_intervals: dict[str, list[list[int]]] = {"HV1": [], "HV2": [], "HV3": []}
    all_sample_flags: list[str] = []
    all_variant_flags: dict[str, list[str]] = {}
    merged_info: dict = {}
    batch_id: str | None = None
    region_sources: dict[str, str] = {}
    all_warnings: list[str] = []

    for group_name, samples in region_samples.items():
        region_state = RegionGroupState(
            all_intervals=all_intervals,
            all_sample_flags=all_sample_flags,
            all_variant_flags=all_variant_flags,
            merged_info=merged_info,
            batch_id=batch_id,
        )
        winning_variants, _winning_tool, batch_id, group_sources, qc_warnings = _process_region_group(
            group_name,
            samples,
            region_state,
        )
        all_variants.extend(winning_variants)
        region_sources.update(group_sources)
        all_warnings.extend(qc_warnings)

    pos_map = _check_variant_conflicts(all_variants)

    merged_intervals = {region: merge_intervals(all_intervals[region]) for region in ("HV1", "HV2", "HV3")}
    merged_intervals = {k: v for k, v in merged_intervals.items() if v}

    merged_info["source_tool"] = Tool.UNIFIED.value
    merged_info["region_sources"] = region_sources
    merged_info.pop("merged_from", None)

    result = Sample(
        sample_id=sample_id,
        variants=sorted(pos_map.values(), key=lambda v: pos_sort_key(v.pos)),
        source_tool=Tool.UNIFIED,
        intervals=merged_intervals or None,
        sample_flags=all_sample_flags,
        variant_flags=all_variant_flags,
        information=merged_info or None,
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

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_sample_json(result, sample_id, output_dir, merger.ref_path)

    return MergeResult(sample=result, warnings=all_warnings)


# ----------------------------------------------------------------------
# CLI entry point
# ----------------------------------------------------------------------


def main() -> None:
    parser = ArgumentParser(description="Merge 2+ mtDNA JSON profiles (single-pass).")
    parser.add_argument(
        "--files",
        type=Path,
        nargs="+",
        required=True,
        help="Input sample_mtDNA.json files",
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        required=True,
        help="Authoritative sample ID; all input files must match",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(),
        help="Directory for output files",
    )
    parser.add_argument(
        "--ref-path",
        type=Path,
        default=None,
        help="Path to rCRS.fasta (overrides default from settings)",
    )
    args = parser.parse_args()

    merger = MtDnaMerger(ref_path=args.ref_path)
    try:
        merger.merge(args.files, args.sample_id, args.output_dir)
    except MergeValidationError as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc
    except Exception as exc:
        logger.error(f"Merge failed: {exc}")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
