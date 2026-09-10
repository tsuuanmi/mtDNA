"""mtDNA variant and sample flagging.

Produces per-variant flags (``flag_variants``) and aggregate sample-level
flag reasons + severity levels (``SampleFlagger.analyze``/``classify``)
from a normalized variant set. Most flags are *presence* flags tied to a
specific variant; the ``"No 315.1 variant"`` flag is a special *absence*
check (a variant missing from the list) and is derived from the variant
set + coverage intervals via :func:`recompute_no_315_1_flag`.
"""

import re
from dataclasses import dataclass, field
from typing import Any, ClassVar

from src.core.models import Variant
from src.core.variants import (
    get_hv_region_for_position,
    is_special_position,
    normalize_position,
    pos_base,
    pos_sort_key,
    variant_to_dict,
)

# Configuration constants
IUPAC_AMBIGUITY_CODES = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
FLAGGING_REGION_START = 16180
FLAGGING_REGION_END = 16193
AC_REPEAT_REGION_START = 513
AC_REPEAT_REGION_END = 525
INSERTION_C_EXCEPTION_POSITIONS = {309, 573}
# Well-known mtDNA positions used in flagging logic
POS_249 = 249
POS_309 = 309
POS_459 = 459
POS_513 = 513
POS_523 = 523
POS_524 = 524
POS_16192 = 16192
POS_16258 = 16258
POS_16263 = 16263
POS_16293 = 16293

# Per-variant flag pattern data (shared between _apply_specific_variant_flags and _check_specific_patterns)
_POSITION_VARIANT_FLAGS: list[tuple[str, str, str]] = [
    ("309", "T", "Has 309T variant"),
    ("309", "Y", "Has 309Y variant"),
]
_POSITION_FLAGS: list[tuple[str, str]] = [
    ("310", "Has 310 variant"),
    ("460", "Has 460 variant"),
]
_TRANSVERSION_CHECKS: list[tuple[int, str, str, str]] = [
    (POS_16293, "A", "C", "Has 16293A-C variant"),
    (POS_16258, "A", "C", "Has 16258A-C variant"),
    (POS_16263, "T", "C", "Has 16263T-C variant"),
]


def _has_consecutive_run(values: list[int], count: int) -> bool:
    """Return True if `values` contains `count` consecutive integers.

    Duplicates are collapsed, since a consecutive run requires distinct
    adjacent integers. Used for both deletion base positions (e.g.
    100, 101, 102, 103, 104) and insertion indices at one base (e.g.
    1, 2, 3, 4, 5 for 100.1, 100.2, ..., 100.5).

    """
    if len(values) < count:
        return False
    unique = sorted(set(values))
    run = 1
    for i in range(1, len(unique)):
        if unique[i] == unique[i - 1] + 1:
            run += 1
            if run >= count:
                return True
        else:
            run = 1
    return run >= count


@dataclass
class VariantFlag:
    """Per-variant flag result.

    Associates a specific variant (identified by position, ref, alt) with
    the list of flag reasons that apply to it.

    Attributes:
        pos: Variant position (normalized).
        ref: Reference base.
        seq: Alternate sequence base.
        flags: List of flag reason strings for this variant.

    """

    pos: int | str
    ref: str
    seq: str
    flags: list[str] = field(default_factory=list)

    def key(self) -> str:
        """Stable identity key for this variant (pos|ref|alt)."""
        return f"{normalize_position(self.pos)}|{self.ref}|{self.seq}"


class VariantAnalyzer:
    """Analyze variants for flagging patterns.

    Accepts Variant objects, dicts, or a mix. Internally normalizes
    to simple dicts with 'pos', 'ref', 'seq' keys for flagging logic.
    Dict input is deprecated; use Variant objects instead.
    """

    def __init__(self, variants: Variant | dict | list[dict[str, Any]]) -> None:
        """Initialize analyzer with variants.

        Args:
            variants: List of Variant objects, list of dicts with
                'pos', 'ref', 'seq' keys, or a single Variant/dict.

        """
        if isinstance(variants, list) and variants and isinstance(variants[0], dict) and "pos" in variants[0]:
            # Already in simple dict format — use directly
            self.variants: list[dict[str, Any]] = variants
        else:
            # Normalize to simple dicts via variant_to_dict. A single
            # Variant/dict is wrapped in a list so it is not iterated field-by-field.
            variant_list = variants if isinstance(variants, list) else [variants]
            self.variants = [v for v in (variant_to_dict(v) for v in variant_list) if v is not None]

    def has_consecutive_indels(self, count: int = 5) -> bool:
        """Check for N+ truly consecutive insertion/deletion variants.

        "Truly consecutive" differs by indel type:

        - Deletions: consecutive integer base positions, e.g. deletions at
          100, 101, 102, 103, 104 (100DEL 101DEL ... 104DEL).
        - Insertions: consecutive insertion indices at the same base
          position, e.g. 100.1C 100.2C 100.3C 100.4C 100.5C.

        Indels that are merely adjacent in the sorted variant list but do not
        meet these position/index criteria do not trigger this flag.

        Args:
            count: Minimum run length of consecutive same-type indels.

        """
        if count <= 0 or len(self.variants) < count:
            return False

        del_positions: list[int] = []
        ins_by_base: dict[int, list[int]] = {}
        for v in self.variants:
            ref, seq = v.get("ref", ""), v.get("seq", "")
            pos: str | int = v.get("pos", 0)
            is_deletion = seq == "-"
            is_insertion = ref == "-" or (isinstance(pos, str) and "." in str(pos))
            if is_deletion:
                del_positions.append(pos_base(pos))
            elif is_insertion:
                base, idx = pos_sort_key(pos)
                ins_by_base.setdefault(base, []).append(idx)

        if _has_consecutive_run(del_positions, count):
            return True
        return any(_has_consecutive_run(indices, count) for indices in ins_by_base.values())

    def has_lowercase_variants(self) -> bool:
        """Check for lowercase letters (low-confidence calls)."""
        return any(v.get("seq") and any(c.islower() for c in str(v.get("seq", ""))) for v in self.variants)

    def check_combined_conditions(self) -> list[str]:
        """Check for complex combined conditions."""
        reasons = []

        # 459 deletion + 16192 issue
        has_459_del = any(pos_base(v.get("pos", 0)) == POS_459 and v.get("seq", "") == "-" for v in self.variants)
        if has_459_del:
            for v in self.variants:
                pos: str | int = v.get("pos", 0)
                seq: str = v.get("seq", "")
                if pos_base(pos) == POS_16192 and (any(c in str(seq) for c in IUPAC_AMBIGUITY_CODES) or seq == "T"):
                    reasons.append("Complex: 459 deletion with 16192 issue")
                    break
        return reasons

    def has_variant_at(self, position: int | str, seq: str | None = None) -> bool:
        """Check if variant exists at position (optionally with specific seq)."""
        for v in self.variants:
            if normalize_position(v.get("pos", 0)) == normalize_position(position) and (
                seq is None or v.get("seq", "") == seq
            ):
                return True
        return False

    def has_523_524_ac_deletion(self) -> bool:
        """Check for paired AC deletion at 523-524."""
        has_523 = any(
            pos_base(v.get("pos", 0)) == POS_523 and v.get("ref") == "A" and v.get("seq") == "-" for v in self.variants
        )
        has_524 = any(
            pos_base(v.get("pos", 0)) == POS_524 and v.get("ref") == "C" and v.get("seq") == "-" for v in self.variants
        )
        return has_523 and has_524


class SampleFlagger:
    """Flag and classify samples based on variant patterns."""

    # Flag level constants
    LEVEL_NONE = 0
    LEVEL_LOW = 1
    LEVEL_MEDIUM = 2
    LEVEL_HIGH = 3
    LEVEL_CRITICAL = 4

    # Flag reason pattern to severity level mapping
    FLAG_REASON_MAP: ClassVar[dict[str, int]] = {
        # Level 4 (Critical)
        "Low depth": 4,
        "Failed position": 4,
        "Possible wrong consensus": 4,
        # Level 1 (Low)
        "Deletion at 523": 1,
        "Deletion at 524": 1,
        "16180-16193 region": 1,
        "309T": 1,
        # Level 2 (Medium)
        # AC repeat region SNP flags (positions 513-525)
        "Has 513": 2,
        "Has 514": 2,
        "Has 515": 2,
        "Has 516": 2,
        "Has 517": 2,
        "Has 518": 2,
        "Has 519": 2,
        "Has 520": 2,
        "Has 521": 2,
        "Has 522": 2,
        "Has 523": 2,
        "Has 524": 2,
        "Has 525": 2,
        " at 16296.1": 2,
        " at 290.1": 2,
        " at 290.2": 2,
        " at 16182.1": 2,
        " at 16193.1": 2,
        " at 16261.1": 2,
        " at 16192.1": 2,
        "309Y": 2,
        "460 variant": 2,
        # Level 3 (High) - specific first
        "310 variant": 3,
        "Deletion at 315": 3,
        "No 315.1": 3,
        "Deletion at 459": 3,
        "Deletion at 513": 3,
        "Deletion at 514": 3,
        "Deletion at 16193": 3,
        "Deletion at 16182": 3,
        "16293A-C": 3,
        "16258A-C": 3,
        "16263T-C": 3,
        " at 16188.1": 3,
        "Complex: 459 deletion": 3,
        "Consecutive indels": 3,
        "Lowercase variant": 3,
        # Level 3 (High) - general (MUST be last)
        "Heteroplasmy": 3,
        "Insertion": 3,
        "Deletion": 3,
    }

    # Classification patterns
    QUALITY_KEYWORDS: ClassVar[list[str]] = [
        "Low depth",
        "Failed position",
        "Possible wrong consensus",
    ]
    POLYC_PATTERNS: ClassVar[list[str]] = [
        r"Deletion at 523",
        r"Deletion at 524",
        r"16180-16193 region",
        r"Insertion \S+ at 16193\.1",
        r"Deletion.*16193",
    ]
    AC_REPEAT_SNP_PATTERN = r"Has \d+[ACGT]-[A-Z] variant"
    NOMENCLATURE_PATTERNS: ClassVar[list[str]] = [
        r"Insertion \S+ at 290\.1.*Insertion \S+ at 290\.2",
        r"Deletion at 513.*Deletion at 514",
        r"Insertion \S+ at 16182\.1.*Deletion at 16193",
        r"Deletion at 16182.*Insertion \S+ at 16193\.1",
        r"Insertion.*16261\.1",
    ]
    COMPLEX_PATTERNS: ClassVar[list[str]] = [
        r"Has 310 variant",
        r"Deletion at 315",
        r"No 315\.1 variant",
        r"Insertion \S+ at 16188\.1.*Insertion \S+ at 16193\.1",
        r"Has 16258A-C variant",
        r"Has 16263T-C variant",
    ]

    def __init__(
        self,
        variants: Variant | dict | list,
        intervals: dict[str, list[list[int]]] | None = None,
    ) -> None:
        """Initialize flagger with variants.

        Args:
            variants: List of dicts with 'pos', 'ref', 'seq' keys, or
                list of Variant objects, or a single Variant/dict.
            intervals: Optional coverage intervals keyed by region
                (e.g. ``{"HV1": [[16024, 16365]], "HV2": [[73, 340]]}``).
                When provided, absence checks such as "No 315.1 variant"
                only fire when the region containing the relevant position
                was analyzed (non-empty intervals for that region), so
                HV1-only / unsequenced-region samples are not flagged. When
                None, full coverage is assumed (legacy behavior), preserving
                the original flag output.

        """
        self.analyzer = VariantAnalyzer(variants)
        if isinstance(variants, list) and variants and isinstance(variants[0], dict) and "pos" in variants[0]:
            self.variants = variants or []
        else:
            variant_list = variants if isinstance(variants, list) else [variants]
            self.variants = [v for v in (variant_to_dict(v) for v in variant_list) if v is not None]
        self.intervals = intervals

    @classmethod
    def level_from_reasons(cls, flag_reasons: list[str]) -> int:
        """Compute maximum flag level from flag reasons."""
        if not flag_reasons:
            return cls.LEVEL_NONE

        max_level = cls.LEVEL_NONE
        for reason in flag_reasons:
            for pattern, level in cls.FLAG_REASON_MAP.items():
                if pattern in reason:
                    max_level = max(max_level, level)
                    break
        return max_level

    def _is_known_deletion_position(self, pos: str | int) -> bool:
        """Return True if deletion at this position is known and should be skipped (249 or 523-524 AC deletion)."""
        return pos_base(pos) == POS_249 or (pos_base(pos) in (523, 524) and self.analyzer.has_523_524_ac_deletion())

    def _is_309_deletion(self, pos: str | int, seq: str) -> bool:
        """Return True if this is a deletion at the base position 309 (309DEL).

        Targets the integer base position 309 only, not insertion sub-positions
        such as ``309.1``. The 309 deletion is a real variant event (e.g. the
        Tracy ``309C>T`` + ``310T>C`` conversion) that must be flagged even
        though 309 is otherwise a special polyC position whose non-insertion
        variants are normally skipped.
        """
        return pos_base(pos) == POS_309 and seq == "-" and not (isinstance(pos, str) and "." in pos)

    def _check_deletion(self, pos: str | int, seq: str, reasons: list[str]) -> bool:
        """Check for deletion flags. Returns True if deletion should be skipped (249/523-524)."""
        if seq != "-":
            return False
        if self._is_known_deletion_position(pos):
            return True
        reasons.append(f"Deletion at {pos}")
        return False

    def _check_ac_repeat_snp(self, pos: str | int, ref: str, seq: str) -> str | None:
        """Check for AC repeat region SNP. Returns flag string or None."""
        base_pos = pos_base(pos)
        if not (AC_REPEAT_REGION_START <= base_pos <= AC_REPEAT_REGION_END):
            return None
        is_snp = ref != "-" and seq != "-" and len(ref) == 1 and len(seq) == 1
        if is_snp and not (base_pos == POS_513 and ref == "A" and self.analyzer.has_523_524_ac_deletion()):
            return f"Has {pos}{ref}-{seq} variant"
        return None

    def _analyze_variant_reasons(self, v: dict) -> tuple[list[str], list[str]]:
        """Collect flag reasons and region positions for a single variant.

        Returns (reasons, region_positions).
        """
        reasons: list[str] = []
        region_positions: list[str] = []
        pos: str | int = v.get("pos", 0)
        ref: str = v.get("ref", "")
        seq: str = v.get("seq", "")

        if is_special_position(pos, ref, seq) and ref != "-" and not self._is_309_deletion(pos, seq):
            return reasons, region_positions

        if ref == "-":
            # Only the expected 315.1C insertion (normal rCRS reference
            # insertion) is suppressed; extra insertions at 315.2/315.3 and
            # other 315.x are flagged.
            if not (normalize_position(pos) == "315.1" and seq.upper() == "C") and not (
                seq.upper() == "C" and pos_base(pos) in INSERTION_C_EXCEPTION_POSITIONS
            ):
                reasons.append(f"Insertion {seq} at {pos}")
        elif self._check_deletion(pos, seq, reasons):
            return reasons, region_positions

        if any(c in seq for c in IUPAC_AMBIGUITY_CODES):
            reasons.append(f"Heteroplasmy at {pos}")

        if FLAGGING_REGION_START <= pos_base(pos) <= FLAGGING_REGION_END:
            region_positions.append(str(pos))

        snp_flag = self._check_ac_repeat_snp(pos, ref, seq)
        if snp_flag:
            reasons.append(snp_flag)

        return reasons, region_positions

    def analyze(self) -> tuple[bool, list[str], int]:
        """Analyze sample and return flagging results.

        Returns:
            Tuple of (should_flag, flag_reasons, flag_level).

        """
        if not self.variants:
            return False, [], self.LEVEL_NONE

        reasons: list[str] = []
        region_positions: list[str] = []

        for v in self.variants:
            var_reasons, var_region_positions = self._analyze_variant_reasons(v)
            reasons.extend(var_reasons)
            region_positions.extend(var_region_positions)

        if region_positions:
            unique = sorted(set(region_positions))
            reasons.append(f"16180-16193 region ({', '.join(unique)})")

        reasons.extend(self._check_specific_patterns())
        reasons.extend(self.analyzer.check_combined_conditions())

        if self.analyzer.has_consecutive_indels(5):
            reasons.append("Consecutive indels (5+)")
        if self.analyzer.has_lowercase_variants():
            reasons.append("Lowercase variant detected")

        flag_level = self.level_from_reasons(reasons)
        return bool(reasons), reasons, flag_level

    def _collect_variant_flags(self, v: dict) -> list[str] | None:
        """Collect per-variant flags for a single variant dict.

        Returns list of flag strings, or None if variant should be skipped
        (special non-insertion position).
        """
        pos: str | int = v.get("pos", 0)
        ref: str = v.get("ref", "")
        seq: str = v.get("seq", "")

        # Special positions: skip special-position deletions (and 309DEL is kept).
        # Insertions (ref == "-") are NOT skipped here — they flow to the
        # insertion check where only 309.1C/573.1C are suppressed; non-C
        # insertions at 309/573 (and 455.1/463.1) are flagged.
        if is_special_position(pos, ref, seq) and ref != "-" and not self._is_309_deletion(pos, seq):
            return None

        flags: list[str] = self._build_variant_flags(pos, ref, seq)
        return flags

    def _build_variant_flags(self, pos: str | int, ref: str, seq: str) -> list[str]:
        """Build flag list for a variant given its position, ref, and seq."""
        flags: list[str] = []

        # Only the expected 315.1C insertion is suppressed; other 315.x
        # insertions (315.2, 315.3, ...) are flagged.
        if ref == "-":
            if not (normalize_position(pos) == "315.1" and seq.upper() == "C") and not (
                seq.upper() == "C" and pos_base(pos) in INSERTION_C_EXCEPTION_POSITIONS
            ):
                flags.append(f"Insertion {seq} at {pos}")

        # Deletions (skip known-irrelevant deletions at 249/523-524)
        elif seq == "-" and not self._is_known_deletion_position(pos):
            flags.append(f"Deletion at {pos}")

        # Heteroplasmy
        if any(c in seq for c in IUPAC_AMBIGUITY_CODES):
            flags.append(f"Heteroplasmy at {pos}")

        # 16180-16193 region
        if FLAGGING_REGION_START <= pos_base(pos) <= FLAGGING_REGION_END:
            flags.append("16180-16193 region")

        # AC repeat region SNP
        snp_flag = self._check_ac_repeat_snp(pos, ref, seq)
        if snp_flag:
            flags.append(snp_flag)

        return flags

    def flag_variants(self) -> list[VariantFlag]:
        """Flag each variant individually and return per-variant flag results.

        Unlike analyze() which returns aggregate sample-level reasons, this
        method returns one VariantFlag per variant with the specific flags
        that apply to that variant.

        Returns:
            List of VariantFlag objects, one per input variant.

        """
        if not self.variants:
            return []

        variant_flags: list[VariantFlag] = []

        for v in self.variants:
            pos: str | int = v.get("pos", 0)
            ref: str = v.get("ref", "")
            seq: str = v.get("seq", "")

            flags = self._collect_variant_flags(v)
            if flags is None:
                variant_flags.append(VariantFlag(pos=pos, ref=ref, seq=seq, flags=[]))
                continue

            variant_flags.append(VariantFlag(pos=pos, ref=ref, seq=seq, flags=flags))

        # Sample-level flags that don't map to individual variants
        # (These are handled by analyze() for sample-level aggregation,
        #  but are not per-variant since they span multiple variants or are
        #  pattern-level checks.)

        # No 315.1 variant — this is an absence check, not a per-variant flag
        # 309T, 309Y — per-variant (position-specific)
        # 310, 460 — per-variant (position-specific)
        # 16293A-C — per-variant
        # Consecutive indels — pattern-level, not per-variant
        # Lowercase — per-variant (but detected on seq field)
        # Complex conditions — pattern-level

        # Add per-variant flags for specific patterns
        self._apply_specific_variant_flags(variant_flags)

        return variant_flags

    def _apply_specific_variant_flags(self, variant_flags: list[VariantFlag]) -> None:
        """Apply position-specific and lowercase flags to existing VariantFlags."""
        specific_position_flags = _POSITION_VARIANT_FLAGS
        specific_pos_checks = _POSITION_FLAGS
        transversion_checks = _TRANSVERSION_CHECKS

        for vf in variant_flags:
            pos_str = str(vf.pos)

            # 309T, 309Y
            for pos_key, seq_val, label in specific_position_flags:
                if pos_str == pos_key and vf.seq == seq_val and label not in vf.flags:
                    vf.flags.append(label)

            # 310, 460 — emit "Has {pos} variant"
            for pos_key, label in specific_pos_checks:
                if pos_str == pos_key and label not in vf.flags:
                    vf.flags.append(label)

            # Transversion checks
            vf_base = pos_base(vf.pos)
            for target_pos, target_ref, target_seq, label in transversion_checks:
                if vf_base == target_pos and vf.ref == target_ref and vf.seq == target_seq and label not in vf.flags:
                    vf.flags.append(label)

            # Lowercase
            if vf.seq and any(c.islower() for c in str(vf.seq)) and "Lowercase variant detected" not in vf.flags:
                vf.flags.append("Lowercase variant detected")

    def classify(
        self,
        fis_flag: str,
        fis_flag_info: str = "",
        flag_reasons: list[str] | None = None,
        genotypes_str: str = "",
    ) -> int:
        """Classify sample into flag level (0-4).

        Args:
            fis_flag: 'Y' or 'N' flag status.
            fis_flag_info: QC flag info string.
            flag_reasons: List of flag reasons (uses self.analyze() if None).
            genotypes_str: Space-separated genotype string.

        Returns:
            Highest applicable flag level.

        """
        if fis_flag == "N":
            return self.LEVEL_NONE

        if flag_reasons is None:
            _, flag_reasons, _ = self.analyze()

        flag_reasons_str = ", ".join(flag_reasons) if flag_reasons else ""

        levels = []
        if self._has_quality_issues(fis_flag_info):
            levels.append(self.LEVEL_CRITICAL)
        if self._has_complex_variants(flag_reasons_str, genotypes_str):
            levels.append(self.LEVEL_HIGH)
        if self._has_technical_issues(flag_reasons_str):
            levels.append(self.LEVEL_MEDIUM)
        if self._has_polyc_variants(flag_reasons_str):
            levels.append(self.LEVEL_LOW)

        return max(levels) if levels else self.LEVEL_LOW

    def check_no_315_1(self) -> bool:
        """Return True if the "No 315.1 variant" absence flag applies.

        The expected 315.1C insertion sits at base position 315, inside the
        HV2 region (polyC stretch 303-315). The flag fires when 315.1 is
        absent AND the region containing position 315 (HV2) was part of the
        sample's analysis — i.e. that region has non-empty coverage intervals.

        Region-level (rather than exact-position) gating is deliberate: a
        sample that sequenced HV2 but left a gap over 303-315 (a data problem)
        still gets flagged, surfacing the missing 315.1C, while a sample that
        did not analyze HV2 at all (e.g. HV1-only) is left alone. When no
        intervals are supplied (``None``), full coverage is assumed to
        preserve legacy behavior for callers like the FIS comparison paths.
        """
        if self.analyzer.has_variant_at("315.1"):
            return False
        if self.intervals is None:
            return True
        region = get_hv_region_for_position(315)
        return region is not None and bool(self.intervals.get(region))

    def _check_specific_patterns(self) -> list[str]:
        """Check for specific variant patterns."""
        reasons: list[str] = []

        # "No 315.1 variant" is an absence check for the expected 315.1C
        # insertion (position 315 is in the HV2 region). It fires when 315.1
        # is absent and the HV2 region was analyzed (non-empty intervals); see
        # check_no_315_1() for the region-level rationale.
        if self.check_no_315_1():
            reasons.append("No 315.1 variant")

        for pos_str, seq_val, label in _POSITION_VARIANT_FLAGS:
            if self.analyzer.has_variant_at(pos_str, seq_val):
                reasons.append(label)

        for pos_str, label in _POSITION_FLAGS:
            if self.analyzer.has_variant_at(pos_str):
                reasons.append(label)

        transversion_checks = _TRANSVERSION_CHECKS
        for target_pos, target_ref, target_seq, label in transversion_checks:
            if any(
                pos_base(v.get("pos", 0)) == target_pos
                and v.get("ref") == target_ref
                and v.get("seq", "") == target_seq
                for v in self.variants
            ):
                reasons.append(label)

        return reasons

    def _has_quality_issues(self, fis_flag_info: str) -> bool:
        """Check for Level 4 quality issues."""
        if not fis_flag_info:
            return False
        return any(kw.lower() in fis_flag_info.lower() for kw in self.QUALITY_KEYWORDS)

    def _has_polyc_variants(self, flag_reasons_str: str) -> bool:
        """Check for Level 1 polyC region variants."""
        return any(re.search(p, flag_reasons_str) for p in self.POLYC_PATTERNS)

    def _has_technical_issues(self, flag_reasons_str: str) -> bool:
        """Check for Level 2 technical/nomenclature issues."""
        if re.search(self.AC_REPEAT_SNP_PATTERN, flag_reasons_str):
            return True
        return any(re.search(p, flag_reasons_str) for p in self.NOMENCLATURE_PATTERNS)

    def _has_complex_variants(self, flag_reasons_str: str, genotypes_str: str) -> bool:
        """Check for Level 3 complex/difficult variants."""
        if any(re.search(p, flag_reasons_str) for p in self.COMPLEX_PATTERNS):
            return True

        if re.search(r"Deletion at 459", flag_reasons_str) and (
            re.search(r"Heteroplasmy at 16192", flag_reasons_str) or "16192T" in genotypes_str
        ):
            return True

        if re.search(r"Has 460 variant", flag_reasons_str) or "16293C" in genotypes_str:
            return True

        return self.analyzer.has_consecutive_indels(5) or self.analyzer.has_lowercase_variants()


# ---------------------------------------------------------------------------
def flag_variants(variants: Variant | dict | list) -> dict[str, list[str]]:
    """Return per-variant flag reasons as a dict keyed by "pos|ref|alt".

    This is the canonical way to get variant-level flags. Each entry maps
    a variant key to the list of flag reasons specific to that variant
    (e.g. "Insertion C at 309.1", "Heteroplasmy at 16192").

    Args:
        variants: Variant objects, dicts, or list of either.
            Dict input is deprecated; use Variant objects.

    Returns:
        Dict mapping "pos|ref|alt" keys to lists of flag reason strings.
        Only includes variants that have at least one flag.

    """
    flagger = SampleFlagger(variants)
    variant_flag_objs = flagger.flag_variants()
    result: dict[str, list[str]] = {}
    for vf in variant_flag_objs:
        if vf.flags:
            result[vf.key()] = vf.flags
    return result


# ---------------------------------------------------------------------------
def deduplicate_sample_flags(
    sample_flags: list[str],
    variant_flags: dict[str, list[str]],
) -> list[str]:
    """Remove sample-level flags that already exist as per-variant flags.

    Per-variant flags belong only in ``variant_flags``. Any flag string that
    also appears among the per-variant flag values (e.g. "Deletion at 16193",
    "Insertion C at 309.1", "Heteroplasmy at 16192") is dropped from the
    sample-level list to avoid duplication.

    Aggregated sample-level flags use a different string from their
    per-variant form and are preserved. For example, the per-variant
    "16180-16193 region" stays in ``variant_flags`` while the aggregated
    "16180-16193 region (pos1, pos2)" stays in ``sample_flags``.

    Args:
        sample_flags: Sample-level flag reasons (e.g. from ``analyze()``).
        variant_flags: Per-variant flag reasons keyed by "pos|ref|alt".

    Returns:
        Filtered sample-level flag list, preserving order and dropping
        exact duplicates of any per-variant flag string.

    """
    per_variant: set[str] = set()
    for reasons in variant_flags.values():
        per_variant.update(reasons)
    # The aggregated "16180-16193 region (pos1, pos2, ...)" sample-level flag is
    # a summary of the per-variant "16180-16193 region" flags. The comparison
    # report consolidates the per-variant entries into the same aggregated
    # form for the Variant Flags column, so drop the aggregated form from
    # sample_flags whenever the per-variant region flag is present.
    drop_region_aggregate = "16180-16193 region" in per_variant
    result: list[str] = []
    for flag in sample_flags:
        if flag in per_variant:
            continue
        if drop_region_aggregate and flag.startswith("16180-16193 region"):
            continue
        result.append(flag)
    return result


# ---------------------------------------------------------------------------
def recompute_no_315_1_flag(
    sample_flags: list[str],
    variants: Variant | dict | list,
    intervals: dict[str, list[list[int]]] | None,
) -> list[str]:
    """Recompute the "No 315.1 variant" absence flag from a variant list.

    ``"No 315.1 variant"`` is a special *absence* flag: unlike per-variant
    flags it reports a variant MISSING from the list, so it must be derived
    from the actual variant set rather than trusted from per-file/per-region
    ETL output (which can be stale — e.g. a per-file read that covered
    position 315 but did not call 315.1C, or a per-region JSON that would
    otherwise inherit the parent sample's flag). It is called at JSON-write
    time (``sample_to_dict``) so each output JSON — the combined 3-region
    JSON and each per-region (HV1 / HV2-3) JSON — reports the flag against
    its own variant set + coverage intervals.

    Drops any existing ``"No 315.1 variant"`` and re-derives it via
    :meth:`SampleFlagger.check_no_315_1` (region-level gating: fires when
    315.1C is absent and the region containing position 315 — HV2 — was
    analyzed). ``intervals=None`` assumes full coverage (legacy behavior);
    an empty intervals dict (no region analyzed) does not flag.

    Args:
        sample_flags: Incoming sample-level flag reasons (may already carry
            a stale "No 315.1 variant").
        variants: Variant set to evaluate (Variant objects, dicts, or a list).
        intervals: Coverage intervals keyed by region, or ``None`` to assume
            full coverage.

    Returns:
        Corrected sample-level flag list with the absence flag recomputed.
    """
    result = [f for f in sample_flags if f != "No 315.1 variant"]
    flagger = SampleFlagger(variants, intervals=intervals)
    if flagger.check_no_315_1():
        result.append("No 315.1 variant")
    return result


# ---------------------------------------------------------------------------
def _is_core_sample_flag(flag: str) -> bool:
    """Return True if *flag* is a core flag produced by :meth:`SampleFlagger.analyze`.

    Uses ``SampleFlagger.FLAG_REASON_MAP`` (the authoritative severity-classification
    patterns) to identify core flags dynamically. Quality keywords
    (``SampleFlagger.QUALITY_KEYWORDS`` — e.g. "Low depth", "Failed position")
    are excluded because they originate from external QC, not ``analyze()``.

    This avoids a separate hardcoded pattern list: when a new flag is added to
    ``analyze()`` it must also be added to ``FLAG_REASON_MAP`` for severity
    classification, so ``recompute_sample_flags`` automatically picks it up.
    """
    if any(qk in flag for qk in SampleFlagger.QUALITY_KEYWORDS):
        return False
    return any(pattern in flag for pattern in SampleFlagger.FLAG_REASON_MAP)


def recompute_sample_flags(
    sample_flags: list[str],
    variants: Variant | dict | list,
    intervals: dict[str, list[list[int]]] | None,
    variant_flags: dict[str, list[str]],
) -> list[str]:
    """Recompute core sample-level flags against the current variant set.

    Per-region JSONs inherit the combined sample's ``sample_flags`` via
    :func:`filter_sample_by_regions` (which carries them forward unchanged).
    Cross-region flags like ``"Complex: 459 deletion with 16192 issue"``
    (requires 459DEL in HV2-3 + 16192T in HV1) are invalid for a per-region
    JSON that only has one of the two variants. This function drops stale
    core flags and replaces them with freshly computed ones, so each JSON
    reports flags against its own variant set.

    Tool-specific flags (not produced by :meth:`SampleFlagger.analyze`) are
    preserved unchanged.

    Args:
        sample_flags: Incoming sample-level flag reasons (may carry stale
            core flags from a different variant set).
        variants: Variant set to evaluate (Variant objects, dicts, or a list).
        intervals: Coverage intervals keyed by region, or ``None`` to assume
            full coverage.
        variant_flags: Per-variant flags for the same variant set (used to
            deduplicate sample-level flags that duplicate per-variant flags).

    Returns:
        Corrected sample-level flag list with core flags recomputed and
        tool-specific flags preserved.
    """
    result = [f for f in sample_flags if not _is_core_sample_flag(f)]
    flagger = SampleFlagger(variants, intervals=intervals)
    if flagger.variants:
        _, analyze_reasons, _ = flagger.analyze()
        core_flags = deduplicate_sample_flags(analyze_reasons, variant_flags)
        for flag in core_flags:
            if flag not in result:
                result.append(flag)
    return result
