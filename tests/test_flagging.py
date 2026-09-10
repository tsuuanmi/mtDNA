"""Tests for src.core.flagging — TASKS.md flags, superset verification, and per-variant coverage."""

from src.core.flagging import (
    SampleFlagger,
    VariantAnalyzer,
    VariantFlag,
    deduplicate_sample_flags,
    flag_variants,
    recompute_no_315_1_flag,
)
from src.core.models import Sample, Tool, Variant
from src.core.sample import filter_sample_by_regions, sample_to_dict

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def make_variant(pos, ref, seq):
    """Create a simple variant dict."""
    return {"pos": pos, "ref": ref, "seq": seq}


def analyze_reasons(variants):
    """Run SampleFlagger.analyze() and return flag reasons."""
    flagger = SampleFlagger(variants)
    _, reasons, _ = flagger.analyze()
    return reasons


# ---------------------------------------------------------------------------
# A. Insertion-C exception (TASKS.md #2)
# ---------------------------------------------------------------------------


class TestInsertionCException:
    """Insertion C at 309/573 (special positions) is suppressed; 16193 C-insertion IS flagged (legacy comparator)."""

    def test_insertion_C_at_16193_flagged(self):
        """16193.1C: insertion C at 16193 IS flagged (matches legacy comparator)."""
        v = make_variant("16193.1", "-", "C")
        reasons = analyze_reasons([v])
        assert "Insertion C at 16193.1" in reasons

    def test_insertion_T_at_16193_flagged(self):
        v = make_variant("16193.1", "-", "T")
        reasons = analyze_reasons([v])
        assert "Insertion T at 16193.1" in reasons

    def test_insertion_A_at_16193_flagged(self):
        v = make_variant("16193.1", "-", "A")
        reasons = analyze_reasons([v])
        assert "Insertion A at 16193.1" in reasons

    def test_insertion_C_at_309_special_position(self):
        """309.1C: insertion C at 309 is excluded per TASKS.md #2 exception."""
        v = make_variant("309.1", "-", "C")
        reasons = analyze_reasons([v])
        assert "Insertion C at 309.1" not in reasons

    def test_insertion_T_at_309_flagged(self):
        """309.1T: only the C insertion is suppressed; other bases are flagged."""
        v = make_variant("309.1", "-", "T")
        reasons = analyze_reasons([v])
        assert "Insertion T at 309.1" in reasons

    def test_insertion_C_at_573_special_position(self):
        """573.1C: insertion C at 573 is excluded per TASKS.md #2 exception."""
        v = make_variant("573.1", "-", "C")
        reasons = analyze_reasons([v])
        assert "Insertion C at 573.1" not in reasons

    def test_insertion_T_at_573_flagged(self):
        """573.1T: only the C insertion is suppressed; other bases are flagged."""
        v = make_variant("573.1", "-", "T")
        reasons = analyze_reasons([v])
        assert "Insertion T at 573.1" in reasons

    def test_insertion_C_at_other_position_flagged(self):
        """Insertion C at positions NOT in the exception set is still flagged."""
        v = make_variant("290.1", "-", "C")
        reasons = analyze_reasons([v])
        assert "Insertion C at 290.1" in reasons


# ---------------------------------------------------------------------------
# B. 16258A-C transversion (TASKS.md #3)
# ---------------------------------------------------------------------------


class Test16258ACTransversion:
    def test_16258AC_flagged(self):
        v = make_variant(16258, "A", "C")
        reasons = analyze_reasons([v])
        assert "Has 16258A-C variant" in reasons

    def test_16258AC_not_flagged_other_seq(self):
        v = make_variant(16258, "A", "G")
        reasons = analyze_reasons([v])
        assert "Has 16258A-C variant" not in reasons

    def test_16258AC_not_flagged_other_ref(self):
        v = make_variant(16258, "C", "T")
        reasons = analyze_reasons([v])
        assert "Has 16258A-C variant" not in reasons

    def test_16258AC_in_per_variant_flags(self):
        v = make_variant(16258, "A", "C")
        result = flag_variants([v])
        key = "16258|A|C"
        assert key in result
        assert "Has 16258A-C variant" in result[key]


# ---------------------------------------------------------------------------
# C. 16263T-C transition (TASKS.md #4)
# ---------------------------------------------------------------------------


class Test16263TCTransition:
    def test_16263TC_flagged(self):
        v = make_variant(16263, "T", "C")
        reasons = analyze_reasons([v])
        assert "Has 16263T-C variant" in reasons

    def test_16263TC_not_flagged_other_seq(self):
        v = make_variant(16263, "T", "A")
        reasons = analyze_reasons([v])
        assert "Has 16263T-C variant" not in reasons

    def test_16263TC_not_flagged_other_ref(self):
        v = make_variant(16263, "C", "T")
        reasons = analyze_reasons([v])
        assert "Has 16263T-C variant" not in reasons

    def test_16263TC_in_per_variant_flags(self):
        v = make_variant(16263, "T", "C")
        result = flag_variants([v])
        key = "16263|T|C"
        assert key in result
        assert "Has 16263T-C variant" in result[key]


# ---------------------------------------------------------------------------
# D. Superset verification: every expected flag reason should be
#    produced by SampleFlagger.analyze()
# ---------------------------------------------------------------------------


class TestSuperset:
    """Verify that SampleFlagger.analyze() produces every expected flag reason."""

    def test_variant_at_310(self):
        """Only 'Has 310 variant' is emitted (duplicate 'Variant at position 310' removed)."""
        v = make_variant(310, "T", "C")
        reasons = analyze_reasons([v])
        assert "Has 310 variant" in reasons
        assert "Variant at position 310" not in reasons

    def test_variant_at_460(self):
        """Only 'Has 460 variant' is emitted (duplicate 'Variant at position 460' removed)."""
        v = make_variant(460, "T", "C")
        reasons = analyze_reasons([v])
        assert "Has 460 variant" in reasons
        assert "Variant at position 460" not in reasons

    def test_insertion_at_non_special(self):
        """Insertion at a non-special position is flagged."""
        v = make_variant("290.1", "-", "T")
        reasons = analyze_reasons([v])
        assert "Insertion T at 290.1" in reasons

    def test_deletion(self):
        v = make_variant(200, "A", "-")
        reasons = analyze_reasons([v])
        assert "Deletion at 200" in reasons

    def test_deletion_249_excluded(self):
        v = make_variant(249, "A", "-")
        reasons = analyze_reasons([v])
        assert "Deletion at 249" not in reasons

    def test_heteroplasmy(self):
        v = make_variant(100, "A", "R")
        reasons = analyze_reasons([v])
        assert "Heteroplasmy at 100" in reasons

    def test_16180_16193_region(self):
        """SNP in the 16180-16193 region is flagged."""
        v = make_variant(16185, "A", "G")
        reasons = analyze_reasons([v])
        assert any(r.startswith("16180-16193 region") for r in reasons)

    def test_no_315_1_variant(self):
        v = make_variant(100, "A", "G")
        reasons = analyze_reasons([v])
        assert "No 315.1 variant" in reasons

    def test_has_309T(self):
        v = make_variant(309, "A", "T")
        reasons = analyze_reasons([v])
        assert "Has 309T variant" in reasons

    def test_has_309Y(self):
        v = make_variant(309, "A", "Y")
        reasons = analyze_reasons([v])
        assert "Has 309Y variant" in reasons

    def test_has_310_variant(self):
        v = make_variant(310, "A", "G")
        reasons = analyze_reasons([v])
        assert "Has 310 variant" in reasons

    def test_has_460_variant(self):
        v = make_variant(460, "A", "G")
        reasons = analyze_reasons([v])
        assert "Has 460 variant" in reasons

    def test_has_16293AC_variant(self):
        v = make_variant(16293, "A", "C")
        reasons = analyze_reasons([v])
        assert "Has 16293A-C variant" in reasons


# ---------------------------------------------------------------------------
# E. flag_variants() per-variant coverage
# ---------------------------------------------------------------------------


class TestFlagVariantsPerVariant:
    def test_flag_variants_16258AC(self):
        v = make_variant(16258, "A", "C")
        result = flag_variants([v])
        key = "16258|A|C"
        assert key in result
        assert "Has 16258A-C variant" in result[key]

    def test_flag_variants_16263TC(self):
        v = make_variant(16263, "T", "C")
        result = flag_variants([v])
        key = "16263|T|C"
        assert key in result
        assert "Has 16263T-C variant" in result[key]

    def test_flag_variants_insertion_C_16193(self):
        """Insertion C at 16193.1 IS flagged (matches legacy comparator)."""
        v = make_variant("16193.1", "-", "C")
        result = flag_variants([v])
        key = "16193.1|-|C"
        assert key in result
        assert "Insertion C at 16193.1" in result[key]

    def test_flag_variants_insertion_T_16193(self):
        """Insertion T at 16193.1 should have Insertion flag."""
        v = make_variant("16193.1", "-", "T")
        result = flag_variants([v])
        key = "16193.1|-|T"
        assert key in result
        assert "Insertion T at 16193.1" in result[key]


# ---------------------------------------------------------------------------
# F. Regression — existing flags still work
# ---------------------------------------------------------------------------


class TestRegression:
    def test_existing_flags_unchanged(self):
        """Comprehensive test that existing flag reasons still work."""
        variants = [
            make_variant(310, "T", "C"),  # Has 310 variant
            make_variant("290.1", "-", "T"),  # Insertion at 290.1
            make_variant(200, "A", "-"),  # Deletion at 200
            make_variant(100, "A", "R"),  # Heteroplasmy at 100
            make_variant(16185, "A", "G"),  # 16180-16193 region
            make_variant(309, "A", "T"),  # Has 309T variant
            make_variant(16293, "A", "C"),  # Has 16293A-C variant
            make_variant(460, "A", "G"),  # Has 460 variant
        ]
        reasons = analyze_reasons(variants)

        assert "Has 310 variant" in reasons
        assert "Insertion T at 290.1" in reasons
        assert "Deletion at 200" in reasons
        assert "Heteroplasmy at 100" in reasons
        assert any(r.startswith("16180-16193 region") for r in reasons)
        assert "Has 309T variant" in reasons
        assert "Has 16293A-C variant" in reasons
        assert "Has 460 variant" in reasons

    def test_ac_repeat_snp_flagging(self):
        """SNP in AC repeat region (513-525) is flagged."""
        v = make_variant(519, "A", "G")
        reasons = analyze_reasons([v])
        assert "Has 519A-G variant" in reasons

    def test_523_524_ac_deletion_exception(self):
        """Deletions at 523+524 are excluded (AC deletion exception)."""
        variants = [
            make_variant(523, "A", "-"),
            make_variant(524, "C", "-"),
        ]
        reasons = analyze_reasons(variants)
        assert "Deletion at 523" not in reasons
        assert "Deletion at 524" not in reasons

    def test_513a_with_523_524_exception(self):
        """SNP at 513 with ref=A excluded when 523-524 AC deletion is present."""
        variants = [
            make_variant(513, "A", "G"),
            make_variant(523, "A", "-"),
            make_variant(524, "C", "-"),
        ]
        reasons = analyze_reasons(variants)
        # 513A SNP should be excluded when paired with 523-524 deletion
        assert not any(r.startswith("Has 513") and "variant" in r for r in reasons)

    def test_region_flag_works(self):
        """16180-16193 region flag is produced for variants in that range."""
        v = make_variant(16185, "A", "G")
        reasons = analyze_reasons([v])
        assert any(r.startswith("16180-16193 region") for r in reasons)


# ---------------------------------------------------------------------------
# G. Comprehensive superset test (critical for AC#5)
# ---------------------------------------------------------------------------


class TestComprehensiveSuperset:
    """Every expected flag reason must also be produced by SampleFlagger.analyze()."""

    def test_analyze_produces_all_expected_flag_strings(self):
        """Test that analyze() produces all expected flag strings."""
        variants = [
            make_variant(310, "T", "C"),  # Has 310 variant
            make_variant(460, "A", "G"),  # Has 460 variant
            make_variant("290.1", "-", "T"),  # Insertion at 290.1
            make_variant(200, "A", "-"),  # Deletion at 200
            make_variant(100, "A", "R"),  # Heteroplasmy at 100
            make_variant(16185, "A", "G"),  # 16180-16193 region
            make_variant(309, "A", "T"),  # Has 309T variant
            make_variant(16293, "A", "C"),  # Has 16293A-C variant
        ]
        reasons = analyze_reasons(variants)

        # Check all expected flag reason patterns
        assert "Has 310 variant" in reasons
        assert "Has 460 variant" in reasons
        assert "Insertion T at 290.1" in reasons
        assert "Deletion at 200" in reasons
        assert "Heteroplasmy at 100" in reasons
        assert any(r.startswith("16180-16193 region") for r in reasons)
        assert "Has 309T variant" in reasons
        assert "Has 16293A-C variant" in reasons
        assert "No 315.1 variant" in reasons

    def test_16258AC_and_16263TC_in_comprehensive(self):
        """New TASKS.md flags also present."""
        variants = [
            make_variant(16258, "A", "C"),
            make_variant(16263, "T", "C"),
        ]
        reasons = analyze_reasons(variants)
        assert "Has 16258A-C variant" in reasons
        assert "Has 16263T-C variant" in reasons


# ---------------------------------------------------------------------------
# H. Recompute "No 315.1 variant" on merged samples (per-file union fix)
# ---------------------------------------------------------------------------


class TestNo3151RegionGating:
    """The "No 315.1 variant" absence flag is gated by region, not by exact
    position coverage: it fires when 315.1 is absent AND the HV2 region (which
    contains position 315) was analyzed."""

    def test_flagged_when_hv2_analyzed_no_315_1(self):
        flagger = SampleFlagger([make_variant(73, "A", "G")], intervals={"HV2": [[73, 340]]})
        assert flagger.check_no_315_1() is True

    def test_flagged_when_hv2_analyzed_with_gap_at_315(self):
        # Full-region sample that left a gap over 303-315 (a data problem):
        # the reviewer still wants "No 315.1 variant" flagged.
        flagger = SampleFlagger([make_variant(73, "A", "G")], intervals={"HV2": [[73, 302], [316, 340]]})
        assert flagger.check_no_315_1() is True

    def test_not_flagged_when_hv2_not_analyzed(self):
        # HV1-only sample: HV2 was never sequenced, so 315.1 cannot be expected.
        flagger = SampleFlagger([make_variant(16223, "A", "T")], intervals={"HV1": [[16024, 16365]]})
        assert flagger.check_no_315_1() is False

    def test_not_flagged_when_hv2_empty(self):
        flagger = SampleFlagger(
            [make_variant(16223, "A", "T")],
            intervals={"HV1": [[16024, 16365]], "HV2": [], "HV3": []},
        )
        assert flagger.check_no_315_1() is False

    def test_not_flagged_when_315_1_present(self):
        flagger = SampleFlagger(
            [make_variant("315.1", "-", "C")],
            intervals={"HV2": [[73, 340]]},
        )
        assert flagger.check_no_315_1() is False

    def test_flagged_when_intervals_none(self):
        # Legacy/FIS path: no intervals => full coverage assumed.
        flagger = SampleFlagger([make_variant(100, "A", "G")])
        assert flagger.check_no_315_1() is True

    def test_not_flagged_when_no_region_analyzed(self):
        # Empty intervals (no region analyzed) must not flag — it is not the
        # same as the legacy intervals=None "full coverage" path.
        flagger = SampleFlagger([make_variant(100, "A", "G")], intervals={})
        assert flagger.check_no_315_1() is False


class TestRecomputeNo3151Flag:
    """The absence flag must be derived from the merged variant set, not the
    per-file flag union."""

    def test_false_positive_removed_when_315_1_present(self):
        """Per-file union leaves a stale 'No 315.1 variant' from a read that
        covered 315 but did not call 315.1C; the merged set has 315.1C."""
        variants = [make_variant("315.1", "-", "C"), make_variant(73, "A", "G")]
        intervals = {"HV2": [[73, 340]]}
        result = recompute_no_315_1_flag(["No 315.1 variant"], variants, intervals)
        assert result == []

    def test_flag_kept_when_hv2_analyzed_and_315_1_absent(self):
        variants = [make_variant(73, "A", "G")]
        intervals = {"HV2": [[73, 340]]}
        result = recompute_no_315_1_flag([], variants, intervals)
        assert result == ["No 315.1 variant"]

    def test_flag_kept_when_hv2_analyzed_with_gap_at_315(self):
        # AA3838-style manual range: HV2 analyzed but with a 303-315 gap.
        variants = [make_variant(73, "A", "G")]
        intervals = {"HV2": [[73, 302], [316, 340]]}
        result = recompute_no_315_1_flag([], variants, intervals)
        assert result == ["No 315.1 variant"]

    def test_flag_not_added_when_hv2_not_analyzed(self):
        variants = [make_variant(16223, "A", "T")]
        intervals = {"HV1": [[16024, 16365]]}
        result = recompute_no_315_1_flag(["No 315.1 variant"], variants, intervals)
        assert result == []

    def test_other_flags_preserved(self):
        variants = [make_variant("315.1", "-", "C"), make_variant(309, "A", "T")]
        intervals = {"HV2": [[73, 340]]}
        result = recompute_no_315_1_flag(["No 315.1 variant", "Has 309T variant"], variants, intervals)
        assert "No 315.1 variant" not in result
        assert "Has 309T variant" in result

    def test_not_added_when_no_region_analyzed(self):
        # Empty intervals dict (no region analyzed) != intervals=None (full
        # coverage): must not add the absence flag.
        variants = [make_variant(100, "A", "G")]
        result = recompute_no_315_1_flag(["No 315.1 variant"], variants, {})
        assert result == []


# ---------------------------------------------------------------------------
# I. sample_to_dict computes "No 315.1 variant" per JSON (combined vs region)
# ---------------------------------------------------------------------------


class TestSampleToDictNo3151:
    """The absence flag must be derived from each output JSON's own variant
    list, not inherited from the parent sample or the per-file flag union."""

    @staticmethod
    def _sample(sample_id, variant_specs, intervals, sample_flags=None) -> Sample:
        variants = [Variant(pos=p, ref=r, seq=s) for (p, r, s) in variant_specs]
        return Sample(
            sample_id=sample_id,
            source_tool=Tool.TRACY,
            variants=variants,
            intervals=intervals,
            sample_flags=list(sample_flags or []),
            variant_flags={},
            information=None,
            batch_id=None,
            hv1=None,
            hv2=None,
            hv3=None,
            no_snps=None,
            no_ins=None,
            no_dels=None,
            no_identicals=None,
            no_unread=None,
        )

    def test_combined_json_flags_when_hv2_analyzed_no_315_1(self, ref_sequence):
        s = self._sample(
            "S1",
            [(73, "A", "G"), (16223, "C", "T")],
            {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
        )
        d = sample_to_dict(s, ref_sequence)
        assert "No 315.1 variant" in d["sample_flags"]

    def test_combined_json_flags_with_gap_at_315(self, ref_sequence):
        """AA3838-style: HV2 analyzed but with a 303-315 gap -> still flagged."""

        s = self._sample(
            "S1",
            [(73, "A", "G")],
            {"HV1": [[16024, 16365]], "HV2": [[73, 302], [316, 340]], "HV3": [[438, 576]]},
            sample_flags=["No 315.1 variant"],
        )
        d = sample_to_dict(s, ref_sequence)
        assert "No 315.1 variant" in d["sample_flags"]

    def test_combined_json_not_flagged_when_315_1_present(self, ref_sequence):
        s = self._sample(
            "S1",
            [("315.1", "-", "C"), (73, "A", "G")],
            {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
            sample_flags=["No 315.1 variant"],
        )
        d = sample_to_dict(s, ref_sequence)
        assert "No 315.1 variant" not in d["sample_flags"]

    def test_hv1_region_json_not_flagged(self, ref_sequence):
        s = self._sample(
            "S1",
            [(73, "A", "G"), (16223, "C", "T")],
            {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
            sample_flags=["No 315.1 variant"],
        )
        hv1 = filter_sample_by_regions(s, ["HV1"])
        assert hv1 is not None
        d = sample_to_dict(hv1, ref_sequence)
        assert "No 315.1 variant" not in d["sample_flags"]

    def test_hv2_3_region_json_flags_when_no_315_1(self, ref_sequence):
        s = self._sample(
            "S1",
            [(73, "A", "G"), (490, "T", "C")],
            {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
            sample_flags=["No 315.1 variant"],
        )
        hv23 = filter_sample_by_regions(s, ["HV2", "HV3"])
        assert hv23 is not None
        d = sample_to_dict(hv23, ref_sequence)
        assert "No 315.1 variant" in d["sample_flags"]

    def test_hv1_only_sample_not_flagged(self, ref_sequence):
        s = self._sample("S1", [(16223, "C", "T")], {"HV1": [[16024, 16365]]})
        d = sample_to_dict(s, ref_sequence)
        assert "No 315.1 variant" not in d["sample_flags"]


# ---------------------------------------------------------------------------
# J. classify() severity-level routing
# ---------------------------------------------------------------------------


class TestClassify:
    """SampleFlagger.classify() maps flag reasons + FIS state to severity 0-4."""

    def test_fis_n_returns_none(self):
        # "N" sample never flagged regardless of reasons.
        assert SampleFlagger([]).classify("N", flag_reasons=["Has 460 variant"]) == SampleFlagger.LEVEL_NONE

    def test_quality_issue_is_critical(self):
        # QC keywords in fis_flag_info -> Level 4.
        for kw in ("Low depth", "Failed position", "Possible wrong consensus"):
            assert SampleFlagger([]).classify("Y", fis_flag_info=kw) == SampleFlagger.LEVEL_CRITICAL

    def test_quality_issue_case_insensitive(self):
        assert SampleFlagger([]).classify("Y", fis_flag_info="low DEPTH here") == SampleFlagger.LEVEL_CRITICAL

    def test_quality_info_empty_not_critical(self):
        assert SampleFlagger([]).classify("Y", fis_flag_info="") == SampleFlagger.LEVEL_LOW

    def test_complex_460_is_high(self):
        # "Has 460 variant" matches COMPLEX_PATTERNS -> Level 3.
        assert SampleFlagger([]).classify("Y", flag_reasons=["Has 460 variant"]) == SampleFlagger.LEVEL_HIGH

    def test_complex_310_is_high(self):
        assert SampleFlagger([]).classify("Y", flag_reasons=["Has 310 variant"]) == SampleFlagger.LEVEL_HIGH

    def test_complex_no_315_1_is_high(self):
        assert SampleFlagger([]).classify("Y", flag_reasons=["No 315.1 variant"]) == SampleFlagger.LEVEL_HIGH

    def test_complex_459_16192_genotype_is_high(self):
        # 459 deletion in reasons + 16192T genotype string -> complex.
        f = SampleFlagger([make_variant(459, "A", "-")])
        assert f.classify("Y", genotypes_str="16192T") == SampleFlagger.LEVEL_HIGH

    def test_complex_459_16192_heteroplasmy_is_high(self):
        f = SampleFlagger([make_variant(459, "A", "-"), make_variant(16192, "T", "Y")])
        assert f.classify("Y") == SampleFlagger.LEVEL_HIGH

    def test_complex_16293c_genotype_is_high(self):
        assert SampleFlagger([]).classify("Y", flag_reasons=[], genotypes_str="16293C") == SampleFlagger.LEVEL_HIGH

    def test_technical_ac_repeat_snp_is_medium(self):
        # AC repeat SNP (519) -> technical Level 2 (isolated from complex).
        f = SampleFlagger([make_variant(519, "A", "G"), make_variant("315.1", "-", "C")])
        assert f.classify("Y") == SampleFlagger.LEVEL_MEDIUM

    def test_technical_nomenclature_290_insertions_is_medium(self):
        f = SampleFlagger(
            [make_variant("290.1", "-", "C"), make_variant("290.2", "-", "C"), make_variant("315.1", "-", "C")]
        )
        assert f.classify("Y") == SampleFlagger.LEVEL_MEDIUM

    def test_polyc_region_is_low(self):
        # 16180-16193 region only (with 315.1 present to avoid complex) -> Level 1.
        f = SampleFlagger([make_variant(16185, "A", "G"), make_variant("315.1", "-", "C")])
        assert f.classify("Y") == SampleFlagger.LEVEL_LOW

    def test_default_when_no_reasons_is_low(self):
        # Flagged sample with no recognizable patterns -> baseline Level 1.
        assert SampleFlagger([make_variant("315.1", "-", "C")]).classify("Y") == SampleFlagger.LEVEL_LOW

    def test_max_level_wins_critical_over_high(self):
        assert (
            SampleFlagger([]).classify("Y", fis_flag_info="Low depth", flag_reasons=["Has 460 variant"])
            == SampleFlagger.LEVEL_CRITICAL
        )


# ---------------------------------------------------------------------------
# K. level_from_reasons() direct mapping
# ---------------------------------------------------------------------------


class TestLevelFromReasons:
    def test_empty_is_none(self):
        assert SampleFlagger.level_from_reasons([]) == SampleFlagger.LEVEL_NONE

    def test_critical_keywords(self):
        for kw in ("Low depth", "Failed position", "Possible wrong consensus"):
            assert SampleFlagger.level_from_reasons([kw]) == SampleFlagger.LEVEL_CRITICAL

    def test_low_polyc(self):
        assert SampleFlagger.level_from_reasons(["16180-16193 region (16185)"]) == SampleFlagger.LEVEL_LOW
        assert SampleFlagger.level_from_reasons(["Deletion at 523"]) == SampleFlagger.LEVEL_LOW
        assert SampleFlagger.level_from_reasons(["309T"]) == SampleFlagger.LEVEL_LOW

    def test_medium_ac_repeat(self):
        assert SampleFlagger.level_from_reasons(["Has 519A-G variant"]) == SampleFlagger.LEVEL_MEDIUM
        assert SampleFlagger.level_from_reasons(["Insertion C at 290.1"]) == SampleFlagger.LEVEL_MEDIUM

    def test_high_specific(self):
        assert SampleFlagger.level_from_reasons(["Deletion at 315"]) == SampleFlagger.LEVEL_HIGH
        assert SampleFlagger.level_from_reasons(["Deletion at 459"]) == SampleFlagger.LEVEL_HIGH
        assert SampleFlagger.level_from_reasons(["No 315.1 variant"]) == SampleFlagger.LEVEL_HIGH
        assert SampleFlagger.level_from_reasons(["Consecutive indels (5+)"]) == SampleFlagger.LEVEL_HIGH
        assert SampleFlagger.level_from_reasons(["Lowercase variant detected"]) == SampleFlagger.LEVEL_HIGH

    def test_general_deletion_insertion_is_high(self):
        # Generic insertion/deletion at a position NOT in the Level-2 list
        # matches the general "Insertion"/"Deletion" patterns -> Level 3.
        assert SampleFlagger.level_from_reasons(["Insertion C at 500.1"]) == SampleFlagger.LEVEL_HIGH
        assert SampleFlagger.level_from_reasons(["Deletion at 200"]) == SampleFlagger.LEVEL_HIGH

    def test_max_across_multiple(self):
        # Critical keyword dominates a high-level deletion.
        assert SampleFlagger.level_from_reasons(["Low depth", "Deletion at 200"]) == SampleFlagger.LEVEL_CRITICAL


# ---------------------------------------------------------------------------
# L. deduplicate_sample_flags helper
# ---------------------------------------------------------------------------


class TestDeduplicateSampleFlags:
    def test_drops_flags_present_in_variant_flags(self):
        sample_flags = ["Deletion at 16193", "Insertion C at 309.1", "Has 310 variant"]
        variant_flags = {
            "16193|A|-": ["Deletion at 16193"],
            "309.1|-|C": ["Insertion C at 309.1"],
        }
        result = deduplicate_sample_flags(sample_flags, variant_flags)
        assert "Deletion at 16193" not in result
        assert "Insertion C at 309.1" not in result
        assert "Has 310 variant" in result

    def test_drops_aggregated_region_when_per_variant_region_present(self):
        sample_flags = ["16180-16193 region (16185, 16190)", "Has 310 variant"]
        variant_flags = {"16185|A|G": ["16180-16193 region"]}
        result = deduplicate_sample_flags(sample_flags, variant_flags)
        assert not any(r.startswith("16180-16193 region") for r in result)
        assert "Has 310 variant" in result

    def test_keeps_aggregated_region_when_no_per_variant_region(self):
        sample_flags = ["16180-16193 region (16185)", "Has 310 variant"]
        variant_flags = {"16185|A|G": ["Heteroplasmy at 16185"]}
        result = deduplicate_sample_flags(sample_flags, variant_flags)
        assert any(r.startswith("16180-16193 region") for r in result)

    def test_preserves_order_of_non_duplicate_flags(self):
        sample_flags = ["Has 310 variant", "Has 460 variant", "Has 309T variant"]
        variant_flags = {"16193|A|-": ["Deletion at 16193"]}
        result = deduplicate_sample_flags(sample_flags, variant_flags)
        assert result == ["Has 310 variant", "Has 460 variant", "Has 309T variant"]

    def test_empty_inputs(self):
        assert deduplicate_sample_flags([], {}) == []


# ---------------------------------------------------------------------------
# M. VariantAnalyzer direct unit tests
# ---------------------------------------------------------------------------


class TestVariantAnalyzerConsecutiveIndels:
    def test_five_consecutive_deletions(self):
        a = VariantAnalyzer([make_variant(p, "A", "-") for p in range(100, 105)])
        assert a.has_consecutive_indels(5) is True

    def test_four_consecutive_deletions_false(self):
        a = VariantAnalyzer([make_variant(p, "A", "-") for p in range(100, 104)])
        assert a.has_consecutive_indels(5) is False

    def test_non_consecutive_deletions_false(self):
        a = VariantAnalyzer([make_variant(100, "A", "-"), make_variant(102, "A", "-"), make_variant(104, "A", "-")])
        assert a.has_consecutive_indels(5) is False

    def test_five_consecutive_insertions_same_base(self):
        a = VariantAnalyzer([make_variant(f"100.{i}", "-", "C") for i in range(1, 6)])
        assert a.has_consecutive_indels(5) is True

    def test_insertions_different_bases_false(self):
        a = VariantAnalyzer([make_variant(f"{b}.1", "-", "C") for b in (100, 101, 102, 103, 104)])
        assert a.has_consecutive_indels(5) is False

    def test_count_le_zero_false(self):
        assert VariantAnalyzer([make_variant(100, "A", "-")]).has_consecutive_indels(0) is False

    def test_empty_false(self):
        assert VariantAnalyzer([]).has_consecutive_indels(5) is False

    def test_duplicates_collapsed_in_run(self):
        # Duplicates should not break a run (collapsed to unique set).
        a = VariantAnalyzer(
            [
                make_variant(100, "A", "-"),
                make_variant(100, "A", "-"),
                make_variant(101, "A", "-"),
                make_variant(102, "A", "-"),
                make_variant(103, "A", "-"),
                make_variant(104, "A", "-"),
            ]
        )
        assert a.has_consecutive_indels(5) is True


class TestVariantAnalyzerLowercase:
    def test_lowercase_seq_detected(self):
        assert VariantAnalyzer([make_variant(100, "A", "g")]).has_lowercase_variants() is True

    def test_uppercase_not_detected(self):
        assert VariantAnalyzer([make_variant(100, "A", "G")]).has_lowercase_variants() is False

    def test_empty_not_detected(self):
        assert VariantAnalyzer([]).has_lowercase_variants() is False

    def test_mixed_case_detected(self):
        a = VariantAnalyzer([make_variant(100, "A", "G"), make_variant(200, "T", "c")])
        assert a.has_lowercase_variants() is True


class TestVariantAnalyzerHasVariantAt:
    def test_match_without_seq(self):
        assert VariantAnalyzer([make_variant(309, "A", "T")]).has_variant_at(309) is True

    def test_match_with_seq(self):
        assert VariantAnalyzer([make_variant(309, "A", "T")]).has_variant_at(309, "T") is True

    def test_no_match_with_seq(self):
        assert VariantAnalyzer([make_variant(309, "A", "T")]).has_variant_at(309, "G") is False

    def test_no_match_position(self):
        assert VariantAnalyzer([make_variant(309, "A", "T")]).has_variant_at(310) is False

    def test_insertion_subposition_match(self):
        assert VariantAnalyzer([make_variant("315.1", "-", "C")]).has_variant_at("315.1") is True


class TestVariantAnalyzerInitForms:
    def test_single_variant_object(self):
        v = Variant(pos=309, ref="A", seq="T")
        a = VariantAnalyzer(v)
        assert a.has_variant_at(309, "T") is True

    def test_single_dict(self):
        a = VariantAnalyzer(make_variant(309, "A", "T"))
        assert a.has_variant_at(309, "T") is True

    def test_list_of_dicts(self):
        a = VariantAnalyzer([make_variant(310, "T", "C"), make_variant(460, "A", "G")])
        assert a.has_variant_at(310) is True
        assert a.has_variant_at(460) is True


# ---------------------------------------------------------------------------
# N. Complex combined condition: 459 deletion + 16192 issue
# ---------------------------------------------------------------------------


class TestComplexCombinedCondition:
    def test_459_del_with_16192_heteroplasmy(self):
        reasons = analyze_reasons([make_variant(459, "A", "-"), make_variant(16192, "T", "Y")])
        assert "Complex: 459 deletion with 16192 issue" in reasons

    def test_459_del_with_16192T_snp(self):
        reasons = analyze_reasons([make_variant(459, "A", "-"), make_variant(16192, "T", "T")])
        assert "Complex: 459 deletion with 16192 issue" in reasons

    def test_459_del_without_16192_issue_no_complex(self):
        reasons = analyze_reasons([make_variant(459, "A", "-"), make_variant(16223, "C", "T")])
        assert "Complex: 459 deletion with 16192 issue" not in reasons

    def test_16192_issue_without_459_del_no_complex(self):
        reasons = analyze_reasons([make_variant(16192, "T", "Y")])
        assert "Complex: 459 deletion with 16192 issue" not in reasons


# ---------------------------------------------------------------------------
# O. 309DEL (base-position 309 deletion) edge case
# ---------------------------------------------------------------------------


class Test309Deletion:
    """The 309 base-position deletion (Tracy 309C>T + 310T>C) is a real variant
    event and must be flagged, unlike 309.1 insertions which are suppressed."""

    def test_309_deletion_flagged(self):
        reasons = analyze_reasons([make_variant(309, "A", "-")])
        assert "Deletion at 309" in reasons

    def test_309_1_insertion_still_suppressed(self):
        reasons = analyze_reasons([make_variant("309.1", "-", "C")])
        assert "Insertion C at 309.1" not in reasons

    def test_309_deletion_in_per_variant_flags(self):
        result = flag_variants([make_variant(309, "A", "-")])
        key = "309|A|-"
        assert key in result
        assert "Deletion at 309" in result[key]


# ---------------------------------------------------------------------------
# P. Insertion at base position 315 (315.x) suppression
# ---------------------------------------------------------------------------


class Test315InsertionSuppression:
    """Only the expected 315.1C insertion (normal rCRS reference insertion) is
    suppressed. Extra insertions at 315.2/315.3 and other 315.x are flagged."""

    def test_315_1C_suppressed(self):
        assert "Insertion C at 315.1" not in analyze_reasons([make_variant("315.1", "-", "C")])

    def test_315_1T_flagged(self):
        # Only 315.1C is suppressed; a different base at 315.1 is flagged.
        assert "Insertion T at 315.1" in analyze_reasons([make_variant("315.1", "-", "T")])

    def test_315_2_insertion_flagged(self):
        assert "Insertion C at 315.2" in analyze_reasons([make_variant("315.2", "-", "C")])

    def test_315_2T_insertion_flagged(self):
        assert "Insertion T at 315.2" in analyze_reasons([make_variant("315.2", "-", "T")])

    def test_315_3_insertion_flagged(self):
        assert "Insertion C at 315.3" in analyze_reasons([make_variant("315.3", "-", "C")])

    def test_315_base_insertion_flagged(self):
        assert "Insertion C at 315" in analyze_reasons([make_variant(315, "-", "C")])

    def test_316_insertion_still_flagged(self):
        # Neighbouring position 316 is not suppressed.
        assert "Insertion C at 316.1" in analyze_reasons([make_variant("316.1", "-", "C")])

    def test_315_2C_in_per_variant_flags(self):
        result = flag_variants([make_variant("315.2", "-", "C")])
        assert "315.2|-|C" in result
        assert "Insertion C at 315.2" in result["315.2|-|C"]

    def test_315_1C_not_in_per_variant_flags(self):
        result = flag_variants([make_variant("315.1", "-", "C")])
        assert "315.1|-|C" not in result


# ---------------------------------------------------------------------------
# Q. AC repeat region SNP edge cases (513-525)
# ---------------------------------------------------------------------------


class TestACRepeatSNPEdges:
    def test_snp_at_513_boundary(self):
        assert "Has 513A-G variant" in analyze_reasons([make_variant(513, "A", "G")])

    def test_snp_at_525_boundary(self):
        assert "Has 525A-G variant" in analyze_reasons([make_variant(525, "A", "G")])

    def test_snp_outside_region_not_flagged(self):
        # 512 and 526 are outside 513-525.
        assert "Has 512" not in " ".join(analyze_reasons([make_variant(512, "A", "G")]))
        assert "Has 526" not in " ".join(analyze_reasons([make_variant(526, "A", "G")]))

    def test_insertion_in_region_no_snp_flag(self):
        # Insertion in 513-525 produces an Insertion flag, not a SNP flag.
        reasons = analyze_reasons([make_variant("513.1", "-", "C")])
        assert "Insertion C at 513.1" in reasons
        assert not any(r.startswith("Has 513") and "variant" in r for r in reasons)

    def test_deletion_in_region_no_snp_flag(self):
        # Deletion at 513 (without 523/524 AC pair) produces Deletion flag, not SNP.
        reasons = analyze_reasons([make_variant(513, "A", "-")])
        assert "Deletion at 513" in reasons
        assert not any(r.startswith("Has 513") and "variant" in r for r in reasons)

    def test_513a_suppressed_with_523_524_ac_deletion(self):
        variants = [make_variant(513, "A", "G"), make_variant(523, "A", "-"), make_variant(524, "C", "-")]
        reasons = analyze_reasons(variants)
        assert not any(r.startswith("Has 513") and "variant" in r for r in reasons)


# ---------------------------------------------------------------------------
# R. Empty variant set
# ---------------------------------------------------------------------------


class TestEmptyVariantSet:
    def test_analyze_empty(self):
        should_flag, reasons, level = SampleFlagger([]).analyze()
        assert should_flag is False
        assert reasons == []
        assert level == SampleFlagger.LEVEL_NONE

    def test_flag_variants_empty(self):
        assert flag_variants([]) == {}

    def test_analyze_no_variants_no_no_315_1(self):
        # With intervals=None, empty variants still flags "No 315.1 variant"
        # because 315.1 is absent and full coverage assumed. This confirms
        # the empty-list short-circuit in analyze() runs before the absence
        # check (analyze returns early for no variants).
        _should_flag, reasons, _ = SampleFlagger([]).analyze()
        assert "No 315.1 variant" not in reasons


# ---------------------------------------------------------------------------
# S. 16180-16193 region aggregation
# ---------------------------------------------------------------------------


class TestRegionAggregation:
    def test_multiple_positions_aggregated_sorted(self):
        reasons = analyze_reasons([make_variant(16190, "C", "T"), make_variant(16185, "A", "G")])
        assert "16180-16193 region (16185, 16190)" in reasons

    def test_single_position_region(self):
        reasons = analyze_reasons([make_variant(16185, "A", "G")])
        assert "16180-16193 region (16185)" in reasons

    def test_outside_region_not_aggregated(self):
        reasons = analyze_reasons([make_variant(16179, "A", "G"), make_variant(16194, "C", "T")])
        assert not any(r.startswith("16180-16193 region") for r in reasons)

    def test_region_boundaries_inclusive(self):
        assert any(r.startswith("16180-16193 region") for r in analyze_reasons([make_variant(16180, "A", "G")]))
        assert any(r.startswith("16180-16193 region") for r in analyze_reasons([make_variant(16193, "A", "G")]))


# ---------------------------------------------------------------------------
# T. VariantFlag.key() identity
# ---------------------------------------------------------------------------


class TestVariantFlagKey:
    def test_key_format(self):
        assert VariantFlag(pos="309.1", ref="-", seq="C").key() == "309.1|-|C"

    def test_key_normalizes_integer_pos(self):
        assert VariantFlag(pos=309, ref="A", seq="T").key() == "309|A|T"


# ---------------------------------------------------------------------------
# U. Special-position indels fully skipped (no flags at all)
# ---------------------------------------------------------------------------


class TestSpecialPositionIndelSkip:
    """Special-position insertions are NOT skipped: only the expected C
    insertion at 309/573 is suppressed, and 455.1/463.1 insertions are
    flagged. Special-position *deletions* (573/455/463) are skipped; the
    309 base deletion is still flagged. Non-indel SNPs at special positions
    are flagged via the specific-pattern checks (e.g. Has 309T)."""

    def test_309_1Y_insertion_flagged_with_heteroplasmy(self):
        # 309.1Y: non-C insertion at a special position -> insertion flag +
        # heteroplasmy (Y is an IUPAC ambiguity code).
        result = flag_variants([make_variant("309.1", "-", "Y")])
        assert "309.1|-|Y" in result
        assert "Insertion Y at 309.1" in result["309.1|-|Y"]
        assert "Heteroplasmy at 309.1" in result["309.1|-|Y"]

    def test_573_1T_insertion_flagged(self):
        result = flag_variants([make_variant("573.1", "-", "T")])
        assert "573.1|-|T" in result
        assert "Insertion T at 573.1" in result["573.1|-|T"]

    def test_455_1_insertion_flagged(self):
        # 455 is a special position but NOT in the C-exception set, so even
        # a C insertion is flagged.
        result = flag_variants([make_variant("455.1", "-", "C")])
        assert "455.1|-|C" in result
        assert "Insertion C at 455.1" in result["455.1|-|C"]

    def test_463_1_insertion_flagged(self):
        result = flag_variants([make_variant("463.1", "-", "C")])
        assert "463.1|-|C" in result
        assert "Insertion C at 463.1" in result["463.1|-|C"]

    def test_573_base_deletion_skipped(self):
        # A special-position *deletion* (seq == "-") is skipped, not flagged.
        reasons = analyze_reasons([make_variant(573, "A", "-")])
        assert "Deletion at 573" not in reasons

    def test_455_base_deletion_skipped(self):
        reasons = analyze_reasons([make_variant(455, "A", "-")])
        assert "Deletion at 455" not in reasons

    def test_309_base_deletion_flagged(self):
        # The 309 base deletion (Tracy 309C>T) is a real event and is flagged.
        reasons = analyze_reasons([make_variant(309, "A", "-")])
        assert "Deletion at 309" in reasons

    def test_309_snp_still_flagged(self):
        # Non-indel SNP at base 309 (A->T) is still flagged as Has 309T.
        reasons = analyze_reasons([make_variant(309, "A", "T")])
        assert "Has 309T variant" in reasons
