"""Comprehensive tests for src.core.variants.

Covers every public function and the key private helpers:
position helpers, variant sequence generation, I/O, classification,
statistics, formatting, region classification, and model conversion.
"""

from pathlib import Path

import pytest

from src.config import get_settings
from src.core.models import Position, Sample, Sequence, Tool, Variant, validate_position
from src.core.sample import (
    _apply_variants_to_dict,
    _extract_region,
    _load_reference_sequence,
    _mark_non_region_as_n,
    generate_sequence,
)
from src.core.variants import (
    format_variants_list,
    format_variants_simplified,
    get_hv_region_for_position,
    get_overlap,
    is_position_in_intervals,
    is_special_position,
    iupac_bases_compatible,
    normalize_position,
    parse_variant_position_allele,
    pos_base,
    pos_sort_key,
    split_variant_types,
    statistic_variants,
    variant_to_dict,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ref_path() -> str:
    """Path to the bundled rCRS reference FASTA (16569 bp)."""
    return str(Path(__file__).parent.parent / "ref" / "rCRS.fasta")


@pytest.fixture
def ref_dict(ref_path):
    """Reference sequence loaded into a 1-based position->base dict."""
    return _load_reference_sequence(ref_path)


@pytest.fixture
def synthetic_ref():
    """Small synthetic reference dict (positions 1-6)."""
    return {1: "A", 2: "C", 3: "G", 4: "T", 5: "A", 6: "C"}


# ---------------------------------------------------------------------------
# Re-exports
# ---------------------------------------------------------------------------


class TestReexports:
    """Position and validate_position are re-exported from src.core.variants."""

    def test_position_alias_is_int_or_str(self):
        assert Position == int | str

    def test_validate_position_callable(self):
        assert callable(validate_position)
        assert validate_position(73) == 73
        assert validate_position("309.1") == "309.1"

    def test_validate_position_rejects_float(self):
        with pytest.raises(TypeError):
            validate_position(309.1)


# ---------------------------------------------------------------------------
# pos_base
# ---------------------------------------------------------------------------


class TestPosBase:
    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            (0, 0),
            (5, 5),
            (-10, -10),
            (309, 309),
            (5.0, 5),
            (309.1, 309),
            (309.9, 309),
            ("309", 309),
            ("309.1", 309),
            ("217.10", 217),
        ],
    )
    def test_pos_base(self, pos, expected):
        assert pos_base(pos) == expected


# ---------------------------------------------------------------------------
# IUPAC compatibility
# ---------------------------------------------------------------------------


class TestIupacBasesCompatible:
    @pytest.mark.parametrize(
        ("left", "right", "expected"),
        [
            ("A", "A", True),
            ("R", "G", True),
            ("R", "M", True),
            ("R", "Y", False),
            ("N", "N", True),
            ("n", "c", True),
            ("?", "?", False),
        ],
    )
    def test_uses_possible_base_intersection(self, left, right, expected):
        assert iupac_bases_compatible(left, right) is expected


# ---------------------------------------------------------------------------
# normalize_position
# ---------------------------------------------------------------------------


class TestNormalizePosition:
    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            # ints / integer-valued -> int
            (309, 309),
            (309.0, 309),
            (0, 0),
            (-1, -1),
            # decimal floats -> str (precision preserved from str())
            (309.1, "309.1"),
            # strings
            ("217.10", "217.10"),
            ("1,234.10", "1234.10"),
            ("309", 309),
            (" 309 ", 309),
            ("309.0", 309),
        ],
    )
    def test_canonical(self, pos, expected):
        assert normalize_position(pos) == expected

    def test_returns_int_type_for_bases(self):
        assert isinstance(normalize_position(309), int)
        assert isinstance(normalize_position("309"), int)
        assert isinstance(normalize_position(309.0), int)

    def test_returns_str_type_for_insertions(self):
        assert isinstance(normalize_position("309.1"), str)
        assert isinstance(normalize_position(309.1), str)

    def test_non_numeric_string_returned_cleaned(self):
        assert normalize_position(" abc ") == "abc"
        assert normalize_position("abc") == "abc"

    def test_strips_commas_and_whitespace(self):
        assert normalize_position("1,234") == 1234
        assert normalize_position("  309  ") == 309

    @pytest.mark.parametrize("bad", [None, [309], (309,), object()])
    def test_unsupported_type_raises_type_error(self, bad):
        with pytest.raises(TypeError, match="Cannot normalize position"):
            normalize_position(bad)

    def test_float_decimal_not_recoverable_as_ten(self):
        # float(217.10) == 217.1 -> ".10" cannot be recovered from a float
        assert normalize_position(217.10) == "217.1"


# ---------------------------------------------------------------------------
# pos_sort_key
# ---------------------------------------------------------------------------


class TestPosSortKey:
    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            (309, (309, -1)),
            ("309.1", (309, 1)),
            ("217.10", (217, 10)),
            ("217.9", (217, 9)),
            ("1,234.5", (1234, 5)),
            ("  309  ", (309, -1)),
            ("309", (309, -1)),
        ],
    )
    def test_key(self, pos, expected):
        assert pos_sort_key(pos) == expected

    def test_base_sorts_before_insertions(self):
        ordered = sorted([309, "309.3", "309.1", 308, "309.2"], key=pos_sort_key)
        assert ordered == [308, 309, "309.1", "309.2", "309.3"]

    def test_multi_digit_insertion_index_ordering(self):
        # 217.10 must sort AFTER 217.9 (numeric, not lexicographic)
        ordered = sorted(["217.10", "217.9", "217.1"], key=pos_sort_key)
        assert ordered == ["217.1", "217.9", "217.10"]

    def test_non_numeric_insertion_falls_back(self):
        # ".x" cannot be parsed as int -> insertion index -1 (base-like sort)
        assert pos_sort_key("309.x") == (309, -1)

    def test_trailing_dot_treated_as_base(self):
        assert pos_sort_key("309.") == (309, -1)

    def test_non_numeric_no_dot_raises(self):
        # pos_base raises on a non-numeric string with no decimal point before the
        # final defensive return is reached.
        with pytest.raises(ValueError, match="invalid literal"):
            pos_sort_key("abc")


# ---------------------------------------------------------------------------
# _load_reference_sequence
# ---------------------------------------------------------------------------


class TestLoadReferenceSequence:
    def test_loads_real_reference(self, ref_path):
        ref = _load_reference_sequence(ref_path)
        assert len(ref) == 16569
        assert ref[1] == "G"  # rCRS position 1
        assert ref[73] == "A"
        assert ref[16569] == "G"
        # 1-based, no position 0
        assert 0 not in ref
        assert all(isinstance(k, int) for k in ref)

    def test_empty_fasta_returns_empty_dict(self, tmp_path):
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        assert _load_reference_sequence(str(empty)) == {}

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            _load_reference_sequence("/no/such/reference.fasta")


# ---------------------------------------------------------------------------
# _apply_variants_to_dict
# ---------------------------------------------------------------------------


class TestApplyVariantsToDict:
    def test_snp_replaces_base(self, synthetic_ref):
        result = _apply_variants_to_dict(synthetic_ref, [{"pos": 3, "ref": "G", "seq": "A"}])
        assert result[3] == "A"
        # original dict untouched
        assert synthetic_ref[3] == "G"

    def test_insertion_uses_normalized_position_key(self, synthetic_ref):
        result = _apply_variants_to_dict(synthetic_ref, [{"pos": 3.1, "ref": "-", "seq": "C"}])
        # float 3.1 normalizes to str "3.1"
        assert result["3.1"] == "C"
        assert result[3] == "G"

    def test_deletion_marked_with_dash(self, synthetic_ref):
        result = _apply_variants_to_dict(synthetic_ref, [{"pos": 4, "ref": "T", "seq": "-"}])
        assert result[4] == "-"

    def test_accepts_dict_format_input(self, synthetic_ref):
        variants = {3: {"ref": "G", "seq": "A"}}
        result = _apply_variants_to_dict(synthetic_ref, variants)
        assert result[3] == "A"

    def test_accepts_list_format_input(self, synthetic_ref):
        result = _apply_variants_to_dict(synthetic_ref, [{"pos": 3, "ref": "G", "seq": "A"}])
        assert result[3] == "A"

    def test_does_not_mutate_input_ref_dict(self, synthetic_ref):
        _apply_variants_to_dict(synthetic_ref, [{"pos": 1, "ref": "A", "seq": "T"}])
        assert synthetic_ref[1] == "A"

    def test_multi_base_insertion_applied_as_string(self, synthetic_ref):
        result = _apply_variants_to_dict(synthetic_ref, [{"pos": 2.1, "ref": "-", "seq": "CCCC"}])
        assert result["2.1"] == "CCCC"

    def test_ref_mismatch_still_applies_variant(self, synthetic_ref):
        # ref base "X" does not match the reference genome base "G" at pos 3;
        # a warning is logged but the variant is still applied.
        result = _apply_variants_to_dict(synthetic_ref, [{"pos": 3, "ref": "X", "seq": "A"}])
        assert result[3] == "A"


# ---------------------------------------------------------------------------
# _mark_non_region_as_n
# ---------------------------------------------------------------------------


class TestMarkNonRegionAsN:
    def test_marks_outside_positions_as_n(self):
        seq_dict = {1: "A", 2: "C", 3: "G", 4: "T", 5: "A"}
        result = _mark_non_region_as_n(seq_dict, [[2, 3]])
        assert result[1] == "N"
        assert result[2] == "C"
        assert result[3] == "G"
        assert result[4] == "N"
        assert result[5] == "N"

    def test_insertion_keys_not_touched(self):
        seq_dict = {1: "A", "2.1": "C", 10: "G"}
        result = _mark_non_region_as_n(seq_dict, [[1, 1]])
        # insertion stays
        assert result["2.1"] == "C"
        # out-of-region base marked N
        assert result[10] == "N"

    def test_empty_regions_returns_unchanged(self):
        seq_dict = {1: "A", 2: "C"}
        assert _mark_non_region_as_n(seq_dict, []) is seq_dict
        assert _mark_non_region_as_n(seq_dict, []) == {1: "A", 2: "C"}

    def test_does_not_mutate_input(self):
        seq_dict = {1: "A", 10: "G"}
        _mark_non_region_as_n(seq_dict, [[1, 1]])
        assert seq_dict[10] == "G"

    def test_multiple_regions(self):
        seq_dict = {1: "A", 5: "G", 10: "T", 20: "C"}
        result = _mark_non_region_as_n(seq_dict, [[1, 2], [10, 12]])
        assert result[1] == "A"
        assert result[5] == "N"
        assert result[10] == "T"
        assert result[20] == "N"


# ---------------------------------------------------------------------------
# _extract_region
# ---------------------------------------------------------------------------


class TestExtractRegion:
    def test_extracts_in_order(self):
        seq_dict = {1: "A", 2: "C", 3: "G", 4: "T"}
        assert _extract_region(seq_dict, 1, 4) == "ACGT"

    def test_skips_deletion_markers(self):
        seq_dict = {1: "A", 2: "-", 3: "G"}
        assert _extract_region(seq_dict, 1, 3) == "AG"

    def test_expands_multi_base_insertions(self):
        seq_dict = {2: "C", "2.1": "XYZ", 3: "G"}
        # insertion "2.1" sorts between 2 and 3, expanded as individual bases
        assert _extract_region(seq_dict, 2, 3) == "CXYZG"

    def test_insertion_sorting_by_index(self):
        seq_dict = {"2.1": "X", "2.10": "Y", "2.2": "Z"}
        # numeric insertion ordering: .1 < .2 < .10 -> X(2.1), Z(2.2), Y(2.10)
        assert _extract_region(seq_dict, 2, 2) == "XZY"

    def test_excludes_out_of_range_keys(self):
        seq_dict = {1: "A", 2: "C", 3: "G"}
        assert _extract_region(seq_dict, 2, 2) == "C"

    def test_empty_when_no_keys_in_range(self):
        assert _extract_region({1: "A"}, 100, 200) == ""


# ---------------------------------------------------------------------------
# generate_sequence
# ---------------------------------------------------------------------------


class TestGenerateSequence:
    def _make_sample(self, variants=None, intervals=None) -> Sample:
        return Sample(
            sample_id="S1",
            source_tool=Tool.UNIFIED,
            variants=variants or [],
            intervals=intervals,
            hv1=None,
            hv2=None,
            hv3=None,
            no_snps=None,
            no_ins=None,
            no_dels=None,
            no_identicals=None,
            no_unread=None,
            batch_id=None,
            information=None,
        )

    def test_returns_sequence_with_regions(self, ref_path):
        sample = self._make_sample(
            variants=[Variant(pos=73, ref="A", seq="G", peaks=None)],
            intervals={"HV2": [[73, 340]]},
        )
        seq = generate_sequence(sample, ref_path)
        assert isinstance(seq, Sequence)
        assert len(seq.reference_seq) == 16569
        assert seq.consensus_seq[72] == "G"
        assert set(seq.hv_seqs) == {"HV1", "HV2", "HV3"}
        assert set(seq.hv_seq_refs) == {"HV1", "HV2", "HV3"}

    def test_empty_variants_consensus_equals_reference(self, ref_path):
        sample = self._make_sample(variants=[], intervals={"HV2": [[73, 340]]})
        seq = generate_sequence(sample, ref_path)
        # No variants and no deletions -> reference and consensus are identical
        # (both get the same N-marking outside the sequenced region).
        assert seq.reference_seq == seq.consensus_seq

    def test_no_intervals(self, ref_path):
        sample = self._make_sample(variants=[], intervals=None)
        seq = generate_sequence(sample, ref_path)
        assert isinstance(seq, Sequence)
        assert len(seq.reference_seq) == 16569

    def test_failed_load_returns_empty_sequence(self, tmp_path):
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        sample = self._make_sample(variants=[], intervals={"HV2": [[73, 340]]})
        seq = generate_sequence(sample, str(empty))
        assert seq.reference_seq == ""
        assert seq.consensus_seq == ""
        assert seq.hv_seqs == {}

    def test_default_ref_path_from_config(self):
        # No ref_path -> uses config default (ref/rCRS.fasta).
        sample = self._make_sample(variants=[], intervals={"HV2": [[73, 340]]})
        seq = generate_sequence(sample)
        assert isinstance(seq, Sequence)
        assert len(seq.reference_seq) == 16569

    def test_nested_region_format(self, monkeypatch):
        # REGIONS may be nested ([[start, end], ...]); ensure that format is handled.
        nested = {
            "HV1": [[16024, 16365]],
            "HV2": [[73, 340]],
            "HV3": [[438, 576]],
        }
        monkeypatch.setattr(get_settings().regions, "REGIONS", nested)
        sample = self._make_sample(variants=[], intervals={"HV2": [[73, 340]]})
        seq = generate_sequence(sample)
        assert seq.hv_seq_refs["HV1"]["start"] == 16024
        assert seq.hv_seq_refs["HV1"]["end"] == 16365
        assert seq.hv_seqs["HV2"]["start"] == 73


# ---------------------------------------------------------------------------
# split_variant_types
# ---------------------------------------------------------------------------


class TestSplitVariantTypes:
    def test_splits_and_normalizes(self):
        variants = [
            {"pos": "73", "ref": "A", "seq": "G"},  # snp
            {"pos": "309.1", "ref": "-", "seq": "C"},  # insertion
            {"pos": "333", "ref": "T", "seq": "-"},  # deletion
        ]
        snps, insertions, deletions = split_variant_types(variants)
        assert len(snps) == 1
        assert len(insertions) == 1
        assert len(deletions) == 1
        # positions normalized in returned dicts
        assert snps[0]["pos"] == 73
        assert isinstance(snps[0]["pos"], int)
        assert insertions[0]["pos"] == "309.1"
        assert deletions[0]["pos"] == 333

    def test_sorted_by_pos_sort_key(self):
        variants = [
            {"pos": 333, "ref": "T", "seq": "-"},
            {"pos": 73, "ref": "A", "seq": "G"},
            {"pos": 100, "ref": "A", "seq": "G"},
        ]
        snps, _ins, _dels = split_variant_types(variants)
        assert [v["pos"] for v in snps] == [73, 100]

    def test_empty_input(self):
        snps, insertions, deletions = split_variant_types([])
        assert snps == insertions == deletions == []


# ---------------------------------------------------------------------------
# statistic_variants
# ---------------------------------------------------------------------------


class TestStatisticVariants:
    def test_counts_snps_ins_dels(self):
        variants = {
            73: {"ref": "A", "seq": "G"},  # snp
            309: {"ref": "-", "seq": "C"},  # insertion
            333: {"ref": "T", "seq": "-"},  # deletion
            100: {"ref": "A", "seq": "A"},  # identical (not snp/ins/del)
        }
        hv_seq_refs = {"HV1": {"seq": ["A", "C", "G"]}}
        hv_seqs = {"HV1": {"seq": ["A", "C", "G"]}}
        snps, ins, dels, identicals, unread = statistic_variants(variants, hv_seq_refs, hv_seqs)
        assert snps == 1
        assert ins == 1
        assert dels == 1
        assert identicals == 3
        assert unread == 0

    def test_counts_unread_n_positions(self):
        variants = {73: {"ref": "A", "seq": "G"}}
        hv_seq_refs = {"HV1": {"seq": ["A", "C", "G"]}}
        hv_seqs = {"HV1": {"seq": ["A", "N", "G"]}}
        _snps, _ins, _dels, identicals, unread = statistic_variants(variants, hv_seq_refs, hv_seqs)
        assert unread == 1
        assert identicals == 2  # idx 0 and 2 match

    def test_no_variants(self):
        hv_seq_refs = {"HV1": {"seq": ["A"]}}
        hv_seqs = {"HV1": {"seq": ["A"]}}
        snps, ins, dels, identicals, unread = statistic_variants({}, hv_seq_refs, hv_seqs)
        assert (snps, ins, dels) == (0, 0, 0)
        assert identicals == 1
        assert unread == 0


# ---------------------------------------------------------------------------
# format_variants_simplified
# ---------------------------------------------------------------------------


class TestFormatVariantsSimplified:
    def test_empty_returns_none_string(self):
        assert format_variants_simplified([]) == "None"

    def test_snps_sorted(self):
        variants = [
            {"pos": "263", "ref": "A", "seq": "G"},
            {"pos": "73", "ref": "A", "seq": "G"},
        ]
        assert format_variants_simplified(variants) == "73G 263G"

    def test_insertion_format(self):
        variants = [{"pos": "315.1", "ref": "-", "seq": "C"}]
        assert format_variants_simplified(variants) == "315.1C"

    def test_deletion_format(self):
        variants = [{"pos": "524", "ref": "A", "seq": "-"}]
        assert format_variants_simplified(variants) == "524DEL"

    def test_mixed_sorted_together(self):
        variants = [
            {"pos": "524", "ref": "A", "seq": "-"},
            {"pos": "73", "ref": "A", "seq": "G"},
            {"pos": "315.1", "ref": "-", "seq": "C"},
        ]
        # sorted by pos_sort_key: 73 < 315.1(base 315) < 524
        assert format_variants_simplified(variants) == "73G 315.1C 524DEL"


# ---------------------------------------------------------------------------
# format_variants_list
# ---------------------------------------------------------------------------


class TestFormatVariantsList:
    def test_default_space_separator(self):
        assert format_variants_list(["73G", "263G", "315.1C"]) == "73G 263G 315.1C"

    def test_explicit_space_separator(self):
        assert format_variants_list(["73G", "263G"], separator="space") == "73G 263G"

    def test_comma_separator(self):
        assert format_variants_list(["73G", "263G"], separator="comma") == "73G, 263G"

    def test_custom_separator(self):
        assert format_variants_list(["73G", "263G"], separator=";") == "73G;263G"

    def test_empty_returns_default_empty_value(self):
        assert format_variants_list([]) == "None"

    def test_empty_custom_value(self):
        assert format_variants_list([], empty_value="-") == "-"


# ---------------------------------------------------------------------------
# parse_variant_position_allele
# ---------------------------------------------------------------------------


class TestParseVariantPositionAllele:
    @pytest.mark.parametrize(
        ("variant", "expected"),
        [
            ("152C", ("152", "C")),
            ("315.1C", ("315.1", "C")),
            ("524DEL", ("524", "DEL")),
            ("309.1DEL", ("309.1", "DEL")),
            ("73G", ("73", "G")),
        ],
    )
    def test_parse(self, variant, expected):
        assert parse_variant_position_allele(variant) == expected

    def test_no_alpha_returns_position_only(self):
        assert parse_variant_position_allele("12345") == ("12345", "")

    def test_pure_del(self):
        assert parse_variant_position_allele("DEL") == ("", "DEL")


# ---------------------------------------------------------------------------
# get_hv_region_for_position
# ---------------------------------------------------------------------------


class TestGetHvRegionForPosition:
    def test_default_regions_hv1(self):
        assert get_hv_region_for_position(16200) == "HV1"

    def test_default_regions_hv2(self):
        assert get_hv_region_for_position(100) == "HV2"

    def test_default_regions_hv3(self):
        assert get_hv_region_for_position(500) == "HV3"

    def test_outside_any_region_returns_none(self):
        assert get_hv_region_for_position(1000) is None

    def test_insertion_position_string(self):
        assert get_hv_region_for_position("309.1") == "HV2"

    def test_custom_flat_regions(self):
        regions = {"R1": [10, 20], "R2": [50, 60]}
        assert get_hv_region_for_position(15, regions) == "R1"
        assert get_hv_region_for_position(55, regions) == "R2"
        assert get_hv_region_for_position(30, regions) is None

    def test_custom_nested_regions(self):
        regions = {"R1": [[10, 20], [30, 40]]}
        assert get_hv_region_for_position(15, regions) == "R1"
        assert get_hv_region_for_position(35, regions) == "R1"
        assert get_hv_region_for_position(25, regions) is None


# ---------------------------------------------------------------------------
# get_overlap
# ---------------------------------------------------------------------------


class TestGetOverlap:
    def test_partial_overlap(self):
        assert get_overlap([[10, 20]], [[15, 25]]) == [[15, 20]]

    def test_no_overlap(self):
        assert get_overlap([[10, 20]], [[30, 40]]) == []

    def test_contained(self):
        assert get_overlap([[10, 30]], [[15, 25]]) == [[15, 25]]

    def test_multiple_intervals(self):
        a = [[10, 20], [30, 40]]
        b = [[15, 35]]
        assert get_overlap(a, b) == [[15, 20], [30, 35]]

    def test_empty_inputs(self):
        assert get_overlap([], [[1, 2]]) == []
        assert get_overlap([[1, 2]], []) == []


# ---------------------------------------------------------------------------
# is_position_in_intervals
# ---------------------------------------------------------------------------


class TestIsPositionInIntervals:
    def test_inside_interval(self):
        assert is_position_in_intervals("152", [(73, 340), (16024, 16365)]) is True

    def test_outside_intervals(self):
        assert is_position_in_intervals("500", [(73, 340), (16024, 16365)]) is False

    def test_int_position(self):
        assert is_position_in_intervals(309, [(73, 340)]) is True

    def test_insertion_string_uses_base_within_interval(self):
        assert is_position_in_intervals("309.1", [(73, 340)]) is True

    def test_base_at_interval_end_is_covered(self):
        assert is_position_in_intervals(16193, [(16024, 16193)]) is True

    def test_insertion_after_interval_end_is_not_covered(self):
        assert is_position_in_intervals("16193.1", [(16024, 16193)]) is False

    def test_insertion_before_interval_end_is_covered(self):
        assert is_position_in_intervals("16192.1", [(16024, 16193)]) is True

    def test_empty_intervals_returns_false(self):
        assert is_position_in_intervals(100, []) is False

    def test_unparseable_returns_false(self):
        # pos_base("abc") raises ValueError -> caught -> False
        assert is_position_in_intervals("abc", [(1, 2)]) is False


# ---------------------------------------------------------------------------
# is_special_position
# ---------------------------------------------------------------------------


class TestIsSpecialPosition:
    @pytest.mark.parametrize("pos", ["455", "463", "573", "309"])
    def test_special_without_ref_seq(self, pos):
        assert is_special_position(pos) is True

    @pytest.mark.parametrize("pos", ["100", "16024", "524"])
    def test_non_special(self, pos):
        assert is_special_position(pos) is False

    @pytest.mark.parametrize("pos", ["309.1", "573.1", "455.2", "463.10"])
    def test_special_insertion_variant(self, pos):
        assert is_special_position(pos) is True

    def test_snp_at_special_position_not_filtered(self):
        # ref/seq both provided and not an INDEL -> False
        assert is_special_position("309", ref="A", seq="C") is False

    def test_insertion_at_special_position_filtered(self):
        assert is_special_position("309", ref="-", seq="C") is True

    def test_deletion_at_special_position_filtered(self):
        assert is_special_position("309", ref="A", seq="-") is True

    def test_length_mismatch_filtered(self):
        assert is_special_position("309", ref="A", seq="CC") is True

    def test_only_one_of_ref_seq_provided(self):
        # ref provided, seq None -> not the INDEL branch -> True (special)
        assert is_special_position("309", ref="A") is True
        assert is_special_position("309", seq="C") is True


# ---------------------------------------------------------------------------
# variant_to_dict
# ---------------------------------------------------------------------------


class TestVariantToDict:
    def test_variant_model_to_dict(self):
        v = Variant(pos=73, ref="A", seq="G", files=["a.ab1"], quality=[20], peaks=None)
        result = variant_to_dict(v)
        assert result == {
            "pos": 73,
            "ref": "A",
            "seq": "G",
            "file": ["a.ab1"],
            "quality": [20],
        }

    def test_variant_model_with_peaks(self):
        v = Variant(pos=73, ref="A", seq="G", peaks=[[1, 2, 3, 4]])
        result = variant_to_dict(v)
        assert result is not None
        assert result["peaks"] == [[1, 2, 3, 4]]

    def test_variant_model_empty_files_quality(self):
        v = Variant(pos=73, ref="A", seq="G", peaks=None)
        result = variant_to_dict(v)
        assert result is not None
        assert result["file"] == []
        assert result["quality"] == []
        assert "peaks" not in result

    def test_dict_input_normalizes_pos(self):
        result = variant_to_dict({"pos": "73", "ref": "A", "seq": "G"})
        assert result is not None
        assert result["pos"] == 73
        assert isinstance(result["pos"], int)
        assert result["file"] == []
        assert result["quality"] == []

    def test_dict_input_position_key_alias(self):
        result = variant_to_dict({"position": 73, "ref": "A", "seq": "G"})
        assert result is not None
        assert result["pos"] == 73

    def test_dict_input_alt_mapped_to_seq(self):
        result = variant_to_dict({"pos": 73, "ref": "A", "alt": "G"})
        assert result is not None
        assert result["seq"] == "G"
        assert "alt" not in result

    def test_dict_input_with_peaks(self):
        result = variant_to_dict({"pos": 73, "ref": "A", "seq": "G", "peaks": [[1.0]]})
        assert result is not None
        assert result["peaks"] == [[1.0]]

    def test_dict_input_none_position_returns_none(self):
        assert variant_to_dict({"ref": "A", "seq": "G"}) is None

    def test_dict_input_explicit_none_pos_returns_none(self):
        assert variant_to_dict({"pos": None, "ref": "A", "seq": "G"}) is None

    def test_dict_input_insertion_position_preserved(self):
        result = variant_to_dict({"pos": "309.1", "ref": "-", "seq": "C"})
        assert result is not None
        assert result["pos"] == "309.1"
        assert result["ref"] == "-"
        assert result["seq"] == "C"
