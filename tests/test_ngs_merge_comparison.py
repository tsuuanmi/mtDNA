"""Regression tests for NGS comparison interval filtering."""

from src.modules.NGS.merge_comparison import filter_variants_by_intervals


def test_boundary_insertion_is_excluded_from_sanger_coverage():
    variants = {"16193C", "16193.1C"}

    assert filter_variants_by_intervals(variants, [(16024, 16193)]) == {"16193C"}


def test_insertion_is_retained_when_sanger_coverage_extends_past_anchor():
    variants = {"16193.1C"}

    assert filter_variants_by_intervals(variants, [(16024, 16194)]) == variants
