"""Focused tests for the canonical FIS profile adapter."""

from pathlib import Path

from src.core.flagging import SampleFlagger
from src.tools.fis.etl import load_reference_bases, profile_to_variants
from src.tools.fis.utils import _integer_or_none, select_profile

REFERENCE = Path("ref/rCRS.fasta")


def test_corrected_profile_wins_over_raw_lowercase_call():
    profile, source = select_profile("16193.1C", "16193.1c")
    variants = profile_to_variants(profile, load_reference_bases(REFERENCE))

    _, reasons, _ = SampleFlagger(variants).analyze()

    assert source == "corrected"
    assert "Lowercase variant detected" not in reasons


def test_missing_corrected_profile_falls_back_to_raw():
    assert select_profile("-", "73G 315.1C") == ("73G 315.1C", "raw_fallback")
    assert select_profile("/", "73G 315.1C") == ("73G 315.1C", "raw_fallback")


def test_reference_backed_523_524_deletion_pair_is_suppressed():
    variants = profile_to_variants("523DEL 524DEL 315.1C", load_reference_bases(REFERENCE))
    _, reasons, _ = SampleFlagger(variants).analyze()

    assert "Deletion at 523" not in reasons
    assert "Deletion at 524" not in reasons


def test_raw_fallback_preserves_lowercase_for_flagging():
    variants = profile_to_variants("73g", load_reference_bases(REFERENCE), normalize_alleles=False)
    _, reasons, _ = SampleFlagger(variants).analyze()

    assert "Lowercase variant detected" in reasons


def test_raw_integer_column_accepts_excel_float():
    assert _integer_or_none(10.0) == 10
