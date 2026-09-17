"""Regression tests for final FIS-to-Sanger report enrichment."""

from src.core.models import Variant
from src.modules.NGS.enrich_fis_sanger_comparison import _fis_correct, _profile


def test_profile_uses_raw_calls_when_corrected_profile_is_dash():
    assert _profile({"FIS_Nomenclature correction genotype": "-", "FIS_Variants": "73G 315.1C"}) == "73G 315.1C"


def test_profile_uses_raw_calls_when_corrected_profile_is_slash():
    assert _profile({"FIS_Nomenclature correction genotype": "/", "FIS_Variants": "73G 315.1C"}) == "73G 315.1C"


def test_fis_correct_accepts_reference_covering_iupac_unique_call():
    row = {"Concordant": "N", "Variants Unique (fis)": "146Y", "Variants Unique (sanger)": ""}

    assert _fis_correct(row, [Variant(pos=146, ref="T", seq="Y")]) == "Y"


def test_fis_correct_accepts_16193_indel_unique_call():
    row = {"Concordant": "N", "Variants Unique (fis)": "16193.1C", "Variants Unique (sanger)": ""}

    assert _fis_correct(row, [Variant(pos="16193.1", ref="-", seq="C")]) == "Y"


def test_fis_correct_rejects_remaining_or_sanger_unique_calls():
    row = {"Concordant": "N", "Variants Unique (fis)": "146Y 200G", "Variants Unique (sanger)": ""}
    assert _fis_correct(row, [Variant(pos=146, ref="T", seq="Y"), Variant(pos=200, ref="A", seq="G")]) == "N"

    sanger_unique = {"Concordant": "N", "Variants Unique (fis)": "146Y", "Variants Unique (sanger)": "200G"}
    assert _fis_correct(sanger_unique, [Variant(pos=146, ref="T", seq="Y")]) == "N"


def test_fis_correct_preserves_core_concordance_and_na():
    assert _fis_correct({"Concordant": "Y"}, []) == "Y"
    assert _fis_correct({"Concordant": "N/A"}, []) == "N/A"
