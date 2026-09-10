"""Focused regression tests for HCLS/TNLS variant comparison."""

from typing import Any

import pytest

from src.modules.TNLS.verification.compare_1n_class import CompareManager
from src.modules.TNLS.verification.utils import compare_variants


def _sample(
    *,
    snps: list[dict[str, Any]] | None = None,
    insertions: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    return {
        "intervals": {"HV1": [[1, 400]]},
        "variants": {
            "snps": snps or [],
            "insertions": insertions or [],
            "deletions": [],
        },
    }


def test_hcls_n_and_equivalent_insertion_position_do_not_reject_pair() -> None:
    hcls = _sample(
        snps=[
            {"pos": 100, "ref": "A", "seq": "N"},
            {"pos": 200, "ref": "A", "seq": "R"},
        ],
        insertions=[{"pos": "315.1", "ref": "-", "seq": "C"}],
    )
    tnls = _sample(
        snps=[{"pos": 200, "ref": "A", "seq": "G"}],
        insertions=[{"pos": 315.1, "ref": "-", "seq": "C"}],
    )
    rules = {
        "minimum": 150,
        "region": ["HV1"],
        "conclusion_threshold": 0,
        "exclusion_threshold": 2,
    }
    manager = CompareManager(hcls, {"TNLS-1": tnls}, rules)

    stage2 = manager.s2_filter_variants_based(manager.s1_get_overlap_list())
    assert stage2.loc[0, "no_mismatch_pos"] == 0

    result = manager.run()
    assert result["target_sample"].tolist() == ["TNLS-1"]
    assert result.loc[0, "result"] == "CANNOT_EXCLUDE"
    assert result.loc[0, "variant_match"] == "200R 315.1C"
    assert result.loc[0, "variant_mismatch"] == "None"


@pytest.mark.parametrize("target", [{}, {"snp_100": {"ref": "A", "seq": "C"}}])
def test_hcls_n_is_ignored_with_or_without_target_call(target: dict[str, dict[str, str]]) -> None:
    hcls = {"snp_100": {"ref": "A", "seq": "N"}}

    assert compare_variants(hcls, target, [[1, 400]]) == (0, 0, {}, {})


def test_iupac_is_used_only_for_snps() -> None:
    hcls = {
        "snp_100": {"ref": "A", "seq": "R"},
        "ins_200.1": {"ref": "-", "seq": "R"},
    }
    tnls = {
        "snp_100": {"ref": "A", "seq": "M"},
        "ins_200.1": {"ref": "-", "seq": "A"},
    }

    matches, mismatches, match_details, mismatch_details = compare_variants(hcls, tnls, [[1, 400]])

    assert (matches, mismatches) == (1, 1)
    assert set(match_details) == {"snp_100"}
    assert set(mismatch_details) == {"ins_200.1"}
