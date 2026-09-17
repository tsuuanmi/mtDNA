"""Tests for indexed TNLS mtDNA-profile family grouping."""

from typing import Any

from src.modules.TNLS.verification.family_profiles import group_families


def _sample(
    intervals: list[list[int]],
    *,
    snps: list[tuple[int | str, str]] | None = None,
    insertions: list[tuple[int | str, str]] | None = None,
) -> dict[str, Any]:
    return {
        "intervals": {"HV1": intervals},
        "variants": {
            "snps": [{"pos": position, "ref": "A", "seq": sequence} for position, sequence in snps or []],
            "insertions": [{"pos": position, "ref": "-", "seq": sequence} for position, sequence in insertions or []],
            "deletions": [],
        },
    }


def test_groups_exact_duplicates_and_iupac_transitively() -> None:
    samples = {
        "a": _sample([[1, 200]], snps=[(10, "A")]),
        "a_duplicate": _sample([[1, 200]], snps=[(10, "A")]),
        "heteroplasmy": _sample([[1, 200]], snps=[(10, "R")]),
        "g": _sample([[1, 200]], snps=[(10, "G")]),
        "different_position": _sample([[1, 200]], snps=[(11, "A")]),
    }

    assert group_families(samples, ["HV1"], 150) == [
        ["a", "a_duplicate", "g", "heteroplasmy"],
        ["different_position"],
    ]


def test_groups_profiles_by_their_common_coverage() -> None:
    samples = {
        "full_c": _sample([[1, 200]], snps=[(10, "A"), (190, "C")]),
        "full_g": _sample([[1, 200]], snps=[(10, "A"), (190, "G")]),
        "partial": _sample([[1, 180]], snps=[(10, "A")]),
    }

    assert group_families(samples, ["HV1"], 150) == [["full_c", "full_g", "partial"]]


def test_keeps_nonmatching_and_short_coverage_samples_separate() -> None:
    samples = {
        "short_a": _sample([[1, 100]], snps=[(10, "A")]),
        "short_duplicate": _sample([[1, 100]], snps=[(10, "A")]),
        "snp": _sample([[1, 200]], snps=[(10, "A")]),
        "insertion": _sample([[1, 200]], insertions=[(10, "A")]),
        "unknown_a": _sample([[1, 200]], snps=[(10, "?")]),
        "unknown_duplicate": _sample([[1, 200]], snps=[(10, "?")]),
    }

    assert group_families(samples, ["HV1"], 150) == [
        ["insertion"],
        ["short_a"],
        ["short_duplicate"],
        ["snp"],
        ["unknown_a"],
        ["unknown_duplicate"],
    ]
