"""Tracy integration tests for shared polyC policy decisions."""

import pytest

from src.tools.tracy.transforms import _apply_polyc_insertion, _apply_polyc_primer_filters
from src.tools.tracy.utils import PrimerTypeInfo, detect_variant_conditions


def _primers(
    *,
    hv1f: bool = False,
    hv1r: bool = False,
    hv2f: bool = False,
    hv2r: bool = False,
) -> PrimerTypeInfo:
    return PrimerTypeInfo(
        is_forward=hv1f or hv2f,
        is_reverse=hv1r or hv2r,
        is_hv2f=hv2f,
        is_hv2f_hv3f=hv2f,
        is_hv2r_hv3r=hv2r,
        is_hv1f=hv1f,
        is_hv1r=hv1r,
    )


def test_variant_conditions_preserve_hv1_polyc_trigger_compatibility():
    conditions = detect_variant_conditions([{"pos": 16189, "ref": "A", "seq": "C"}])

    assert conditions["has_16189_T_C"] is True


def test_raw_hv2_polyc_insertion_is_replaced_by_canonical_insertion():
    variant = {"pos": "315.1", "ref": "-", "seq": "C"}

    _apply_polyc_insertion(variant)

    assert variant["remove"] is True
    assert variant["reason"] == "HV2 polyC insertion normalized to canonical insertion"


@pytest.mark.parametrize(
    ("position", "primers", "polyc_created", "expected_reason"),
    [
        (
            16190,
            _primers(hv1f=True),
            True,
            "HV1 polyC region variant with 16189T>C present (HV1F primer - forward sequencing artifact)",
        ),
        (
            16188,
            _primers(hv1r=True),
            True,
            "HV1 polyC region variant with 16189T>C present (HV1R primer - reverse sequencing artifact)",
        ),
        (16189, _primers(hv1f=True), True, None),
        (16188, _primers(hv1f=True), True, None),
        (16190, _primers(hv1r=True), True, None),
        (16190, _primers(), True, None),
        (16190, _primers(hv1f=True), False, None),
    ],
)
def test_tracy_applies_shared_hv1_directional_policy(
    position,
    primers,
    polyc_created,
    expected_reason,
):
    variant = {"pos": position, "ref": "A", "seq": "G"}

    _apply_polyc_primer_filters(variant, primers, {"has_16189_T_C": polyc_created})

    assert variant.get("reason") == expected_reason
    assert variant.get("remove", False) is (expected_reason is not None)


@pytest.mark.parametrize(
    ("position", "primers", "expected_reason"),
    [
        (303, _primers(hv2f=True), None),
        (304, _primers(hv2f=True), "HV2 directional polyC variant (forward trace)"),
        (341, _primers(hv2f=True), "HV2 directional polyC variant (forward trace)"),
        (1, _primers(hv2r=True), "HV2 directional polyC variant (reverse trace)"),
        (310, _primers(hv2r=True), "HV2 directional polyC variant (reverse trace)"),
        (315, _primers(hv2r=True), "HV2 directional polyC variant (reverse trace)"),
        (316, _primers(hv2r=True), None),
        (304, _primers(), None),
    ],
)
def test_tracy_applies_shared_hv2_directional_coverage_policy(position, primers, expected_reason):
    variant = {"pos": position, "ref": "A", "seq": "G"}

    _apply_polyc_primer_filters(variant, primers, {})

    assert variant.get("reason") == expected_reason
    assert variant.get("remove", False) is (expected_reason is not None)
