"""Regression tests for shared directional polyC policy helpers."""

import pytest

from src.core.polyc import (
    directional_hv1_polyc_suppression_reason,
    directional_hv2_polyc_exclusion_reason,
    is_hv1_polyc_created,
)


@pytest.mark.parametrize(
    ("position", "seq", "expected"),
    [
        (16189, "C", True),
        ("16189.1", "C", True),
        (16189, "T", False),
        (16188, "C", False),
    ],
)
def test_hv1_polyc_trigger_uses_called_base_and_position_anchor(position, seq, expected):
    assert is_hv1_polyc_created(position, seq) is expected


@pytest.mark.parametrize(
    ("position", "direction", "expected"),
    [
        (
            16190,
            "forward",
            "HV1 polyC region variant with 16189T>C present (HV1F primer - forward sequencing artifact)",
        ),
        (
            16188,
            "reverse",
            "HV1 polyC region variant with 16189T>C present (HV1R primer - reverse sequencing artifact)",
        ),
        (16189, "forward", None),
        (16189, "reverse", None),
        (16188, "forward", None),
        (16190, "reverse", None),
        (16190, None, None),
    ],
)
def test_hv1_directional_suppression_uses_strict_boundaries(position, direction, expected):
    assert directional_hv1_polyc_suppression_reason(position, direction, polyc_created=True) == expected


def test_hv1_directional_suppression_requires_polyc_trigger():
    assert directional_hv1_polyc_suppression_reason(16190, "forward", polyc_created=False) is None


@pytest.mark.parametrize(
    ("position", "direction", "expected"),
    [
        (303, "forward", False),
        (304, "forward", True),
        ("16569.1", "forward", True),
        (1, "reverse", True),
        ("315.1", "reverse", True),
        (315, "reverse", True),
        (316, "reverse", False),
        (16569, "reverse", False),
    ],
)
def test_hv2_directional_exclusion_has_no_fixed_trace_endpoint(position, direction, expected):
    reason = directional_hv2_polyc_exclusion_reason(position, direction)
    assert (reason is not None) is expected
