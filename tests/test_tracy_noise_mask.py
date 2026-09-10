import pytest

from src.core.flagging import flag_variants
from src.core.models import Sample, Tool, Variant
from src.tools.tracy.noise_mask import apply_noise_mask
from src.tools.tracy.quality_control import NoiseRange
from src.tools.tracy.utils import variant_overlaps_range


def _range(start: int, end: int) -> NoiseRange:
    return NoiseRange(
        start=start,
        end=end,
        alignment_start=0,
        alignment_end=end - start,
        supporting_windows=2,
        median_snr=1.5,
        low_end_snr=1.0,
        median_purity=0.4,
    )


def _sample(variants: list[Variant]) -> Sample:
    return Sample(
        sample_id="S1",
        source_tool=Tool.TRACY,
        variants=variants,
        intervals={"S1_HV1F.json": [[100, 140]]},
        sample_flags=[],
        variant_flags=flag_variants(variants),
    )


@pytest.mark.parametrize(
    ("position", "ref", "start", "end", "expected"),
    [
        (100, "A", 101, 105, False),
        (100, "A", 100, 105, True),
        ("315.1", "-", 315, 315, True),
        ("315.1", "-", 316, 320, False),
        (120, "AC", 120, 120, True),
        (120, "AC", 121, 125, True),
        (120, "AC", 122, 125, False),
    ],
)
def test_variant_range_overlap(position, ref, start, end, expected):
    assert variant_overlaps_range(position, ref, start, end) is expected


def test_apply_noise_mask_excludes_overlaps_and_recomputes_flags():
    variants = [
        Variant(pos=100, ref="A", seq="G", files=["S1_HV1F.json"]),
        Variant(pos="110.1", ref="-", seq="C", files=["S1_HV1F.json"]),
        Variant(pos=120, ref="AC", seq="-", files=["S1_HV1F.json"]),
        Variant(pos=130, ref="T", seq="C", files=["S1_HV1F.json"]),
    ]
    sample = _sample(variants)

    result = apply_noise_mask(sample, (_range(100, 100), _range(110, 110), _range(121, 121)))

    assert result.sample.intervals == {
        "S1_HV1F.json": [[101, 109], [111, 120], [122, 140]],
    }
    assert [variant.pos for variant in result.sample.variants] == [130]
    assert result.sample.variant_flags == flag_variants(result.sample.variants)
    assert [variant.pos for variant in result.excluded_variants] == [100, "110.1", 120]
    assert all(variant.reason == "overlaps likely-noisy range" for variant in result.excluded_variants)
    assert result.excluded_variants[2].noise_range.start == 121


def test_apply_noise_mask_removes_noisy_coverage_without_variants():
    sample = Sample(
        sample_id="S1",
        source_tool=Tool.TRACY,
        intervals={"HV1": [[16024, 16365]], "HV2": [[73, 340]]},
    )

    result = apply_noise_mask(sample, (_range(300, 320),))

    assert result.sample.variants == []
    assert result.sample.intervals == {
        "HV1": [[16024, 16365]],
        "HV2": [[73, 299], [321, 340]],
    }
    assert result.excluded_variants == ()


def test_apply_noise_mask_removes_stale_sample_flags():
    variant = Variant(pos=16189, ref="A", seq="G")
    sample = _sample([variant]).model_copy(
        update={"sample_flags": ["16180-16193 region (16189)"]},
    )

    result = apply_noise_mask(sample, (_range(16189, 16189),))

    assert result.sample.variants == []
    assert "16180-16193 region (16189)" not in result.sample.sample_flags


def test_apply_noise_mask_without_ranges_preserves_sample():
    sample = _sample([Variant(pos=100, ref="A", seq="G")])

    result = apply_noise_mask(sample, ())

    assert result.sample is sample
    assert result.excluded_variants == ()
