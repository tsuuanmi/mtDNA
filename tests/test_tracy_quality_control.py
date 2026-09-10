from typing import Any

import pytest
from pydantic import ValidationError

from src.config import TracySettings
from src.tools.tracy.quality_control import NoiseConfig, analyze_trace_data, build_qc_report


def _config() -> NoiseConfig:
    return NoiseConfig(
        mask_enabled=True,
        window_size=15,
        window_step=5,
        min_valid_bases=10,
        min_supporting_windows=2,
        snr_threshold=3.0,
        low_snr_threshold=2.0,
        purity_threshold=0.60,
        background_threshold=0.35,
        second_peak_threshold=0.50,
        quality_threshold=20.0,
        signal_threshold=75.0,
        bad_base_fraction=0.40,
        suspicious_base_fraction=0.20,
    )


def _trace(length=30, noisy_indices=frozenset(), ref_align=None, ref_start=100) -> dict[str, Any]:
    ref_align = ref_align or "A" * length
    alt_align = "A" * len(ref_align)
    peak_a = []
    peak_c = []
    peak_g = []
    peak_t = []
    qualities = []
    for index in range(len(alt_align)):
        if index in noisy_indices:
            peaks = (100, 80, 70, 60)
            quality = 10
        else:
            peaks = (1000, 10, 10, 10)
            quality = 40
        a, c, g, t = peaks
        peak_a.append(a)
        peak_c.append(c)
        peak_g.append(g)
        peak_t.append(t)
        qualities.append(quality)
    return {
        "ref1pos": ref_start,
        "ref1align": ref_align,
        "alt1align": alt_align,
        "basecallPos": list(range(len(alt_align))),
        "basecallQual": qualities,
        "peakA": peak_a,
        "peakC": peak_c,
        "peakG": peak_g,
        "peakT": peak_t,
    }


def test_noise_window_settings_are_consistent():
    assert TracySettings().noise_mask_enabled is True
    with pytest.raises(ValidationError, match="noise_min_valid_bases"):
        TracySettings(noise_window_size=10, noise_min_valid_bases=11)
    with pytest.raises(ValidationError, match="greater than or equal to 2"):
        TracySettings(noise_min_supporting_windows=1)


def test_clean_trace_has_no_noisy_ranges():
    result = analyze_trace_data(
        _trace(),
        sample_id="S1",
        filename="S1_HV1F.json",
        config=_config(),
    )

    assert result.status == "clean"
    assert result.bases_evaluated == 30
    assert result.ranges == ()
    assert all(window.classification == "clean" for window in result.windows)


def test_single_noisy_window_is_suspicious_without_mask_range():
    result = analyze_trace_data(
        _trace(length=10, noisy_indices=frozenset(range(10))),
        sample_id="S1",
        filename="S1_HV1F.json",
        config=_config(),
    )

    assert result.status == "suspicious"
    assert result.ranges == ()
    assert result.windows[0].classification == "likely_noisy"


def test_overlapping_noisy_windows_are_merged():
    result = analyze_trace_data(
        _trace(noisy_indices=frozenset(range(5, 20))),
        sample_id="S1",
        filename="S1_HV1F.json",
        config=_config(),
    )

    assert result.status == "likely_noisy"
    assert len(result.ranges) == 1
    noise_range = result.ranges[0]
    assert noise_range.start == 100
    assert noise_range.end == 124
    assert noise_range.alignment_start == 0
    assert noise_range.alignment_end == 24
    assert noise_range.supporting_windows == 3

    report = build_qc_report([result], _config())
    assert report["likely_noisy"] == 1
    assert report["analyzed"] == 1
    assert report["traces"][0]["ranges"][0]["start"] == 100


def test_isolated_noisy_base_does_not_create_a_range():
    result = analyze_trace_data(
        _trace(noisy_indices=frozenset({7})),
        sample_id="S1",
        filename="S1_HV1F.json",
        config=_config(),
    )

    assert result.status == "clean"
    assert result.ranges == ()


def test_missing_trace_channel_is_unavailable():
    data = _trace()
    data.pop("peakT")

    result = analyze_trace_data(
        data,
        sample_id="S1",
        filename="S1_HV1F.json",
        config=_config(),
    )

    assert result.status == "unavailable"
    assert result.reason == "Missing Tracy fields: peakT"


def test_reverse_trace_uses_canonical_reference_range():
    result = analyze_trace_data(
        _trace(length=20, noisy_indices=frozenset(range(20)), ref_start=16000),
        sample_id="S1",
        filename="S1_HV1R.json",
        config=_config(),
    )

    assert result.status == "likely_noisy"
    assert result.ranges[0].start == 16000
    assert result.ranges[0].end == 16019
    assert result.strand == "reverse"


def test_reference_gap_does_not_shift_later_coordinates():
    ref_align = "AAAA-AAAAAAAAAAAAAAA"
    result = analyze_trace_data(
        _trace(
            length=len(ref_align),
            noisy_indices=frozenset(range(len(ref_align))),
            ref_align=ref_align,
            ref_start=100,
        ),
        sample_id="S1",
        filename="S1_HV1F.json",
        config=_config(),
    )

    assert result.status == "likely_noisy"
    assert result.ranges[0].start == 100
    assert result.ranges[0].end == 118
