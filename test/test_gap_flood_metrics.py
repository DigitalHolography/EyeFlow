"""Fractional peak-height flooding and spectral excess contracts."""

import warnings

import numpy as np
import pytest

from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.metrics import analyze_phase_spectrum, analyze_spectrum, pack_spectrum_metrics


@pytest.mark.parametrize("gain, sampling_frequency_hz", [(1, 4000), (0.1, 4000), (3, 8000)])
def test_flooded_v_areas_include_between_bin_crossings_and_scale_in_power_khz(gain, sampling_frequency_hz):
    # Peak bins 0 and 2, unit kHz spacing at fs=4000. Each ramp is linear.
    # At level h, clipped V area = 2*h - h*h/8, rather than just h.
    values = np.array([8., 0., 8., 0.]) * gain
    original = values.copy()
    result = analyze_spectrum(values, sampling_frequency_hz=sampling_frequency_hz)
    scale = gain * sampling_frequency_hz / 4000
    assert result["A1"] == 8 * scale
    for percent, area in ((25, 3.5), (50, 6), (75, 7.5)):
        assert result[f"A_{percent}"] == pytest.approx(area * scale)
    np.testing.assert_array_equal(values, original)


def test_flood_height_uses_arithmetic_mean_of_asymmetric_peaks():
    # Mean=6, so h_50=3. The right ramp remains below h, with area=1.
    result = analyze_spectrum([10., 0., 2., 0.], sampling_frequency_hz=4000)
    assert result["A_50"] == pytest.approx(3 - 3**2 / 20 + 1)


def test_phase_flooding_uses_actual_peak_heights_with_fixed_bounds_and_preserves_se():
    values = np.array([[8., 0., 8., 0.], [16., 0., 16., 1000.], [8., 4., 8., 0.]])
    curves = analyze_phase_spectrum(values, index_max1=0, index_max2=2, sampling_frequency_hz=4000)
    # The huge peak outside the interval changes the flood height only.
    np.testing.assert_allclose(curves["A_50"], [6, 16, 8])
    np.testing.assert_allclose(curves["SE_50"], [1/3, 0, 0.5])
    np.testing.assert_allclose(curves["A1"], [8, 16, 12])
    np.testing.assert_allclose(curves["SE"], [0, 0, 0])
    assert set(curves) == {"A0", "A1", "SE", "A_25", "A_50", "A_75", "SE_25", "SE_50", "SE_75",
                           "SEII_25", "SEII_50", "SEII_75"}


def test_pooled_flood_scalar_is_independent_of_phase_summaries_and_has_provenance():
    spectrum = np.array([8., 0., 8., 0.])
    phase_spectrum = np.array([2, 3, 4])[:, None] * spectrum
    products = pack_spectrum_metrics(
        spectrum, sampling_frequency_hz=4000, source_attrs={"unit": "power"},
        phase_spectrum=phase_spectrum,
        phase_attrs={"unit": "power", "phase_samples": np.array([0., 1/3, 2/3])},
    )
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root

    def get(name):
        return products[f"{root}/{name}"]

    assert get("A_50")[0] == 6
    np.testing.assert_allclose(get("A_50_t")[0], [12, 18, 24])
    assert get("A_50_early")[0] == 15
    assert get("A_50_late")[0] == 21
    assert get("delta_A_50")[0] == 6
    for name in ("SEI_50", "SEI_50_early", "SEI_50_late"):
        assert get(name)[0] == pytest.approx(1/3)
        assert get(name)[1]["source_path"] == f"/{root}/SE_50"
        assert get(name)[1]["unit"] == "1"
    assert get("delta_SEI_50")[0] == 0
    assert get("A_50")[1]["flood_height_evaluation"] == "pooled_spectrum"
    assert get("A_50_t")[1]["flood_height_evaluation"] == "per_phase"
    assert get("A_50_t")[1]["flood_height_fraction"] == 0.5
    assert get("A_50_t")[1]["flood_height_reference"] == "arithmetic_mean_of_phasewise_maxima_per_frequency_half"
    assert get("A_50_t")[1]["flood_integration_bounds"] == "fixed_pooled_peak_indices"
    assert get("A_50_t")[1]["unit"] == "power*kHz"
    assert get("SE_50")[1]["reference_area_source_path"] == f"/{root}/A_50_t"
    assert get("SE_50")[1]["denominator_validity_rule"] == "finite_and_strictly_positive_A_50"
    assert get("SEI_50")[1]["valid_sample_count"] == 3


def test_flood_references_require_finite_nonnegative_interval_and_ratios_require_positive_area():
    values = np.array([
        [8., np.nan, 8., 0.], [8., np.inf, 8., 0.], [np.nan, 0., 8., 0.],
        [8., -1., 8., 0.], [0., 0., 0., 0.], [8., 0., 8., np.nan],
    ])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        curves = analyze_phase_spectrum(values, index_max1=0, index_max2=2, sampling_frequency_hz=4000)
    for percent in (25, 50, 75):
        assert np.isnan(curves[f"A_{percent}"][:4]).all()
        assert curves[f"A_{percent}"][4] == 0
        assert np.isnan(curves[f"SE_{percent}"][:5]).all()
        assert np.isfinite(curves[f"SE_{percent}"][5])


def test_moving_peaks_do_not_move_integration_bounds_or_change_original_areas():
    # Peaks move from bins 1/6 to 2/5 with the same maxima. Integrating the
    # second spectrum on moving bounds would give 14 instead of 24.
    values = np.array([
        [0, 8, 2, 3, 3, 2, 8, 0],
        [0, 2, 8, 3, 3, 8, 2, 0],
    ], dtype=float)
    curves = analyze_phase_spectrum(values, index_max1=1, index_max2=6, sampling_frequency_hz=8000)
    # Exact clipped ramp areas at h=4, retaining all five bin intervals.
    np.testing.assert_allclose(curves["A_50"], [46/3, 272/15])
    np.testing.assert_allclose(curves["A1"], [18, 24])
    np.testing.assert_allclose(curves["A0"], [15, 15])
    np.testing.assert_allclose(curves["SE"], [0.2, 0.6])


@pytest.mark.parametrize("gain, sampling_frequency_hz", [(1, 4000), (3, 8000)])
def test_seii_uses_flooded_area_over_valley_reference_and_summarizes_ratios(gain, sampling_frequency_hz):
    spectra = gain * np.array([[8., 0., 8., 0.], [8., 4., 8., 0.], [8., 6., 8., 0.]])
    products = pack_spectrum_metrics(
        spectra[0], sampling_frequency_hz=sampling_frequency_hz, source_attrs={"unit": "power"},
        phase_spectrum=spectra, phase_attrs={"unit": "power", "phase_samples": np.array([0., 1/3, 2/3])},
    )
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    expected = {25: [-9/16, -2/3, -5/7], 50: [-1/4, -1/3, -3/7], 75: [-1/16, -1/12, -1/7]}
    for percent, curve in expected.items():
        name = f"SEII_{percent}"
        actual, attrs = products[f"{root}/{name}_t"]
        np.testing.assert_allclose(actual, curve)
        assert attrs["reference_area_source_path"] == f"/{root}/A0_t"
        assert attrs["signal_area_source_path"] == f"/{root}/A_{percent}_t"
        assert attrs["denominator_validity_rule"] == "finite_and_strictly_positive_A0"
        assert attrs["unit"] == "1" and actual.dtype == np.float32
        early, late = np.mean(curve[:2]), np.mean(curve[1:])
        for suffix, value in (("", np.median(curve)), ("_early", early), ("_late", late)):
            actual, attrs = products[f"{root}/{name}{suffix}"]
            assert actual == pytest.approx(value)
            assert attrs["source_path"] == f"/{root}/{name}_t"
        assert products[f"{root}/delta_{name}"][0] == pytest.approx(late - early)


def test_seii_requires_positive_valley_reference_even_with_finite_flooded_area():
    spectra = np.array([[0., 8., 1., 0., 0., 1., 8., 0.]])
    curves = analyze_phase_spectrum(spectra, index_max1=1, index_max2=6, sampling_frequency_hz=8000)
    assert curves["A0"][0] == 0
    for percent in (25, 50, 75):
        assert curves[f"A_{percent}"][0] > 0
        assert np.isnan(curves[f"SEII_{percent}"][0])
