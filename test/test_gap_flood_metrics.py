"""Fractional peak-height flooding and spectral excess contracts."""

import warnings

import numpy as np
import pytest

from pipelines.low_freq_spectra.metrics import analyze_phase_spectrum, analyze_spectrum


def test_flooded_v_areas_include_between_bin_crossings_and_scale_in_power_khz():
    gain, sampling_frequency_hz = 3, 8000
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
