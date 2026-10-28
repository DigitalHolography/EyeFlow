"""Scientific references and regular-run contracts for selected shape metrics."""

from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from calculations.math.spectra import describe_spectra, interval_area, inverse_cdf
from input_output.h5_access import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis
from pipelines.low_freq_spectra.metrics import pack_spectrum_metrics
from pipelines.low_freq_spectra.shape_metrics import pack_distribution_metrics, pack_vertical_metrics
from test_spectrogram_frequency import _write_frequency_table


ROOT = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
VERTICAL_NAMES = {
    name for p in (5, 25, 50) for name in (
        f"Av_{p}", f"Av_{p}_t", f"Av_{p}_early", f"Av_{p}_late", f"delta_Av_{p}",
        f"SEv_{p}", f"SEIv_{p}", f"SEIv_{p}_early", f"SEIv_{p}_late", f"delta_SEIv_{p}",
    )
}
DISTRIBUTION_NAMES = {
    "PowerQuantilePercent", "DistributionFrequencyMagnitude", "ValleyMinimum_t",
    "ValleyNormalizedSpectrum_t", "SpectralProbability_t", "WassersteinEarlyLate",
    "WassersteinFrequencyMagnitude",
    *(f"Wasserstein{window}{quantity}" for window in ("Early", "Late")
      for quantity in ("Probability", "CDF", "QuantileFrequency")),
    *(f"{prefix}{name}{suffix}" for name in ("SpectralEntropy", "PowerQuantileFrequency", "F50", "F90")
      for prefix, suffix in (("", "_t"), ("", ""), ("", "_early"), ("", "_late"), ("delta_", ""))),
}


def value(metrics, name):
    return metrics[f"{ROOT}/{name}"][0]


def distributions(phases=None, phase_samples=None, **kwargs):
    phases = np.asarray([[3., 1., 1., 3.], [1., 3., 3., 1.]] if phases is None else phases)
    return pack_distribution_metrics(
        phases, pooled_spectrum=[4., 2., 2., 4.],
        phase_attrs={"unit": "power", "phase_samples": [0., .5] if phase_samples is None else phase_samples},
        frequencies_khz=[1., 2., -2., -1.], peak_indices=[0, 3], **kwargs,
    )


def test_fractional_trapezoids_consult_only_support_and_match_linear_integral():
    # f(x)=2x+1 integrates to x*x+x, even across fractional endpoints.
    signal = 2 * np.arange(8.) + 1
    signal[-1] = np.nan
    assert interval_area(signal, 1.25, 5.5, .4) == pytest.approx(((5.5**2+5.5)-(1.25**2+1.25))*.4)
    assert interval_area(signal, 5, 6, 1) == 12
    assert np.isnan(interval_area(signal, 5, 6.5, 1))


def test_vertical_half_of_each_side_is_not_the_whole_peak_interval():
    pooled = np.ones(64)
    pooled[3], pooled[60] = 10, 9
    phases = np.stack([pooled, pooled * 7])
    attrs = {"unit": "power", "phase_samples": np.array([0., .5])}
    base = pack_spectrum_metrics(
        pooled, source_attrs={"unit": "power"}, sampling_frequency_hz=64000,
        phase_spectrum=phases, phase_attrs=attrs,
        frequency_calibration_attrs={"fft_bin_count": 256, "frequency_binning_factor": 4},
    )
    result = pack_vertical_metrics(pooled, phases, phase_attrs=attrs, base_metrics=base)
    assert {path.rsplit("/", 1)[1] for path in result} == VERTICAL_NAMES
    np.testing.assert_allclose(result[f"{ROOT}/Av_50"][1]["integration_intervals_bin_indices"],
                               [[3, 17.25], [45.75, 60]])
    # Baseline integral 28.5 plus two half-triangle excesses 4.5+4.
    assert value(result, "Av_50") == 37
    assert value(result, "SEIv_50") == pytest.approx(37 / 65.5)
    np.testing.assert_allclose(value(result, "SEv_50"), 37 / 65.5)
    assert value(result, "delta_SEIv_50") == pytest.approx(0)
    assert result[f"{ROOT}/Av_50"][1]["unit"] == "power*kHz"
    assert result[f"{ROOT}/SEIv_50"][1]["unit"] == "1"


def test_entropy_uses_probability_after_valley_normalization_and_is_gain_invariant():
    a = describe_spectra([2., 4., 6., 8.], [1, 2, -2, -1])
    b = describe_spectra(np.array([2., 4., 6., 8.]) * 125, [1, 2, -2, -1])
    np.testing.assert_allclose(a["ValleyNormalizedSpectrum"], [1, 2, 3, 4])
    np.testing.assert_allclose(a["SpectralProbability"], [.1, .2, .3, .4])
    assert a["SpectralEntropy"] == pytest.approx(-np.sum(np.array([.1, .2, .3, .4])*np.log2([.1, .2, .3, .4])))
    for name in a:
        if name != "ValleyMinimum":
            np.testing.assert_allclose(a[name], b[name])
    assert describe_spectra([7.] * 4, [1, 2, -2, -1])["SpectralEntropy"] == 2


def test_quantiles_fold_duplicate_magnitudes_and_keep_discrete_bin_resolution():
    np.testing.assert_array_equal(inverse_cdf([-4, -1, 1, 4], [.1, .2, .3, .4], [0, 25, 50, 75, 100]),
                                  [1, 1, 1, 4, 4])
    np.testing.assert_array_equal(inverse_cdf([0, 1, 2], [0, 1, 0]), 1)


@pytest.mark.parametrize("signal", [[0, 0, 0], [1, np.nan, 2], [1, np.inf, 2], [1, -1, 2]])
def test_invalid_distribution_power_does_not_silently_drop_bins(signal):
    assert all(np.isnan(value).all() for value in describe_spectra(signal, [1, 2, -1]).values())


def test_zero_valley_preserves_defined_entropy_without_infinite_ratios():
    described = describe_spectra([0., 1., 1.], [0, 1, -1])
    assert np.isnan(described["ValleyNormalizedSpectrum"]).all()
    assert described["SpectralEntropy"] == 1
    np.testing.assert_array_equal(described["PowerQuantileFrequency"], 1)


def test_wasserstein_detects_redistribution_at_equal_entropy_and_records_axes():
    result = distributions()
    assert {path.rsplit("/", 1)[1] for path in result} == DISTRIBUTION_NAMES
    assert len(result) == 33
    assert value(result, "delta_SpectralEntropy") == 0
    assert value(result, "F50_early") == 1
    assert value(result, "F50_late") == 2
    assert value(result, "delta_F50") == 1
    assert value(result, "WassersteinEarlyLate") == .5
    x, early, late = (value(result, name) for name in
                      ("WassersteinFrequencyMagnitude", "WassersteinEarlyCDF", "WassersteinLateCDF"))
    assert np.sum(np.abs(early[:-1]-late[:-1])*np.diff(x)) == .5
    assert value(result, "PowerQuantileFrequency_t").shape == (2, 101)
    np.testing.assert_array_equal(value(result, "PowerQuantileFrequency_t")[:, 50], value(result, "F50_t"))
    assert result[f"{ROOT}/SpectralEntropy"][1]["unit"] == "bits"
    assert result[f"{ROOT}/WassersteinEarlyLate"][1]["unit"] == "kHz"


def test_custom_closed_phase_windows_include_shared_boundary_and_allow_empty_windows():
    signal = [[3., 1., 1., 3.], [2., 2., 2., 2.], [1., 3., 3., 1.]]
    result = distributions(signal, np.array([0., .5, .75], dtype=np.float32),
                           early_interval=(0, .5), late_interval=(.5, 1))
    for window in ("early", "late"):
        assert result[f"{ROOT}/SpectralEntropy_{window}"][1]["source_phase_sample_count"] == 2
    empty = distributions(signal, [.5, .6, .7])
    assert np.isnan(value(empty, "SpectralEntropy_early"))
    assert np.isnan(value(empty, "WassersteinEarlyLate"))


@pytest.mark.parametrize("invalid", [False, True])
def test_regular_gap_runner_writes_exact_selected_schema_and_calibrated_provenance(tmp_path, invalid):
    schema = EyeFlowOutputPaths.active()
    spectrum = np.ones(64)
    spectrum[3], spectrum[60] = 10, 9
    if invalid:
        spectrum[:] = np.nan
    with h5py.File(tmp_path / "output.h5", "w") as output:
        _write_frequency_table(output, 64, 64000)
        output.create_dataset(schema.spectrograms_masked, data=np.broadcast_to(spectrum, (9, 64)))
        output[schema.spectrograms_masked].attrs["unit"] = "power"
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(output)))
        moments = SimpleNamespace(
            cycle_boundary_indexes=np.array([0, 2, 8]), cycle_durations_seconds=np.array([.2, .6]),
            cardiac_cycle_source="artery", continuous=lambda vessel, raw: np.zeros(9),
        )
        result = run_gap_analysis(ctx, moments, early_interval=(0, .5), late_interval=(.5, 1))
        ctx.output.h5.write_many(result)
        group = output[ROOT]
        selected = VERTICAL_NAMES | DISTRIBUTION_NAMES
        assert len(selected) == 63 and len(group) == 127
        assert selected <= set(group)
        assert "SEIv_10" not in group
        for name in selected:
            ds = group[name]
            assert ds.dtype == np.float32
            assert len(ds.attrs["dimDesc"]) == ds.ndim
            if "source_path" in ds.attrs:
                assert ds.attrs["source_path"] in output
            assert ds.attrs["fft_bin_count"] == 256
            assert ds.attrs["frequency_binning_factor"] == 4
            assert ds.attrs["sampling_frequency_hz"] == 64000
        assert group["ValleyNormalizedSpectrum_t"].shape == (128, 64)
        assert group["PowerQuantileFrequency_t"].shape == (128, 101)
        for name in ("SEIv_50_early", "SpectralEntropy_early"):
            np.testing.assert_array_equal(group[name].attrs["phase_interval"], [0, .5])
        if invalid:
            for name in ("SEIv_50", "SpectralEntropy", "F50", "WassersteinEarlyLate"):
                assert np.isnan(group[name][()])
        else:
            assert group["SEIv_50"][()] == pytest.approx(37 / 65.5)
            assert group["WassersteinEarlyLate"][()] == 0
