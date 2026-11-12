"""Scientific references and regular-run contracts for selected shape metrics."""

from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from calculations.math.spectra import describe_spectra, interval_area
from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis
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


def test_fractional_trapezoids_consult_only_support_and_match_linear_integral():
    # f(x)=2x+1 integrates to x*x+x, even across fractional endpoints.
    signal = 2 * np.arange(8.) + 1
    signal[-1] = np.nan
    assert interval_area(signal, 1.25, 5.5, .4) == pytest.approx(((5.5**2+5.5)-(1.25**2+1.25))*.4)
    assert interval_area(signal, 5, 6, 1) == 12
    assert np.isnan(interval_area(signal, 5, 6.5, 1))


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


def test_regular_gap_runner_writes_exact_selected_schema_and_calibrated_provenance(tmp_path):
    invalid = False
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
