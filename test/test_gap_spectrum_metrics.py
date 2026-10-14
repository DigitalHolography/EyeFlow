"""Peak selection, missing-value contracts, and scalar gap-spectrum HDF5 outputs."""

import warnings
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.velocity.gaps.metrics import SPECTRUM_METRIC_DEFINITIONS, analyze_spectrum
from pipelines.velocity.gaps.runner import run_gap_analysis


def _spectrum():
    values = np.ones(32, dtype=np.float32)
    values[2], values[28], values[15], values[16] = 10, 6, 2, 4
    return values


def test_metrics_match_hand_calculated_inclusive_areas_and_leave_input_unchanged():
    spectrum = _spectrum()
    original = spectrum.copy()
    result = analyze_spectrum(spectrum)
    assert result == {
        "index_max1": 2,
        "max1": 10,
        "index_max2": 28,
        "max2": 6,
        "index_min1": 15,
        "min1": 2,
        "index_min2": 16,
        "min2": 4,
        "h0": 8,
        "h1": 3,
        "rho_h": 0.375,
        "A0": 216,
        "A1": 45,
        "rho_A": 45 / 216,
    }
    np.testing.assert_array_equal(spectrum, original)


def test_peak_selection_uses_one_half_each_lowest_tie_index_and_includes_endpoints():
    spectrum = _spectrum()
    spectrum[0] = spectrum[2] = 100
    spectrum[31] = 6
    result = analyze_spectrum(spectrum)
    # Both globally highest samples are in the first half; still use one per half.
    assert result["index_max1"] == 0
    assert result["index_max2"] == 28
    spectrum[31] = 7
    assert analyze_spectrum(spectrum)["index_max2"] == 31
    spectrum[15], spectrum[16] = 200, 300
    result = analyze_spectrum(spectrum)
    assert (result["index_max1"], result["index_max2"]) == (15, 16)
    assert result["A0"] == result["A1"] == 500
    assert result["rho_A"] == 1


def test_interval_nan_propagates_area_while_peak_selection_ignores_it():
    spectrum = _spectrum()
    spectrum[8] = np.nan
    result = analyze_spectrum(spectrum)
    assert (result["index_max1"], result["index_max2"]) == (2, 28)
    assert result["h0"] == 8
    assert result["rho_h"] == 0.375
    assert result["A0"] == 216
    assert np.isnan(result["A1"])
    assert np.isnan(result["rho_A"])
    spectrum[8], spectrum[0] = 1, np.nan
    assert analyze_spectrum(spectrum)["A1"] == 45  # NaNs outside the interval do not affect it.


def test_missing_valley_does_not_use_a_one_sided_mean():
    spectrum = _spectrum()
    spectrum[15] = np.nan
    result = analyze_spectrum(spectrum)
    assert np.isnan(result["min1"])
    assert np.isnan(result["h1"])
    assert np.isnan(result["rho_h"])
    assert result["h0"] == 8
    assert result["A0"] == 216
    assert np.isnan(result["A1"])


@pytest.mark.parametrize("half", [slice(0, 16), slice(16, 32), slice(0, 32)])
def test_empty_peak_half_propagates_missing_peak_without_warnings(half):
    spectrum = _spectrum()
    spectrum[half] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = analyze_spectrum(spectrum)
    missing_peak = "1" if half.start == 0 else "2"
    assert np.isnan(result[f"index_max{missing_peak}"])
    assert np.isnan(result[f"max{missing_peak}"])
    for name in ("h0", "h1", "rho_h", "A0", "A1", "rho_A"):
        assert np.isnan(result[name])
    assert result["index_min1"] == 15
    assert result["index_min2"] == 16


def test_zero_baseline_returns_nan_ratio_without_division_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = analyze_spectrum(np.zeros(32))
    assert np.isnan(result["rho_h"])
    assert np.isnan(result["rho_A"])
    assert result["h0"] == result["h1"] == result["A0"] == result["A1"] == 0


@pytest.mark.parametrize("shape", [(0,), (2,), (31,), (33,), (2, 16)])
def test_unsupported_shape_returns_all_fourteen_nan_metrics_with_one_warning(shape):
    with patch("pipelines.velocity.gaps.metrics.Logger.log_warning") as warning:
        result = analyze_spectrum(np.ones(shape))
    assert set(result) == set(SPECTRUM_METRIC_DEFINITIONS)
    assert len(result) == 14
    assert all(np.isnan(value) for value in result.values())
    warning.assert_called_once()
    assert "32 bins" in warning.call_args.args[0]


@pytest.mark.parametrize("bin_count", [32, 3])
def test_runner_persists_scalar_metrics_from_the_raw_aggregate_spectrum(tmp_path, bin_count):
    schema = EyeFlowOutputPaths.active()
    spectrum = _spectrum() if bin_count == 32 else np.ones(bin_count)
    with h5py.File(tmp_path / "output.h5", "w") as output:
        source = output.create_dataset(
            schema.spectrograms_masked,
            data=np.broadcast_to(spectrum, (9, bin_count)),
        )
        source.attrs["unit"] = "power_unit"
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(output)))
        moments = SimpleNamespace(
            cycle_boundary_indexes=np.array([0, 2, 8]),
            cycle_durations_seconds=np.array([0.2, 0.6]),
            cardiac_cycle_source="artery",
            continuous=lambda vessel, raw: np.zeros(9),
        )
        with patch("pipelines.velocity.gaps.metrics.Logger.log_warning") as warning:
            metrics = run_gap_analysis(ctx, moments)
        ctx.output.h5.write_many(metrics)
        median = output[schema.spectrograms_masked_time_aggregated_median]
        np.testing.assert_array_equal(median[()], spectrum)
        group = output[schema.spectrograms_masked_time_aggregated_metrics_root]
        assert set(group) == set(SPECTRUM_METRIC_DEFINITIONS)
        for name, dataset in group.items():
            assert dataset.shape == ()
            assert dataset.dtype == np.float32
            assert not dataset.name.endswith("/value")
            assert dataset.attrs["source_path"] == median.name
            assert len(dataset.attrs["dimDesc"]) == 0
            assert dataset.attrs["unit"] == (
                "1" if name.startswith("index_") or name in {"rho_h", "rho_A"} else "power_unit"
            )
            if name.startswith("index_"):
                assert dataset.attrs["index_base"] == 0
            if bin_count == 32:
                assert dataset[()] == np.float32(analyze_spectrum(spectrum)[name])
            else:
                assert np.isnan(dataset[()])
        assert "ProcessingBandRatio/Spectrograms" not in output
        if bin_count == 32:
            warning.assert_not_called()
        else:
            assert any("32 bins" in call.args[0] for call in warning.call_args_list)
