"""Peak selection, missing-value contracts, and scalar gap-spectrum HDF5 outputs."""

import warnings
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.metrics import (
    PHASE_METRIC_NAMES, SPECTRUM_METRIC_DEFINITIONS, analyze_spectrum, pack_spectrum_metrics,
)
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis


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
        "A0": 216,
        "A1": 45,
        "A2": 27,
        "R1": 45 / 216,
        "R2": 27 / 216,
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
    assert result["R1"] == 1


def test_interval_nan_propagates_area_while_peak_selection_ignores_it():
    spectrum = _spectrum()
    spectrum[8] = np.nan
    result = analyze_spectrum(spectrum)
    assert (result["index_max1"], result["index_max2"]) == (2, 28)
    assert result["A0"] == 216
    assert np.isnan(result["A1"])
    assert np.isnan(result["R1"])
    assert result["A2"] == 27
    assert result["R2"] == 27 / 216
    spectrum[8], spectrum[0] = 1, np.nan
    assert analyze_spectrum(spectrum)["A1"] == 45  # NaNs outside the interval do not affect it.


def test_missing_valley_preserves_other_samples_for_the_valley_median():
    spectrum = _spectrum()
    spectrum[15] = np.nan
    result = analyze_spectrum(spectrum)
    assert np.isnan(result["min1"])
    assert result["A0"] == 216
    assert np.isnan(result["A1"])
    assert result["A2"] == 27
    assert result["R2"] == 27 / 216


@pytest.mark.parametrize("length, expected_height", [(2, 3), (4, 1.5), (6, 1), (32, 1)])
def test_valley_neighborhood_length_is_centered_and_uses_same_rectangle_width(length, expected_height):
    result = analyze_spectrum(_spectrum(), valley_neighborhood_length=length)
    assert result["A2"] == 27 * expected_height
    assert result["R2"] == expected_height / 8


def test_valley_height_uses_robust_median_rather_than_mean_or_minimum():
    spectrum = _spectrum()
    spectrum[12:20] = [1, 2, 3, 4, 5, 6, 7, 1000]
    result = analyze_spectrum(spectrum)
    assert result["A2"] == (19 - 2 + 1) * 4.5


def test_all_missing_valley_returns_nan_without_warning_and_preserves_baseline():
    spectrum = _spectrum()
    spectrum[12:20] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = analyze_spectrum(spectrum)
    assert result["A0"] == 216
    for name in ("A2", "R2"):
        assert np.isnan(result[name])


@pytest.mark.parametrize("length", [0, 1, 3, 33, 34, -2, 6.0, True])
def test_invalid_valley_neighborhood_length_is_rejected(length):
    with pytest.raises(ValueError, match="even integer from 2 through 32"):
        analyze_spectrum(_spectrum(), valley_neighborhood_length=length)


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
    for name in ("A0", "A1", "A2", "R1", "R2"):
        assert np.isnan(result[name])
    assert result["index_min1"] == 15
    assert result["index_min2"] == 16


def test_zero_baseline_returns_nan_ratio_without_division_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = analyze_spectrum(np.zeros(32))
    assert np.isnan(result["R1"])
    assert np.isnan(result["R2"])
    assert result["A0"] == result["A1"] == result["A2"] == 0


@pytest.mark.parametrize("shape", [(), (0,), (1,), (2, 16), (3, 1)])
def test_unsupported_shape_returns_all_thirteen_nan_metrics_with_one_warning(shape):
    with patch("pipelines.low_freq_spectra.metrics.Logger.log_warning") as warning:
        result = analyze_spectrum(np.ones(shape))
    assert set(result) == set(SPECTRUM_METRIC_DEFINITIONS)
    assert len(result) == 13
    assert all(np.isnan(value) for value in result.values())
    warning.assert_called_once()
    assert "at least 2 bins" in warning.call_args.args[0]


@pytest.mark.parametrize(
    "bin_count, start, stop, height", [(16, 6, 10, 1.5), (32, 12, 20, 3.5), (64, 24, 40, 7.5)],
)
def test_default_neighborhood_is_quarter_of_actual_domain_and_metadata_matches(bin_count, start, stop, height):
    spectrum = np.ones(bin_count)
    spectrum[start:stop] = np.arange(stop - start)
    spectrum[0], spectrum[-1] = 100, 80
    result = analyze_spectrum(spectrum)
    assert result["index_max1"] == 0
    assert result["index_max2"] == bin_count - 1
    assert result["index_min1"] == bin_count // 2 - 1
    assert result["index_min2"] == bin_count // 2
    assert result["A0"] == bin_count * 90
    assert result["A2"] == bin_count * height
    assert result["R2"] == height / 90
    metrics = pack_spectrum_metrics(spectrum, source_attrs={"unit": "power_unit"})
    for _, attrs in metrics.values():
        assert attrs["valley_neighborhood_length"] == bin_count // 4
        np.testing.assert_array_equal(attrs["valley_neighborhood_indices"], np.arange(start, stop))
        assert attrs["valley_neighborhood_selection"] == "quarter_of_bin_domain"


@pytest.mark.parametrize(
    "bin_count, start, stop", [(2, 0, 2), (3, 0, 2), (12, 4, 8), (31, 11, 19), (33, 12, 20), (63, 23, 39), (65, 24, 40)],
)
def test_finite_non_multiple_of_eight_domains_produce_valid_peak_indices(bin_count, start, stop):
    spectrum = np.arange(bin_count, dtype=float)
    with patch("pipelines.low_freq_spectra.metrics.Logger.log_warning") as warning:
        result = analyze_spectrum(spectrum)
    warning.assert_not_called()
    assert result["index_max1"] == bin_count // 2 - 1
    assert result["index_max2"] == bin_count - 1
    assert result["A2"] == (bin_count - bin_count // 2 + 1) * (start + stop - 1) / 2
    metrics = pack_spectrum_metrics(spectrum, source_attrs={})
    _, attrs = next(iter(metrics.values()))
    np.testing.assert_array_equal(attrs["valley_neighborhood_indices"], np.arange(start, stop))
    assert attrs["source_bin_count"] == bin_count
    assert attrs["peak_half_valid_bin_counts"] == [bin_count // 2, bin_count - bin_count // 2]


def test_missing_first_half_retains_nan_index_and_explains_missing_input():
    spectrum = np.arange(63, dtype=float)
    spectrum[:31] = np.nan
    with patch("pipelines.low_freq_spectra.metrics.Logger.log_warning") as warning:
        metrics = pack_spectrum_metrics(spectrum, source_attrs={})
    path = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    index, attrs = metrics[path + "/index_max1"]
    assert np.isnan(index)
    assert metrics[path + "/index_max2"][0] == 62
    assert attrs["peak_half_valid_bin_counts"] == [0, 32]
    warning.assert_called_once()
    assert "index_max1" in warning.call_args.args[0]
    assert "only NaN" in warning.call_args.args[0]


@pytest.mark.parametrize("bin_count, valley_length, expected_length", [
    (32, None, 8), (32, 2, 2), (32, 6, 6), (64, None, 16),
    (63, None, 16), (65, None, 16), (12, None, 4), (3, None, 2), (3, 2, 2),
])
def test_runner_persists_scalar_metrics_from_the_raw_aggregate_spectrum(
    tmp_path, bin_count, valley_length, expected_length,
):
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
        with patch("pipelines.low_freq_spectra.metrics.Logger.log_warning") as warning:
            metrics = run_gap_analysis(ctx, moments, valley_neighborhood_length=valley_length)
        ctx.output.h5.write_many(metrics)
        median = output[schema.spectrograms_masked_time_aggregated_median]
        assert median.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat"
        np.testing.assert_array_equal(median[()], spectrum)
        group = output[schema.spectrograms_masked_metrics_root]
        expected_keys = set(SPECTRUM_METRIC_DEFINITIONS) | {
            key for name in PHASE_METRIC_NAMES
            for key in (f"{name}_t", f"{name}_early", f"{name}_late", f"delta_{name}")
        }
        assert group.name == "/ProcessingSpectrograms/Gaps/Metrics"
        assert set(group) == expected_keys
        assert not {"h0", "h1", "hv", "rho_h", "rho_v"} & set(group)
        assert "ProcessingSpectrograms/Gaps/TimeAggregated/Metrics" not in output
        for name, dataset in group.items():
            assert dataset.shape == ((128,) if name.endswith("_t") else ())
            assert dataset.dtype == np.float32
            assert not dataset.name.endswith("/value")
            assert list(dataset.attrs["dimDesc"]) == (["phase"] if name.endswith("_t") else [])
            metric_name = name.removeprefix("delta_").split("_")[0]
            assert dataset.attrs["unit"] == (
                "1" if name.startswith("index_") or metric_name in {"R1", "R2"} else "power_unit"
            )
            assert dataset.attrs["valley_neighborhood_length"] == expected_length
            middle = bin_count // 2
            np.testing.assert_array_equal(
                dataset.attrs["valley_neighborhood_indices"],
                np.arange(middle - expected_length // 2, middle + expected_length // 2),
            )
            assert dataset.attrs["valley_aggregation_method"] == "nanmedian"
            if name.startswith("index_"):
                assert dataset.attrs["index_base"] == 0
            expected = analyze_spectrum(spectrum, valley_neighborhood_length=valley_length)
            if name in expected:
                assert dataset.attrs["source_path"] == median.name
                assert dataset[()] == np.float32(expected[name])
            elif name.startswith("delta_"):
                assert dataset[()] == 0  # A constant source has identical early/late medians.
            else:
                np.testing.assert_allclose(dataset[()], expected[metric_name])
        assert "ProcessingBandRatio/Spectrograms" not in output
        warning.assert_not_called()
