"""Fixed spectral geometry, phase-window medians, and shared gap-metric outputs."""

import warnings
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis
from test_spectrogram_frequency import _write_frequency_table
from pipelines.low_freq_spectra.metrics import (
    PHASE_METRIC_OUTPUT_NAMES, PHASE_SUMMARY_OUTPUT_NAMES,
    analyze_phase_spectrum, pack_spectrum_metrics,
)


def _phase_spectrum():
    return np.array([
        [0, 10, 1, 2, 4, 1, 8, 0],
        [1000, 20, 2, 4, 8, 2, 16, 1000],
        [0, 6, 1, 1, 3, 1, 4, 0],
        [0, 4, 2, 3, 5, 2, 2, 0],
        [0, 2, 3, 4, 6, 3, 2, 0],
        [0, 8, 4, 5, 7, 4, 6, 0],
    ], dtype=np.float32)


def _pack(values=None, phase=None, **window_options):
    if values is None:
        values = _phase_spectrum()
    if phase is None:
        phase = np.arange(len(values), dtype=np.float32) / len(values)
    return pack_spectrum_metrics(
        _phase_spectrum()[0], sampling_frequency_hz=8000, source_attrs={"unit": "power_unit"},
        phase_spectrum=values,
        phase_attrs={"unit": "power_unit", "phase_samples": phase, "cycle_boundary_workflow": "doppler_moments"},
        **window_options,
    )


def test_phase_areas_use_fixed_peak_bins_and_leave_input_unchanged():
    values = _phase_spectrum()
    original = values.copy()
    curves = analyze_phase_spectrum(values, sampling_frequency_hz=8000, index_max1=1, index_max2=6)
    expected = {
        "A1": [17, 34, 11, 15, 18, 27],
        "A0": [15, 30, 10, 20, 25, 30],
    }
    for name, result in expected.items():
        np.testing.assert_array_equal(curves[name], result)
    np.testing.assert_allclose(curves["SE"], (np.array(expected["A1"]) - expected["A0"]) / expected["A0"])
    np.testing.assert_array_equal(values, original)
    # Huge maxima outside the fixed interval at phase 1 do not move its bounds.
    assert curves["A1"][1] == 34
    override = analyze_phase_spectrum(values, sampling_frequency_hz=8000, index_max1=1, index_max2=6, valley_neighborhood_length=4)
    assert override["A0"][0] == 7.5


def test_early_late_medians_include_shared_third_boundary_and_deltas_have_signed_difference():
    phase = np.arange(6, dtype=np.float32) / 6
    metrics = _pack(phase=phase)
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    expected = {
        "A1": (17, 16.5, -0.5), "A0": (15, 22.5, 7.5),
        "SE": (2 / 15, -7 / 40, -37 / 120),
    }
    for name, (early, late, delta) in expected.items():
        curve_name = PHASE_METRIC_OUTPUT_NAMES[name]
        summary_name = PHASE_SUMMARY_OUTPUT_NAMES[name]
        value, attrs = metrics[f"{root}/{curve_name}"]
        assert value.shape == (6,)
        assert value.dtype == np.float32
        assert attrs["source_path"].endswith("/TimeAggregated/MedianBeat")
        assert attrs["peak_source_path"].endswith("/TimeAggregated/MedianPhaseBeat")
        assert attrs["fixed_peak_indices"] == [1, 6]
        assert attrs["dimDesc"] == ["phase"]
        assert attrs["cycle_boundary_workflow"] == "doppler_moments"
        assert "aggregation_method" not in attrs
        np.testing.assert_array_equal(attrs["phase_samples"], phase)
        for suffix, target in (("early", early), ("late", late)):
            result, aggregate_attrs = metrics[f"{root}/{summary_name}_{suffix}"]
            assert result == pytest.approx(target)
            assert aggregate_attrs["source_phase_sample_count"] == (3 if suffix == "early" else 4)
            assert aggregate_attrs["source_path"] == f"/{root}/{curve_name}"
            assert aggregate_attrs["dimDesc"] == []
            assert "phase_samples" not in aggregate_attrs
        assert metrics[f"{root}/delta_{summary_name}"][0] == pytest.approx(delta)
        delta_attrs = metrics[f"{root}/delta_{summary_name}"][1]
        assert delta_attrs["source_paths"] == [
            f"/{root}/{summary_name}_late", f"/{root}/{summary_name}_early",
        ]
        assert delta_attrs["definition"] == f"{summary_name}_late minus {summary_name}_early"
        assert delta_attrs["calculation_method"] == "late_median_minus_early_median"
    # Aggregate the ratios themselves, rather than dividing aggregated areas.
    assert metrics[f"{root}/SEI_late"][0] != pytest.approx(16.5 / 22.5 - 1)


def test_runner_shares_custom_intervals_across_spectra_std_and_metric_summaries(tmp_path):
    schema = EyeFlowOutputPaths.active()
    values = _phase_spectrum()
    values[:, 6] *= 2  # Keep the pooled second-half peak at bin 6 rather than bin 4.
    raw = np.stack((values, 2 * values), axis=2)
    timing = SimpleNamespace(
        cycle_boundary_indexes=np.array([0, 4, 8]),
        cycle_durations_seconds=np.array([0.4, 0.4]),
        cardiac_cycle_source="artery", continuous=lambda vessel, raw: np.zeros(9),
    )
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.spectrograms_masked, data=np.ones((9, 8)))
        h5.create_dataset(schema.sampling_frequency_hz, data=8000)
        _write_frequency_table(h5, 8, 8000)
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(h5)))
        with patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat", return_value=raw):
            metrics = run_gap_analysis(ctx, timing, early_interval=(0, 0.25), late_interval=(0.5, 1))
        for path, interval, count in (
            (schema.spectrograms_masked_time_aggregated_median_early, [0, 0.25], 2),
            (schema.spectrograms_masked_time_aggregated_median_late, [0.5, 1], 3),
            (schema.spectrograms_masked_time_aggregated_std_early, [0, 0.25], 2),
            (schema.spectrograms_masked_time_aggregated_std_late, [0.5, 1], 3),
        ):
            assert metrics[path][1]["phase_interval"] == interval
            assert metrics[path][1]["source_phase_sample_count"] == count
        assert metrics[schema.spectrograms_masked_time_aggregated_median_early][0][3] == 4.5
        assert metrics[schema.spectrograms_masked_time_aggregated_median_late][0][3] == 6
        assert metrics[schema.spectrograms_masked_time_aggregated_std_early][0][3] == 1.5
        assert metrics[schema.spectrograms_masked_time_aggregated_std_late][0][3] == 2
        assert metrics[f"{schema.spectrograms_masked_metrics_root}/index_max1"][0] == 1
        assert metrics[f"{schema.spectrograms_masked_metrics_root}/index_max2"][0] == 6
        assert metrics[f"{schema.spectrograms_masked_metrics_root}/delta_A0"][0] == 3.75


def test_phase_nan_and_zero_policies_preserve_independent_area_information():
    values = _phase_spectrum()
    values[0, 1] = np.nan
    values[1, 2] = np.nan
    values[2, 3] = np.nan
    values[3, 3:5] = np.nan
    values[4, [1, 6]] = 0
    values[5, [0, 7]] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        curves = analyze_phase_spectrum(values, sampling_frequency_hz=8000, index_max1=1, index_max2=6)
    assert np.isnan(curves["A_50"][0]) and curves["A0"][0] == 15
    assert np.isnan(curves["A_50"][1]) and np.isnan(curves["A1"][1])
    assert curves["A0"][2] == 15  # One remaining valley sample.
    assert np.isnan(curves["A0"][3]) and np.isnan(curves["A_50"][3])
    # Zero heights at the fixed bounds do not erase the actual interior peaks.
    assert curves["A_50"][4] == pytest.approx(125 / 12) and curves["A0"][4] == 25
    assert curves["SE"][4] == -9 / 25 and curves["SE_50"][4] == pytest.approx(67 / 125)
    for name in ("SE", "SE_25", "SE_50", "SE_75"):
        assert np.isnan(curves[name][3])  # The new A0 denominator is missing.
    assert curves["A1"][5] == 27
