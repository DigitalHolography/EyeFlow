"""Fixed spectral geometry, phase-window medians, and shared gap-metric outputs."""

import warnings
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import PipelineH5Output
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
        "A2": [45, 90, 25, 15, 10, 35],
        "A1": [17, 34, 11, 15, 18, 27],
        "A0": [15, 30, 10, 20, 25, 30],
    }
    for name, result in expected.items():
        np.testing.assert_array_equal(curves[name], result)
    np.testing.assert_allclose(curves["R1"], np.array(expected["A1"]) / expected["A0"])
    np.testing.assert_allclose(curves["R2"], np.array(expected["A2"]) / expected["A0"])
    np.testing.assert_allclose(curves["SE"], (np.array(expected["A1"]) - expected["A0"]) / expected["A0"])
    np.testing.assert_allclose(curves["Z2"], (np.array(expected["A2"]) - expected["A0"]) / expected["A0"])
    np.testing.assert_array_equal(values, original)
    # Huge maxima outside the fixed interval at phase 1 do not move its bounds.
    assert curves["A2"][1] == 90
    override = analyze_phase_spectrum(values, sampling_frequency_hz=8000, index_max1=1, index_max2=6, valley_neighborhood_length=4)
    assert override["A0"][0] == 7.5


def test_phase_areas_integrate_in_khz_and_keep_ratio_units_dimensionless():
    curves = analyze_phase_spectrum(
        _phase_spectrum()[:1], index_max1=1, index_max2=6, sampling_frequency_hz=4000,
    )
    assert curves["A0"][0] == 7.5
    assert curves["A1"][0] == 8.5
    assert curves["A2"][0] == 22.5
    assert curves["R1"][0] == 17 / 15
    assert curves["R2"][0] == 3


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_early_late_medians_include_shared_third_boundary_and_deltas_have_signed_difference(dtype):
    phase = np.arange(6, dtype=dtype) / 6
    metrics = _pack(phase=phase)
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    expected = {
        "A2": (45, 20, -25), "A1": (17, 16.5, -0.5), "A0": (15, 22.5, 7.5),
        "R1": (17 / 15, 33 / 40, -37 / 120),
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
    assert metrics[f"{root}/R1_late"][0] != pytest.approx(16.5 / 22.5)


def test_default_128_phase_grid_uses_43_early_and_85_late_samples():
    metrics = _pack(np.tile(_phase_spectrum()[0], (128, 1)))
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    for name in PHASE_SUMMARY_OUTPUT_NAMES.values():
        assert metrics[f"{root}/{name}_early"][1]["source_phase_sample_count"] == 43
        assert metrics[f"{root}/{name}_late"][1]["source_phase_sample_count"] == 85


def test_custom_intervals_apply_to_every_metric_and_delta_metadata():
    metrics = _pack(early_interval=(0, 0.25), late_interval=(0.5, 1))
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    assert metrics[f"{root}/A0_early"][0] == 22.5
    assert metrics[f"{root}/A0_late"][0] == 25
    assert metrics[f"{root}/delta_A0"][0] == 2.5
    assert metrics[f"{root}/A_early"][0] == 25.5
    assert metrics[f"{root}/A_late"][0] == 18
    assert metrics[f"{root}/delta_A"][0] == -7.5
    for name in PHASE_SUMMARY_OUTPUT_NAMES.values():
        early, early_attrs = metrics[f"{root}/{name}_early"]
        late, late_attrs = metrics[f"{root}/{name}_late"]
        delta, delta_attrs = metrics[f"{root}/delta_{name}"]
        assert early_attrs["phase_interval"] == delta_attrs["early_phase_interval"] == [0, 0.25]
        assert late_attrs["phase_interval"] == delta_attrs["late_phase_interval"] == [0.5, 1]
        assert early_attrs["source_phase_sample_count"] == 2
        assert late_attrs["source_phase_sample_count"] == 3
        assert delta == pytest.approx(float(late) - float(early), abs=1e-6)


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
    assert np.isnan(curves["A2"][0]) and curves["A0"][0] == 15
    assert curves["A2"][1] == 90 and np.isnan(curves["A1"][1])
    assert curves["A0"][2] == 15  # One remaining valley sample.
    assert np.isnan(curves["A0"][3]) and curves["A2"][3] == 15
    assert curves["A2"][4] == 0 and curves["A0"][4] == 25
    assert curves["R1"][4] == 16 / 25 and curves["R2"][4] == 0
    for name in ("R1", "R2", "SE", "Z2"):
        assert np.isnan(curves[name][3])  # The new A0 denominator is missing.
    assert curves["A1"][5] == 27


def test_nan_windows_ignore_missing_phases_and_leave_missing_deltas_nan():
    values = _phase_spectrum()
    values[:3] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        metrics = _pack(values)
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    for name in PHASE_SUMMARY_OUTPUT_NAMES.values():
        assert np.isnan(metrics[f"{root}/{name}_early"][0])
        assert np.isfinite(metrics[f"{root}/{name}_late"][0])
        assert np.isnan(metrics[f"{root}/delta_{name}"][0])
        assert metrics[f"{root}/{name}_early"][1]["valid_sample_count"] == 0
    assert metrics[f"{root}/A2_late"][0] == 15


def test_zero_valley_denominator_keeps_areas_but_makes_ratios_and_offsets_nan():
    values = _phase_spectrum()[:1].copy()
    values[:, 3:5] = 0
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        curves = analyze_phase_spectrum(values, sampling_frequency_hz=8000, index_max1=1, index_max2=6)
    assert curves["A0"][0] == 0
    assert curves["A1"][0] == 11
    assert curves["A2"][0] == 45
    assert all(np.isnan(curves[name][0]) for name in ("R1", "R2", "SE", "Z2"))


def test_missing_fixed_peak_indices_do_not_trigger_a_phasewise_peak_search():
    curves = analyze_phase_spectrum(_phase_spectrum(), sampling_frequency_hz=8000, index_max1=np.nan, index_max2=6)
    assert all(np.isnan(curve).all() for curve in curves.values())


def test_pooled_and_phase_frequency_domains_must_agree():
    with pytest.raises(ValueError, match="same frequency-bin domain"):
        _pack(np.ones((6, 16)))


def test_runner_calculates_curves_after_the_beat_median(tmp_path):
    schema = EyeFlowOutputPaths.active()
    raw = np.ones((6, 8, 3), dtype=np.float32)
    raw[:, [0, 7], :] = 0
    raw[:, [1, 6], :] = 10
    raw[0, 1, :] = [100, 1, 10]
    raw[0, 6, :] = [1, 100, 10]
    timing = SimpleNamespace(
        cycle_boundary_indexes=np.array([0, 2, 4, 8]),
        cycle_durations_seconds=np.array([0.2, 0.2, 0.4]),
        cardiac_cycle_source="artery", continuous=lambda vessel, raw: np.zeros(9),
    )
    with h5py.File(tmp_path / "output.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=8000)
        _write_frequency_table(h5, 8, 8000)
        h5.create_dataset(schema.spectrograms_masked, data=np.ones((9, 8)))
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(h5)))
        with patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat", return_value=raw):
            ctx.output.h5.write_many(run_gap_analysis(ctx, timing))
        group = h5[schema.spectrograms_masked_metrics_root]
        assert group["index_max1"][()] == 1
        assert group["index_max2"][()] == 6
        # Median of per-beat A2 at the first phase would instead be 252.5.
        np.testing.assert_array_equal(group["A2_t"][:], [50] * 6)
        np.testing.assert_array_equal(group["A0_t"][:], [5] * 6)
        np.testing.assert_allclose(group["R1_t"][:], [2.8] * 6)
        np.testing.assert_allclose(group["SE"][:], [1.8] * 6)
        assert group["A2_early"][()] == 50
        assert group["A2_late"][()] == 50
        assert group["delta_A2"][()] == 0
        assert group["A2_t"].attrs["source_path"] == h5[schema.spectrograms_masked_time_aggregated_median_beat].name
        assert not {"A1", "A1_t", "R2", "R2_t", "Z2", "Z2_t"} & set(group)
        assert not {"h0", "h1", "hv", "rho_h", "rho_v"} & set(group)
