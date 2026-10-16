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
from pipelines.low_freq_spectra.metrics import (
    PHASE_METRIC_NAMES, analyze_phase_spectrum, pack_spectrum_metrics,
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


def _pack(values=None, phase=None):
    if values is None:
        values = _phase_spectrum()
    if phase is None:
        phase = np.arange(len(values), dtype=np.float32) / len(values)
    return pack_spectrum_metrics(
        _phase_spectrum()[0], source_attrs={"unit": "power_unit"},
        phase_spectrum=values,
        phase_attrs={"unit": "power_unit", "phase_samples": phase, "cycle_boundary_workflow": "doppler_moments"},
    )


def test_phase_areas_use_fixed_peak_bins_and_leave_input_unchanged():
    values = _phase_spectrum()
    original = values.copy()
    curves = analyze_phase_spectrum(values, index_max1=1, index_max2=6)
    expected = {
        "A0": [54, 108, 30, 18, 12, 42],
        "A1": [26, 52, 16, 18, 20, 34],
        "A2": [18, 36, 12, 24, 30, 36],
    }
    for name, result in expected.items():
        np.testing.assert_array_equal(curves[name], result)
    np.testing.assert_allclose(curves["R1"], np.array(expected["A1"]) / expected["A0"])
    np.testing.assert_allclose(curves["R2"], np.array(expected["A2"]) / expected["A0"])
    np.testing.assert_array_equal(values, original)
    # Huge maxima outside the fixed interval at phase 1 do not move its bounds.
    assert curves["A0"][1] == 108
    override = analyze_phase_spectrum(values, index_max1=1, index_max2=6, valley_neighborhood_length=4)
    assert override["A2"][0] == 9


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_early_late_medians_include_shared_third_boundary_and_deltas_have_signed_difference(dtype):
    phase = np.arange(6, dtype=dtype) / 6
    metrics = _pack(phase=phase)
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    expected = {
        "A0": (54, 24, 30), "A1": (26, 19, 7), "A2": (18, 27, -9),
        "R1": (13 / 27, 19 / 21, -80 / 189), "R2": (1 / 3, 23 / 21, -16 / 21),
    }
    for name, (early, late, delta) in expected.items():
        value, attrs = metrics[f"{root}/{name}_t"]
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
            result, aggregate_attrs = metrics[f"{root}/{name}_{suffix}"]
            assert result == pytest.approx(target)
            assert aggregate_attrs["source_phase_sample_count"] == (3 if suffix == "early" else 4)
            assert aggregate_attrs["source_path"] == f"/{root}/{name}_t"
            assert aggregate_attrs["dimDesc"] == []
            assert "phase_samples" not in aggregate_attrs
        assert metrics[f"{root}/delta_{name}"][0] == pytest.approx(delta)
        assert metrics[f"{root}/delta_{name}"][1]["source_paths"] == [
            f"/{root}/{name}_early", f"/{root}/{name}_late",
        ]
    # Aggregate the ratios themselves, rather than dividing aggregated areas.
    assert metrics[f"{root}/R1_late"][0] != pytest.approx(19 / 24)


def test_default_128_phase_grid_uses_43_early_and_85_late_samples():
    metrics = _pack(np.tile(_phase_spectrum()[0], (128, 1)))
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    for name in PHASE_METRIC_NAMES:
        assert metrics[f"{root}/{name}_early"][1]["source_phase_sample_count"] == 43
        assert metrics[f"{root}/{name}_late"][1]["source_phase_sample_count"] == 85


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
        curves = analyze_phase_spectrum(values, index_max1=1, index_max2=6)
    assert np.isnan(curves["A0"][0]) and curves["A2"][0] == 18
    assert curves["A0"][1] == 108 and np.isnan(curves["A1"][1])
    assert curves["A2"][2] == 18  # One remaining valley sample.
    assert np.isnan(curves["A2"][3]) and curves["A0"][3] == 18
    assert curves["A0"][4] == 0 and curves["A2"][4] == 30
    assert np.isnan(curves["R1"][4]) and np.isnan(curves["R2"][4])
    assert curves["A1"][5] == 34


def test_nan_windows_ignore_missing_phases_and_leave_missing_deltas_nan():
    values = _phase_spectrum()
    values[:3] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        metrics = _pack(values)
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    for name in PHASE_METRIC_NAMES:
        assert np.isnan(metrics[f"{root}/{name}_early"][0])
        assert np.isfinite(metrics[f"{root}/{name}_late"][0])
        assert np.isnan(metrics[f"{root}/delta_{name}"][0])
        assert metrics[f"{root}/{name}_early"][1]["valid_sample_count"] == 0
    assert metrics[f"{root}/A0_late"][0] == 18


def test_missing_fixed_peak_indices_do_not_trigger_a_phasewise_peak_search():
    curves = analyze_phase_spectrum(_phase_spectrum(), index_max1=np.nan, index_max2=6)
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
        h5.create_dataset(schema.spectrograms_masked, data=np.ones((9, 8)))
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(h5)))
        with patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat", return_value=raw):
            ctx.output.h5.write_many(run_gap_analysis(ctx, timing))
        group = h5[schema.spectrograms_masked_metrics_root]
        assert group["index_max1"][()] == 1
        assert group["index_max2"][()] == 6
        # Median of per-beat A0 at the first phase would instead be 303.
        np.testing.assert_array_equal(group["A0_t"][:], [60] * 6)
        np.testing.assert_allclose(group["R1_t"][:], [0.4] * 6)
        assert group["A0_early"][()] == 60
        assert group["A0_late"][()] == 60
        assert group["delta_A0"][()] == 0
        assert group["A0_t"].attrs["source_path"] == h5[schema.spectrograms_masked_time_aggregated_median_beat].name
        assert not {"h0", "h1", "hv", "rho_h", "rho_v"} & set(group)
