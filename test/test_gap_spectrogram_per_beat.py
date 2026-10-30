"""Moments beat alignment and NaN-preserving gap spectra."""

import warnings
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import PipelineDAG
from pipelines import load_pipeline_catalog
from pipelines.low_freq_spectra.metrics import (
    PHASE_METRIC_OUTPUT_NAMES, PHASE_SUMMARY_OUTPUT_NAMES, SPECTRUM_METRIC_OUTPUT_NAMES,
)
from pipelines.low_freq_spectra.per_beat import gap_spectrograms_per_beat
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis
from test_spectrogram_frequency import _write_frequency_table


def test_unequal_beats_align_phase_without_reordering_frequency_or_duplicate_endpoint():
    values = np.column_stack([2 * np.arange(7), 100 + 2 * np.arange(7)])
    result = gap_spectrograms_per_beat(values, np.array([0, 2, 6]), sample_count=4)
    assert result.shape == (4, 2, 2)
    np.testing.assert_array_equal(result[:, 0], [[0, 4], [1, 6], [2, 8], [3, 10]])
    np.testing.assert_array_equal(result[:, 1], result[:, 0] + 100)
    assert result.dtype == np.float32


def test_interpolation_preserves_invalid_gaps_and_exact_valid_samples():
    values = np.array([[1], [np.nan], [3], [4], [np.inf]])
    result = gap_spectrograms_per_beat(values, np.array([0, 4]), sample_count=8)
    np.testing.assert_array_equal(result[:, 0, 0], [1, np.nan, np.nan, np.nan, 3, 3.5, 4, np.nan])
    assert np.isposinf(values[-1, 0])  # Caller input remains unchanged.


@pytest.mark.parametrize("boundaries", [[0], [0, 0], [0, 7], [-1, 4]])
def test_invalid_boundaries_fail_explicitly(boundaries):
    with pytest.raises(ValueError):
        gap_spectrograms_per_beat(np.ones((7, 2)), np.array(boundaries))


def _ctx(output):
    schema = EyeFlowOutputPaths.active()
    if schema.sampling_frequency_hz not in output:
        output.create_dataset(schema.sampling_frequency_hz, data=32000)
    if schema.spectrograms_masked in output:
        _write_frequency_table(output, output[schema.spectrograms_masked].shape[1], 32000)
    return SimpleNamespace(
        output=SimpleNamespace(h5=PipelineH5Output(output)),
    )


def _moments(boundaries, frame_count, source="artery"):
    boundaries = np.asarray(boundaries, dtype=np.int32)
    return SimpleNamespace(
        cycle_boundary_indexes=boundaries,
        cycle_durations_seconds=np.diff(boundaries) * 0.1,
        cardiac_cycle_source=source,
        continuous=lambda vessel, raw: np.zeros(frame_count),
    )


def test_runner_publishes_per_beat_output_with_moments_timing(tmp_path):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        values = np.full((9, 2), 5, dtype=np.float32)
        values[1, 0] = 50
        source = output.create_dataset(schema.spectrograms_masked, data=values)
        source.attrs["unit"] = "power_unit"
        metrics = run_gap_analysis(_ctx(output), _moments([0, 2, 4, 6, 8], 9))
        PipelineH5Output(output).write_many(metrics)
        raw = output[schema.spectrograms_masked_per_beat]
        assert raw.shape == (128, 2, 4)
        assert raw[64, 0, 0] == 50
        assert raw.attrs["unit"] == "power_unit"
        assert raw.attrs["source_path"] == "/ProcessingSpectrograms/Gaps/SpatialMedian"
        assert list(raw.attrs["dimDesc"]) == ["phase", "freq", "beat"]
        np.testing.assert_array_equal(raw.attrs["cycle_boundary_indices"], [0, 2, 4, 6, 8])
        assert raw.attrs["cycle_boundary_workflow"] == "doppler_moments"
        assert raw.attrs["cycle_boundary_source_path"] == (
            "/" + schema.cardiac_cycle.systolic_peak_frame_indices
        )
        svd = output[schema.spectrograms_masked_per_beat_svd_continuous]
        assert svd.name == "/ProcessingSpectrograms/Gaps/PerBeat/SVDContinuous"
        assert svd.shape == (512, 2)
        assert list(svd.attrs["dimDesc"]) == ["time", "freq"]
        assert svd.attrs["time_sampling"] == "uniform_within_each_beat_variable_between_beats"
        assert "phase_samples" not in svd.attrs
        np.testing.assert_array_equal(svd.attrs["within_beat_phase_samples"], raw.attrs["phase_samples"])
        # Two bins cannot support removing two modes: raw concatenation is retained.
        np.testing.assert_array_equal(svd[:128], raw[:, :, 0])
        np.testing.assert_array_equal(svd[128:256], raw[:, :, 1])
        assert svd[64, 0] == 50
        median = output[schema.spectrograms_masked_time_aggregated_median]
        assert median.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat"
        assert "ProcessingSpectrograms/Gaps/TimeAggregated/Median" not in output
        assert median.shape == (2,)
        assert median.dtype == np.float32
        np.testing.assert_array_equal(median[()], [5, 5])
        assert median.attrs["unit"] == "power_unit"
        assert list(median.attrs["dimDesc"]) == ["freq"]
        assert median.attrs["source_path"] == raw.name
        assert "phase_samples" not in median.attrs
        beat_median = output[schema.spectrograms_masked_time_aggregated_median_beat]
        assert beat_median.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianBeat"
        assert beat_median.shape == (128, 2)
        assert beat_median.dtype == np.float32
        np.testing.assert_array_equal(beat_median[:], np.full((128, 2), 5))
        assert list(beat_median.attrs["dimDesc"]) == ["phase", "freq"]
        assert list(beat_median.attrs["aggregation_axes"]) == ["beat"]
        assert beat_median.attrs["source_path"] == raw.name
        assert beat_median.attrs["unit"] == "power_unit"
        np.testing.assert_array_equal(
            beat_median.attrs["phase_samples"], raw.attrs["phase_samples"],
        )
        assert beat_median.attrs["phase_endpoint"] == "excluded"
        for path, count in (
            (schema.spectrograms_masked_time_aggregated_median_early, 43),
            (schema.spectrograms_masked_time_aggregated_median_late, 85),
        ):
            window = output[path]
            np.testing.assert_array_equal(window[:], [5, 5])
            assert window.dtype == np.float32
            assert window.attrs["source_path"] == beat_median.name
            assert window.attrs["source_phase_sample_count"] == count
            assert window.attrs["source_total_phase_sample_count"] == 128
            assert "phase_samples" not in window.attrs
            assert "phase_endpoint" not in window.attrs


def test_time_aggregated_median_pools_samples_ignores_nan_and_preserves_frequency_order(tmp_path):
    schema = EyeFlowOutputPaths.active()
    per_beat = np.array([
        [[0, 10], [1, np.nan], [np.nan, np.nan]],
        [[0, 10], [np.nan, 5], [np.nan, np.nan]],
        [[0, 10], [7, np.nan], [np.nan, np.nan]],
        [[100, 100], [np.nan, 9], [np.nan, np.nan]],
    ], dtype=np.float32)
    original = per_beat.copy()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        source = output.create_dataset(schema.spectrograms_masked, data=np.ones((9, 3)))
        source.attrs["unit"] = "power_unit"
        source.attrs["frequency_bin_centers_hz"] = [30, -10, 20]
        with (
            patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat", return_value=per_beat),
            warnings.catch_warnings(),
        ):
            warnings.simplefilter("error", RuntimeWarning)
            metrics = run_gap_analysis(_ctx(output), _moments([0, 2, 8], 9))
        PipelineH5Output(output).write_many(metrics)
        median = output[schema.spectrograms_masked_time_aggregated_median]
        # Pooled median is 10; median of the two beat medians would be 5.
        np.testing.assert_array_equal(median[()], [10, 6, np.nan])
        np.testing.assert_array_equal(median.attrs["frequency_bin_centers_hz"], [30, -10, 20])
        np.testing.assert_array_equal(median.attrs["valid_sample_count_per_frequency"], [8, 4, 0])
        assert median.attrs["aggregation_method"] == "nanmedian"
        assert list(median.attrs["aggregation_axes"]) == ["phase", "beat"]
        assert median.attrs["source_beat_count"] == 2
        assert median.attrs["source_phase_sample_count"] == 4
        beat_median = output[schema.spectrograms_masked_time_aggregated_median_beat]
        np.testing.assert_array_equal(beat_median[:], [
            [5, 1, np.nan], [5, 5, np.nan], [5, 7, np.nan], [100, 9, np.nan],
        ])
        assert beat_median.shape == (4, 3)
        assert list(beat_median.attrs["dimDesc"]) == ["phase", "freq"]
        assert beat_median.attrs["aggregation_method"] == "nanmedian"
        assert list(beat_median.attrs["aggregation_axes"]) == ["beat"]
        assert beat_median.attrs["aggregation_weighting"] == "equal_weight_per_valid_beat_sample"
        assert beat_median.attrs["source_beat_count"] == 2
        assert beat_median.attrs["source_phase_sample_count"] == 4
        np.testing.assert_array_equal(beat_median.attrs["phase_samples"], [0, .25, .5, .75])
        np.testing.assert_array_equal(beat_median.attrs["frequency_bin_centers_hz"], [30, -10, 20])
        assert "valid_sample_count_per_frequency" not in beat_median.attrs
        early = output[schema.spectrograms_masked_time_aggregated_median_early]
        late = output[schema.spectrograms_masked_time_aggregated_median_late]
        np.testing.assert_array_equal(early[:], [5, 3, np.nan])
        # Taking a pooled raw-sample median over the late window would give 55,
        # whereas the requested median of MedianBeat is 52.5.
        np.testing.assert_array_equal(late[:], [52.5, 8, np.nan])
        for window in (early, late):
            np.testing.assert_array_equal(window.attrs["frequency_bin_centers_hz"], [30, -10, 20])
            assert list(window.attrs["dimDesc"]) == ["freq"]
            assert list(window.attrs["aggregation_axes"]) == ["phase"]
            assert window.attrs["aggregation_method"] == "nanmedian"
            assert window.attrs["aggregation_weighting"] == "equal_weight_per_valid_phase_sample"
            assert window.attrs["source_path"] == beat_median.name
            assert window.attrs["unit"] == "power_unit"
            np.testing.assert_array_equal(window.attrs["valid_sample_count_per_frequency"], [2, 2, 0])
        np.testing.assert_array_equal(per_beat, original)


def test_early_late_spectra_include_one_third_boundary_and_ignore_nan_per_frequency(tmp_path):
    schema = EyeFlowOutputPaths.active()
    beat_spectrum = np.array([
        [0, 1, np.nan], [2, np.nan, np.nan], [10, 5, np.nan],
        [20, 7, np.nan], [30, 9, np.nan], [40, 11, np.nan],
    ], dtype=np.float32)
    raw = np.repeat(beat_spectrum[:, :, None], 2, axis=2)
    original = raw.copy()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        output.create_dataset(schema.spectrograms_masked, data=np.ones((9, 3)))
        with (
            patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat", return_value=raw),
            warnings.catch_warnings(),
        ):
            warnings.simplefilter("error", RuntimeWarning)
            metrics = run_gap_analysis(_ctx(output), _moments([0, 4, 8], 9))
        PipelineH5Output(output).write_many(metrics)
        early = output[schema.spectrograms_masked_time_aggregated_median_early]
        late = output[schema.spectrograms_masked_time_aggregated_median_late]
        assert early.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeatEarly"
        assert late.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeatLate"
        np.testing.assert_array_equal(early[:], [2, 3, np.nan])
        np.testing.assert_array_equal(late[:], [25, 8, np.nan])
        assert early.attrs["source_phase_sample_count"] == 3
        assert late.attrs["source_phase_sample_count"] == 4
        np.testing.assert_array_equal(early.attrs["valid_sample_count_per_frequency"], [3, 2, 0])
        np.testing.assert_array_equal(late.attrs["valid_sample_count_per_frequency"], [4, 4, 0])
        np.testing.assert_allclose(early.attrs["phase_interval"], [0, 1 / 3])
        np.testing.assert_allclose(late.attrs["phase_interval"], [1 / 3, 1])
        assert early.attrs["phase_interval_endpoints"] == late.attrs["phase_interval_endpoints"] == "both_included"
        np.testing.assert_array_equal(raw, original)


@pytest.mark.parametrize("frequency_count", [1, 3])
def test_median_over_a_single_beat_preserves_phase_and_frequency_axes(tmp_path, frequency_count):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        values = np.arange(9)[:, None] + 10 * np.arange(frequency_count)[None, :]
        output.create_dataset(schema.spectrograms_masked, data=values)
        metrics = run_gap_analysis(_ctx(output), _moments([0, 8], 9))
        PipelineH5Output(output).write_many(metrics)
        raw = output[schema.spectrograms_masked_per_beat][:]
        beat_median = output[schema.spectrograms_masked_time_aggregated_median_beat]
        assert beat_median.shape == (128, frequency_count)
        np.testing.assert_array_equal(beat_median[:], raw[:, :, 0])
        for path in (
            schema.spectrograms_masked_time_aggregated_median_early,
            schema.spectrograms_masked_time_aggregated_median_late,
        ):
            assert output[path].shape == (frequency_count,)


@pytest.mark.parametrize("case", [
    "missing_spectrum", "no_detected_cycles", "missing_moments", "two_beats", "no_finite_cells",
])
def test_optional_skips_do_not_discard_existing_spectrum_or_other_products(tmp_path, case):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        output.create_dataset("unrelated", data=42)
        if case != "missing_spectrum":
            values = np.full((9, 2), np.nan) if case == "no_finite_cells" else np.ones((9, 2))
            output.create_dataset(schema.spectrograms_masked, data=values)
        moments = _moments(
            [0, 2, 4, 6, 8] if case == "no_finite_cells" else [0, 4, 8], 9,
            source="none" if case == "no_detected_cycles" else "artery",
        )
        with patch("pipelines.low_freq_spectra.gap_analysis.Logger.log_warning") as warning:
            metrics = run_gap_analysis(
                _ctx(output), None if case == "missing_moments" else moments,
            )
        assert output["unrelated"][()] == 42
        if case in {"two_beats", "no_finite_cells"}:
            # Legacy products survive optional skips; shape outputs are covered
            # by test_gap_shape_metrics, including their all-invalid schema.
            assert set(metrics) >= {
                schema.spectrograms_masked_per_beat,
                schema.spectrograms_masked_per_beat_band_limited,
                schema.spectrograms_masked_per_beat_svd_continuous,
                schema.spectrograms_masked_time_aggregated_median,
                schema.spectrograms_masked_time_aggregated_median_beat,
                schema.spectrograms_masked_time_aggregated_median_early,
                schema.spectrograms_masked_time_aggregated_median_late,
                schema.spectrograms_masked_time_aggregated_std,
                schema.spectrograms_masked_time_aggregated_std_early,
                schema.spectrograms_masked_time_aggregated_std_late,
            } | {
                f"{schema.spectrograms_masked_metrics_root}/{name}"
                for name in SPECTRUM_METRIC_OUTPUT_NAMES.values()
            } | {
                f"{schema.spectrograms_masked_metrics_root}/{key}"
                for key in [
                    "SEI", "SEI_25", "SEI_50", "SEI_75", "SEII_25", "SEII_50", "SEII_75",
                    *PHASE_METRIC_OUTPUT_NAMES.values(),
                    *(f"{name}_{suffix}" for name in PHASE_SUMMARY_OUTPUT_NAMES.values() for suffix in ("early", "late")),
                    *(f"delta_{name}" for name in PHASE_SUMMARY_OUTPUT_NAMES.values()),
                ]
            }
            root = schema.spectrograms_masked_metrics_root
            if case == "no_finite_cells":
                for name in ("index_max1", "index_max2", "A0", "A", "A_50", "SEI_50", "SEI", "delta_A0"):
                    assert np.isnan(metrics[f"{root}/{name}"][0])
                assert any("only NaN" in call.args[0] for call in warning.call_args_list)
            else:
                assert metrics[f"{root}/index_max1"][0] == 0
                assert metrics[f"{root}/index_max2"][0] == 1
                assert metrics[f"{root}/SEI_50"][0] == 1
            assert any("SVD" in call.args[0] for call in warning.call_args_list)
        else:
            assert metrics == {}
            if case != "missing_spectrum":
                warning.assert_called_once()


def test_low_frequency_spectra_are_hidden_and_automatically_follow_velocity():
    available, missing = load_pipeline_catalog()
    assert "gap_analysis" not in {item.name for item in [*available, *missing]}
    descriptors = {item.name: item for item in available}
    assert "gap_per_beat_spectrograms" not in descriptors["velocity"].dag_produces
    spectra = descriptors["low_freq_spectra"]
    assert spectra.visibility == "hidden"
    assert spectra.run_once
    assert spectra.auto_run_after == "velocity"
    assert "gap_per_beat_spectrograms" in spectra.dag_produces
    assert descriptors["velocity"].dag_requires == ()
    dag = PipelineDAG(available)
    assert dag.dependencies_of("velocity") == ()
    assert set(dag.dependencies_of("low_freq_spectra")) == {"velocity", "topology_core"}
    for target in ("velocity_analysis", "spatial_gradient_moment0"):
        names = dag.resolve_targets([target]).names
        assert names.index("velocity") < names.index(target)
        assert names.index("topology_core") < names.index("low_freq_spectra")
        assert names.index("velocity") < names.index("low_freq_spectra")
        assert "gap_analysis" not in names
    assert "low_freq_spectra" not in dag.resolve_targets(["displacement_map"]).names
