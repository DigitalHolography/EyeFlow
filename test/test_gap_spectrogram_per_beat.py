"""Moments beat alignment and NaN-preserving gap spectra."""

from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
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


def test_optional_skips_do_not_discard_existing_spectrum_or_other_products(tmp_path):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        output.create_dataset("unrelated", data=42)
        output.create_dataset(schema.spectrograms_masked, data=np.full((9, 2), np.nan))
        moments = _moments([0, 2, 4, 6, 8], 9, source="artery")
        with patch("pipelines.low_freq_spectra.gap_analysis.Logger.log_warning") as warning:
            metrics = run_gap_analysis(_ctx(output), moments)
        assert output["unrelated"][()] == 42
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
        for name in ("index_max1", "index_max2", "A0", "A", "A_50", "SEI_50", "SEI", "delta_A0"):
            assert np.isnan(metrics[f"{root}/{name}"][0])
        assert any("only NaN" in call.args[0] for call in warning.call_args_list)
        assert any("SVD" in call.args[0] for call in warning.call_args_list)
