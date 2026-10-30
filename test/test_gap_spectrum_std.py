"""Population STD across existing beats after phase-window median reduction."""

import warnings
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis
from pipelines.low_freq_spectra.per_beat import beat_spectrum_std
from calculations.math.phase_windows import early_late_phase_windows


def _raw():
    values = np.full((4, 3, 3), np.nan, dtype=np.float32)
    values[:, 0, :] = [[0, 2, 4], [2, 6, 10], [10, 18, 26], [14, 22, 30]]
    values[:, 1, 0] = [1, np.nan, 3, np.nan]
    return values


def test_std_uses_per_beat_phase_medians_population_weighting_and_nan_counts():
    raw = _raw()
    original = raw.copy()
    windows = {"whole": (np.ones(4, bool), [0, 1]), **early_late_phase_windows(np.arange(4) / 4)}
    for name, expected in (("whole", np.sqrt(24)), ("early", np.sqrt(6)), ("late", np.sqrt(128 / 3))):
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            std, counts = beat_spectrum_std(raw, windows[name][0])
        np.testing.assert_allclose(std[:2], [expected, 0])
        assert np.isnan(std[2])
        np.testing.assert_array_equal(counts, [3, 1, 0])
        assert std.shape == (3,) and std.dtype == np.float32
    np.testing.assert_array_equal(raw, original)


def test_a_single_valid_beat_has_zero_std_and_missing_bins_remain_nan():
    std, counts = beat_spectrum_std(_raw()[:, :, :1], np.ones(4, bool))
    np.testing.assert_array_equal(std, [0, 0, np.nan])
    np.testing.assert_array_equal(counts, [1, 1, 0])


def test_runner_persists_three_raw_power_std_spectra_with_reduction_provenance(tmp_path):
    schema = EyeFlowOutputPaths.active()
    timing = SimpleNamespace(
        cycle_boundary_indexes=np.array([0, 2, 4, 8]),
        cycle_durations_seconds=np.array([0.2, 0.2, 0.4]),
        cardiac_cycle_source="artery", continuous=lambda vessel, raw: np.zeros(9),
    )
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        source = h5.create_dataset(schema.spectrograms_masked, data=np.ones((9, 3)))
        source.attrs["unit"] = "power_unit"
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(h5)))
        with patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat", return_value=_raw()) as align:
            ctx.output.h5.write_many(run_gap_analysis(ctx, timing))
        align.assert_called_once()
        for path, suffix, expected, phases in (
            (schema.spectrograms_masked_time_aggregated_std, "StdBeat", np.sqrt(24), 4),
            (schema.spectrograms_masked_time_aggregated_std_early, "StdBeatEarly", np.sqrt(6), 2),
            (schema.spectrograms_masked_time_aggregated_std_late, "StdBeatLate", np.sqrt(128 / 3), 2),
        ):
            std = h5[path]
            assert std.name == f"/ProcessingSpectrograms/Gaps/TimeAggregated/{suffix}"
            np.testing.assert_allclose(std[:2], [expected, 0])
            assert np.isnan(std[2])
            assert std.attrs["source_path"] == h5[schema.spectrograms_masked_per_beat].name
            assert list(std.attrs["dimDesc"]) == ["freq"]
            assert std.attrs["unit"] == "power_unit"
            assert std.attrs["phase_reduction_method"] == "nanmedian"
            assert std.attrs["aggregation_method"] == "nanstd"
            assert std.attrs["standard_deviation_ddof"] == 0
            assert std.attrs["source_phase_sample_count"] == phases
            assert std.attrs["source_beat_count"] == 3
            np.testing.assert_array_equal(std.attrs["valid_beat_count_per_frequency"], [3, 1, 0])
            assert "phase_samples" not in std.attrs
