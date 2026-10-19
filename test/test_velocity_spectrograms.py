"""Full-frame spectra, optional gap reduction, and new shared HDF5 hierarchy."""

from unittest.mock import patch

import h5py
import numpy as np

from input_output.h5_access import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms
from test_gap_spectrogram_preprocessing import FrameSource


def test_full_frame_mean_median_and_aligned_gap_median_read_each_frame_once(tmp_path):
    first = np.array([[2, 100, 100], [100, 100, 10]], dtype=np.float32)
    values = np.stack([first[None], (first * 2)[None]])
    zones = np.zeros((4, 6), bool)
    zones[:2, :2] = zones[2:, 4:] = True
    zones[:2, 2:4] = True
    source = FrameSource(values)
    schema = EyeFlowOutputPaths.active()
    metrics = preprocess_spectrograms(source, zones)
    assert source.reads == [0, 1]
    with h5py.File(tmp_path / "out.h5", "w") as output:
        PipelineH5Output(output).write_many(metrics)
        mean = output[schema.spectrograms_spatial_average]
        median = output[schema.spectrograms_spatial_median]
        gap = output[schema.spectrograms_masked]
        np.testing.assert_allclose(mean[:], [[412 / 6], [824 / 6]])
        np.testing.assert_array_equal(median[:], [[100], [200]])
        # Included values are 2, 100, 10; the median differs from the gap mean.
        np.testing.assert_array_equal(gap[:], [[10], [20]])
        assert gap.attrs["quantity"] == "spatial_median_band_power"
        assert gap.attrs["spatial_aggregation_method"] == "median"
        assert mean.name == "/ProcessingSpectrograms/Global/SpatialAverage"
        assert median.name == "/ProcessingSpectrograms/Global/SpatialMedian"
        assert "ProcessingSpectrograms/SpatialAverage" not in output
        assert "ProcessingSpectrograms/SpatialMedian" not in output
        assert gap.name == "/ProcessingSpectrograms/Gaps/SpatialMedian"
        assert "ProcessingSpectrograms/Gaps/Raw" not in output
        assert "spectrogramsMasked" not in output
        assert mean.dtype == median.dtype == gap.dtype == np.float32
        assert list(median.attrs["dimDesc"]) == ["time", "freq"]


def test_global_invalid_samples_do_not_poison_valid_gap_support():
    values = np.array([[[[2, np.nan, 100], [np.inf, 100, 10]]]], dtype=np.float32)
    zones = np.zeros((2, 3), bool)
    zones[0, 0] = zones[1, 2] = True
    schema = EyeFlowOutputPaths.active()
    metrics = preprocess_spectrograms(values, zones)
    assert np.isnan(metrics[schema.spectrograms_spatial_average][0]).all()
    assert np.isnan(metrics[schema.spectrograms_spatial_median][0]).all()
    np.testing.assert_array_equal(metrics[schema.spectrograms_masked][0], [[6]])


def test_empty_or_missing_gap_support_keeps_full_frame_outputs():
    schema = EyeFlowOutputPaths.active()
    for zones in (None, np.zeros((4, 6), bool)):
        source = FrameSource(np.ones((2, 3, 2, 3), dtype=np.float32))
        with patch("pipelines.low_freq_spectra.preprocessing.Logger.log_warning"):
            metrics = preprocess_spectrograms(source, zones)
        assert set(metrics) == {schema.spectrograms_spatial_average, schema.spectrograms_spatial_median}
        assert source.reads == [0, 1]


def test_finite_negative_reductions_are_clamped_after_aggregation():
    values = np.array([[[[-10, -2, 7]], [[-10, 2, 11]]]], dtype=np.float32)
    schema = EyeFlowOutputPaths.active()
    metrics = preprocess_spectrograms(values, None)
    np.testing.assert_array_equal(metrics[schema.spectrograms_spatial_average][0], [[0, 1]])
    np.testing.assert_array_equal(metrics[schema.spectrograms_spatial_median][0], [[0, 2]])


def test_missing_spectra_skip_all_reductions_with_warning():
    with patch("pipelines.low_freq_spectra.spectrograms.Logger.log_warning") as warning:
        assert preprocess_spectrograms(None, None) == {}
    warning.assert_called_once()
