"""Masked band-power reduction, spatial coordinates, and gap skip behavior."""

from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.preprocessing import (
    preprocess_gap_spectrograms,
    spatial_median_band_power,
)
from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms


class FrameSource:
    """A lazy source that rejects attempts to read the complete video."""

    def __init__(self, values):
        self.values = values
        self.shape = values.shape
        self.ndim = values.ndim
        self.dtype = values.dtype
        self.reads = []

    def __getitem__(self, index):
        assert isinstance(index, int), "Read only one time frame at a time."
        self.reads.append(index)
        return self.values[index]


def test_nonsquare_xy_mask_uses_inverted_hd_y_and_preserves_time_bin_order():
    # In persisted EF coordinates these samples are (x=0, y=1) and (x=2, y=0).
    first = np.array([[2, 100, 101], [102, 103, 10]], dtype=np.float32)
    values = np.stack([np.stack([first, first * 3]), np.stack([first * 5, first * 7])])
    mask = np.zeros((3, 2), dtype=bool)
    mask[0, 1] = mask[2, 0] = True
    source = FrameSource(values)

    result = spatial_median_band_power(source, mask)

    np.testing.assert_array_equal(result, [[6, 18], [30, 42]])
    assert source.reads == [0, 1]


def test_nonfinite_included_samples_become_nan_excluded_samples_do_not_contribute():
    values = np.array([[
        [[-2, np.nan, 100], [np.inf, 100, 10]],
        [[np.inf, 100, 100], [100, 100, 10]],
        [[2, 100, 100], [100, 100, np.nan]],
        [[2, 100, 100], [100, 100, -np.inf]],
    ]], dtype=np.float32)
    mask = np.zeros((3, 2), dtype=bool)
    mask[0, 1] = mask[2, 0] = True

    np.testing.assert_array_equal(
        spatial_median_band_power(values, mask), [[4, np.nan, np.nan, np.nan]],
    )


def test_negative_medians_are_zeroed_after_reduction_while_nonfinite_medians_stay_nan():
    values = np.array([[
        [[-10, 2]], [[-2, 10]], [[-3, -1]], [[np.nan, 1]], [[-np.inf, 1]],
    ]], dtype=np.float32)
    original = values.copy()
    zones = np.ones((1, 2), dtype=bool)
    metrics = preprocess_gap_spectrograms(values, zones)
    power, attrs = metrics[EyeFlowOutputPaths.active().spectrograms_masked]
    # Clipping individual pixels before the even-sample median would give 1 and 5.
    np.testing.assert_array_equal(power, [[0, 4, 0, np.nan, np.nan]])
    assert attrs["negative_value_policy"] == "clamp_spatial_median_to_zero"
    np.testing.assert_array_equal(values, original)


@pytest.mark.parametrize("preprocess", [preprocess_gap_spectrograms, preprocess_spectrograms])
def test_masked_median_rejects_spatial_outliers_and_preserves_invalid_pairs(preprocess):
    values = np.array([
        [[[1, 3, 1000, np.nan]], [[2, 8, 2000, np.inf]]],
        [[[-100, 2, 4, 99999]], [[-10, -2, -1, 99999]]],
        [[[1, 3, np.nan, 99999]], [[2, np.inf, 4, 99999]]],
    ], dtype=np.float32)
    original = values.copy()
    zones = np.array([[True, True, True, False]])
    source = FrameSource(values)

    power, attrs = preprocess(source, zones)[EyeFlowOutputPaths.active().spectrograms_masked]

    np.testing.assert_array_equal(power, [[3, 8], [2, 0], [np.nan, np.nan]])
    assert power.dtype == np.float64
    assert source.reads == [0, 1, 2]
    assert attrs["quantity"] == "spatial_median_band_power"
    assert attrs["spatial_aggregation_method"] == "median"
    assert attrs["negative_value_policy"] == "clamp_spatial_median_to_zero"
    assert attrs["gap_pixel_count"] == 3
    np.testing.assert_array_equal(values, original)


@pytest.mark.parametrize("preprocess", [preprocess_gap_spectrograms, preprocess_spectrograms])
def test_even_sample_median_uses_float64_to_preserve_large_finite_power(preprocess):
    values = np.full((1, 1, 1, 2), 3e38, dtype=np.float32)
    with np.errstate(over="raise", invalid="raise"):
        power, _ = preprocess(values, np.ones((1, 2), bool))[
            EyeFlowOutputPaths.active().spectrograms_masked
        ]
    np.testing.assert_array_equal(power, [[float(values[0, 0, 0, 0])]])


@pytest.mark.parametrize("factor", [4, 2])
def test_nearest_neighbor_interpolation_and_hdf5_schema(tmp_path, factor):
    native_zones = np.zeros((4, 6), dtype=bool)
    native_zones[:2, :2] = True
    native_zones[2:, 4:] = True
    native_zones[:2, 2:4] = True
    native_zones[1, 3] = False  # The sampled center is excluded; other block pixels are irrelevant.
    if factor == 4:
        native_zones = native_zones.repeat(2, axis=0).repeat(2, axis=1)
    values = np.array([[[[2, 100, 100], [100, 100, 10]]]], dtype=np.float32)
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "source.h5", "w") as source, h5py.File(tmp_path / "out.h5", "w") as out:
        spectra = source.create_dataset("spectrograms", data=values)
        spectra.attrs["unit"] = "power_unit"
        PipelineH5Output(out).write_many(preprocess_gap_spectrograms(spectra, native_zones))
        dataset = out[schema.spectrograms_masked]
        assert dataset.name == "/ProcessingSpectrograms/Gaps/SpatialMedian"
        assert "ProcessingSpectrograms/Gaps/Raw" not in out
        assert dataset.shape == (1, 1)
        assert dataset.dtype == np.float32
        assert dataset[0, 0] == 6
        assert list(dataset.attrs["dimDesc"]) == ["time", "freq"]
        assert dataset.attrs["unit"] == "power_unit"
        assert dataset.attrs["quantity"] == "spatial_median_band_power"
        assert dataset.attrs["spatial_aggregation_method"] == "median"
        assert dataset.attrs["gap_pixel_count"] == 2
        assert dataset.attrs["mask_reduction"] == "nearest_neighbor_interpolation"
        assert dataset.attrs["mask_interpolation_order"] == 0
        assert dataset.attrs["mask_interpolation_grid"] == "pixel_centers"
        assert "/spectrograms" in dataset.attrs["mask_target_shape_source"]
        assert list(dataset.attrs["mask_application_shape"]) == [3, 2]
        assert dataset.attrs["spectrogram_spatial_transform"] == "flip_D2_then_swap_D2_D3"


@pytest.mark.parametrize("reduced_only", [False, True])
def test_empty_native_or_reduced_support_skips_analysis_without_reading_source(reduced_only):
    zones = np.zeros((4, 6), dtype=bool)
    if reduced_only:
        zones[0, 0] = True  # This pixel is not selected by center-based interpolation.
    source = FrameSource(np.ones((2, 3, 2, 3), dtype=np.float32))
    with patch("pipelines.low_freq_spectra.preprocessing.Logger.log_warning") as warning:
        assert preprocess_gap_spectrograms(source, zones) == {}
    assert source.reads == []
    warning.assert_called_once()
    assert "Skipping gap analysis" in warning.call_args.args[0]
    assert "empty" in warning.call_args.args[0]


def test_missing_spectra_omit_output_with_warning():
    with patch("pipelines.low_freq_spectra.preprocessing.Logger.log_warning") as warning:
        assert preprocess_gap_spectrograms(None, np.ones((2, 3), bool)) == {}
    warning.assert_called_once()
    assert "absent" in warning.call_args.args[0]


def test_target_shape_comes_from_hdf5_spectrograms_for_noninteger_nonsquare_resize(tmp_path):
    native_zones = np.zeros((5, 7), dtype=bool)
    native_zones[0, 6] = native_zones[2, 2] = native_zones[4, 0] = True
    original = native_zones.copy()
    image = np.arange(12, dtype=np.float32).reshape(3, 4)
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "HD.h5", "w") as source, h5py.File(tmp_path / "EF.h5", "w") as output:
        spectra = source.create_dataset("spectrograms", data=np.stack([image, image * 3])[None])
        metrics = preprocess_gap_spectrograms(spectra, native_zones)
        PipelineH5Output(output).write_many(metrics)
        result = output[schema.spectrograms_masked]
        # Nearest-neighbor centers select native spectral pixels (0,3), (1,1), (2,0).
        np.testing.assert_array_equal(result[()], [[5, 15]])
        np.testing.assert_array_equal(result.attrs["mask_native_shape"], [5, 7])
        np.testing.assert_array_equal(result.attrs["mask_application_shape"], [4, 3])
        assert result.attrs["gap_pixel_count"] == 3
        assert result.attrs["mask_reduction"] == "nearest_neighbor_interpolation"
    np.testing.assert_array_equal(native_zones, original)
