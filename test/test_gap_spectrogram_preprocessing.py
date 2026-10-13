"""Masked band-power reduction, spatial coordinates, and gap skip behavior."""

from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.velocity.gaps.preprocessing import (
    preprocess_gap_spectrograms,
    spatial_mean_band_power,
)


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

    result = spatial_mean_band_power(source, mask)

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

    np.testing.assert_array_equal(spatial_mean_band_power(values, mask), [[4, np.nan, np.nan, np.nan]])


def test_negative_means_are_zeroed_after_averaging_while_nonfinite_means_stay_nan():
    values = np.array([[
        [[-10, 2]], [[-2, 10]], [[-3, -1]], [[np.nan, 1]], [[-np.inf, 1]],
    ]], dtype=np.float32)
    original = values.copy()
    zones = np.ones((1, 2), dtype=bool)
    metrics = preprocess_gap_spectrograms(values, zones)
    power, attrs = metrics[EyeFlowOutputPaths.active().spectrograms_masked]
    # Clipping individual pixels before averaging would incorrectly give 1 and 5.
    np.testing.assert_array_equal(power, [[0, 4, 0, np.nan, np.nan]])
    assert attrs["negative_value_policy"] == "clamp_spatial_mean_to_zero"
    np.testing.assert_array_equal(values, original)


def test_conservative_reduction_and_hdf5_schema(tmp_path):
    native_zones = np.zeros((4, 6), dtype=bool)
    native_zones[:2, :2] = True
    native_zones[2:, 4:] = True
    native_zones[:2, 2:4] = True
    native_zones[0, 2] = False  # Three true pixels do not qualify this block.
    values = np.array([[[[2, 100, 100], [100, 100, 10]]]], dtype=np.float32)
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "source.h5", "w") as source, h5py.File(tmp_path / "out.h5", "w") as out:
        spectra = source.create_dataset("spectrograms", data=values)
        spectra.attrs["unit"] = "power_unit"
        PipelineH5Output(out).write_many(preprocess_gap_spectrograms(spectra, native_zones))
        dataset = out[schema.spectrograms_masked]
        assert dataset.name == "/Processing/Spectrograms/Gaps/Raw"
        assert dataset.shape == (1, 1)
        assert dataset.dtype == np.float32
        assert dataset[0, 0] == 6
        assert list(dataset.attrs["dimDesc"]) == ["time", "freq"]
        assert dataset.attrs["unit"] == "power_unit"
        assert dataset.attrs["gap_pixel_count"] == 2
        assert dataset.attrs["mask_reduction"] == "all_true_in_2x2_block"
        assert list(dataset.attrs["mask_application_shape"]) == [3, 2]
        assert dataset.attrs["spectrogram_spatial_transform"] == "flip_D2_then_swap_D2_D3"


@pytest.mark.parametrize("reduced_only", [False, True])
def test_empty_native_or_reduced_support_skips_analysis_without_reading_source(reduced_only):
    zones = np.zeros((4, 6), dtype=bool)
    if reduced_only:
        zones[0, 0] = True
    source = FrameSource(np.ones((2, 3, 2, 3), dtype=np.float32))
    with patch("pipelines.velocity.gaps.preprocessing.Logger.log_warning") as warning:
        assert preprocess_gap_spectrograms(source, zones) == {}
    assert source.reads == []
    warning.assert_called_once()
    assert "Skipping gap analysis" in warning.call_args.args[0]
    assert "empty" in warning.call_args.args[0]


def test_missing_spectra_omit_output_with_warning():
    with patch("pipelines.velocity.gaps.preprocessing.Logger.log_warning") as warning:
        assert preprocess_gap_spectrograms(None, np.ones((2, 3), bool)) == {}
    warning.assert_called_once()
    assert "absent" in warning.call_args.args[0]
