"""Lazy HD spectograms input and preservation of the raw exported axes."""

import json
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import RawH5SourceReader
from input_output.schema import (
    HD_SPECTOGRAMS_PATH,
    DopplerViewSource,
    HolodopplerSource,
)
from pipelines.vessel_inputs import load_retinal_source_data


def _hd_source(hd):
    return HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))


def test_reads_exact_root_dataset_lazily_and_preserves_exported_axes(tmp_path):
    values = np.arange(3 * 4 * 5 * 7, dtype=np.float32).reshape(3, 4, 5, 7)
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        exported = hd.create_dataset("spectrograms", data=values)
        with patch.object(RawH5SourceReader, "array", side_effect=AssertionError("eager read")):
            dataset = _hd_source(hd).spectograms_dataset()
        assert isinstance(dataset, h5py.Dataset)
        assert dataset.name == "/spectrograms"
        assert dataset.id == exported.id
        assert dataset.shape == (3, 4, 5, 7)
        np.testing.assert_array_equal(dataset[1:2, 1:3], values[1:2, 1:3])


def test_optional_source_is_none_without_the_exact_dataset(tmp_path):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        hd.create_dataset("band_2_0_6000", shape=(3, 5, 7), dtype=np.float32)
        assert _hd_source(hd).optional_spectograms_dataset() is None
        with pytest.raises(KeyError, match="spectrograms"):
            _hd_source(hd).spectograms_dataset()


def test_optional_source_ignores_legacy_spectrograms_group(tmp_path):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        group = hd.create_group("spectrograms")
        group.create_group("longtimes")
        group.create_group("singlebeat")
        group.create_dataset("spectral_cube_parameters", data="{}")

        assert _hd_source(hd).optional_spectograms_dataset() is None
        with pytest.raises(KeyError, match="spectrograms"):
            _hd_source(hd).spectograms_dataset()


@pytest.mark.parametrize("shape", [(3, 4), (3, 5, 7), (1, 3, 4, 5, 7)])
def test_rejects_incorrect_rank_with_the_external_path(tmp_path, shape):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        hd.create_dataset(HD_SPECTOGRAMS_PATH, shape=shape, dtype=np.float32)
        with pytest.raises(ValueError, match="/spectrograms.*4-D"):
            _hd_source(hd).optional_spectograms_dataset()


@pytest.mark.parametrize("dtype", ["S4", np.complex64, np.bool_])
def test_requires_real_numeric_band_power(tmp_path, dtype):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        hd.create_dataset(HD_SPECTOGRAMS_PATH, shape=(3, 4, 5, 7), dtype=dtype)
        with pytest.raises(TypeError, match="/spectrograms.*real numeric"):
            _hd_source(hd).spectograms_dataset()


@pytest.mark.parametrize("method", ["doppler_moments", "frequency_bands"])
@pytest.mark.parametrize("spectral_shape", [None, "legacy_group", (3, 4, 5, 7), (3, 4, 7, 5)])
def test_workflows_carry_optional_spectograms_without_aligning_the_volume(
    tmp_path, method, spectral_shape
):
    image_shape = (5, 7)
    with (
        h5py.File(tmp_path / "hd.h5", "w") as hd,
        h5py.File(tmp_path / "dv.h5", "w") as dv,
    ):
        for name in ("moment0", "moment2", "band_0_3000_9000", "band_1_9000_18000"):
            hd.create_dataset(name, data=np.ones((3, *image_shape), dtype=np.float32))
        hd.create_dataset(
            "HD_parameters",
            data=json.dumps({
                "pixel_pitch": [0.0002, 0.0002],
                "sampling_freq": 20,
                "batch_stride": 1,
            }),
        )
        # Swapped DV geometry must not cause a transpose of the HD spectograms.
        for name in ("artery_mask", "vein_mask"):
            dv.create_dataset(f"segmentation/Retina/{name}", data=np.zeros(image_shape[::-1]))
        if spectral_shape == "legacy_group":
            hd.create_group(HD_SPECTOGRAMS_PATH).create_group("longtimes")
        elif spectral_shape is not None:
            hd.create_dataset(HD_SPECTOGRAMS_PATH, shape=spectral_shape, dtype=np.float32)

        source = load_retinal_source_data(
            _hd_source(hd),
            DopplerViewSource(RawH5SourceReader(h5file=dv, label="DV")),
            velocity_estimation_method=method,
        )

        assert source.doppler_view.spatial_axes_swapped_to_match_hd
        if spectral_shape is None or spectral_shape == "legacy_group":
            assert source.image_maps.spectograms is None
        else:
            dataset = source.image_maps.spectograms
            assert isinstance(dataset, h5py.Dataset)
            assert dataset.id == hd[HD_SPECTOGRAMS_PATH].id
            assert dataset.shape == spectral_shape
