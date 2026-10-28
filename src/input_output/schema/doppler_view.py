"""DopplerView source adapter for segmentation and analysis settings."""

from __future__ import annotations

import numpy as np
import h5py

from .base import MISSING, SourceFileLayout, TypedSource
from retinal_source_data import OpticDiscMeasurements

DV_CONFIG_DIR_NAME = "json"
DV_CONFIG_FILENAME = "DV_params.json"
DV_RETINAL_VESSEL_LOGITS_PATH = "segmentation/Retina/vessel_segmentation_logits"
DV_RETINAL_VESSEL_LOGITS_THRESHOLD = 0.1
DV_CHOROID_VESSEL_MASK_PATH = "segmentation/Choroid/vessel_mask"
DV_ARTERY_MASK_PATH = "segmentation/Retina/artery_mask"
DV_VEIN_MASK_PATH = "segmentation/Retina/vein_mask"
DV_VESSEL_MASK_PATH = "segmentation/Retina/vessel_mask"
DV_LABELED_VESSELS_PATH = "segmentation/Retina/labeled_vessels"

DOPPLER_VIEW_LAYOUT = SourceFileLayout(
    label="DV",
    companion_suffix="DV",
    h5_folder_name="h5",
    h5_filename_template="{folder}.h5",
    config_dir_name=DV_CONFIG_DIR_NAME,
    config_filename=DV_CONFIG_FILENAME,
)


class DopplerViewSource(TypedSource):
    """Typed access to the DopplerView HDF5 file and sidecar config."""

    def __init__(self, reader, config=None, *, cache_small_arrays: bool = False) -> None:
        super().__init__(reader, config)
        self._cache_small_arrays = cache_small_arrays
        self._small_array_cache: dict[str, np.ndarray | None] = {}

    def _array(self, path: str, *, dtype=None, default=MISSING):
        # DV fields used here are small masks or scalar measurements. Share the
        # source read across velocity methods while isolating callers' edits.
        if not self._cache_small_arrays:
            return super()._array(path, dtype=dtype, default=default)
        if path not in self._small_array_cache:
            self._small_array_cache[path] = super()._array(path, default=None)
        value = self._small_array_cache[path]
        if value is None:
            if default is MISSING:
                return super()._array(path, dtype=dtype)
            if default is None:
                return None
            return np.asarray(default, dtype=dtype)
        return np.array(value, dtype=dtype, copy=True)

    def retinal_mask(self, path: str, *, required: bool = False) -> np.ndarray | None:
        """Read one native retinal mask without choosing pipeline mask policy."""
        found = self._reader.get(path)
        if found is None:
            if required:
                raise KeyError(f"Missing DopplerView mask dataset at '{path}'.")
            return None
        if not isinstance(found, h5py.Dataset):
            raise ValueError(f"DopplerView mask path '{path}' is not a dataset.")
        return self._array(path)

    def retinal_artery_mask(self) -> np.ndarray:
        return self._array(DV_ARTERY_MASK_PATH, dtype=bool)

    def retinal_vein_mask(self) -> np.ndarray:
        return self._array(DV_VEIN_MASK_PATH, dtype=bool)

    def retinal_labeled_vessels(self) -> np.ndarray | None:
        return self._array(
            DV_LABELED_VESSELS_PATH,
            dtype=np.int32,
            default=None,
        )

    def retinal_vessel_mask(self) -> np.ndarray | None:
        """Return native (D0=y, D1=x) vessel support using raw logits > 0.1.

        Preserve row order: native Y is inverted relative to EF's published
        lower-left image frame. Shared input alignment handles legacy axis
        swaps; output serialization handles the Y flip and (x, y) ordering.
        """
        return self._optional_binary_mask(DV_RETINAL_VESSEL_LOGITS_PATH, logits=True)

    def choroid_vessel_mask(self) -> np.ndarray | None:
        """Return binarized choroid support when exported by DopplerView."""
        return self._optional_binary_mask(DV_CHOROID_VESSEL_MASK_PATH)

    def _optional_binary_mask(self, path: str, *, logits: bool = False):
        value = self._array(path, default=None)
        if value is None:
            return None
        value = np.asarray(value)
        if value.ndim != 2:
            raise ValueError(f"DopplerView /{path} must be 2-D, got shape {value.shape}.")
        if value.dtype.kind not in "biuf":
            raise TypeError(f"DopplerView /{path} must contain real numeric or boolean values.")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"DopplerView /{path} must contain finite values.")
        return value > DV_RETINAL_VESSEL_LOGITS_THRESHOLD if logits else value != 0

    def optic_disc_measurements(self) -> OpticDiscMeasurements:
        """Read DopplerView's optic-disc measurements in its native frame."""

        return OpticDiscMeasurements(
            mask=self._array(
                "segmentation/OpticDisc/mask",
                default=None,
            ),
            center=self._array(
                "segmentation/OpticDisc/center",
                dtype=np.float32,
                default=None,
            ),
            width=self._array(
                "segmentation/OpticDisc/width",
                dtype=np.float32,
                default=None,
            ),
            height=self._array(
                "segmentation/OpticDisc/height",
                dtype=np.float32,
                default=None,
            ),
        )

    def local_background_dist(self) -> int:
        return int(self._config_value("VelocityEstimation", "LocalBackgroundDist", 2))
