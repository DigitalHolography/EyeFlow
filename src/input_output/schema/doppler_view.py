"""DopplerView source adapter for segmentation and analysis settings."""

from __future__ import annotations

import numpy as np

from calculations.topology import OpticDisc

from .base import SourceFileLayout, TypedSource

DV_CONFIG_DIR_NAME = "json"
DV_CONFIG_FILENAME = "DV_params.json"
DV_RETINAL_VESSEL_LOGITS_PATH = "segmentation/Retina/vessel_segmentation_logits"
DV_RETINAL_VESSEL_LOGITS_THRESHOLD = 0.1
DV_CHOROID_VESSEL_MASK_PATH = "segmentation/Choroid/vessel_mask"

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

    def retinal_artery_mask(self) -> np.ndarray:
        return self._array("segmentation/Retina/artery_mask", dtype=bool)

    def retinal_vein_mask(self) -> np.ndarray:
        return self._array("segmentation/Retina/vein_mask", dtype=bool)

    def retinal_labeled_vessels(self) -> np.ndarray | None:
        return self._array(
            "segmentation/Retina/labeled_vessels",
            dtype=np.int32,
            default=None,
        )

    def retinal_vessel_mask(self) -> np.ndarray | None:
        """Return all detected retinal vessels, binarizing raw logits at > 0.1."""
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

    def optic_disc(
        self,
        image_shape: tuple[int, int],
        *,
        pixel_pitch_mm: float | None = None,
    ) -> OpticDisc:
        """Return DopplerView's optic-disc measurements in its native frame."""

        return OpticDisc.from_measurements(
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
            image_shape=image_shape,
            pixel_pitch_mm=pixel_pitch_mm,
        )

    def local_background_dist(self) -> int:
        return int(self._config_value("VelocityEstimation", "LocalBackgroundDist", 2))
