"""DopplerView source adapter for segmentation and analysis settings."""

from __future__ import annotations

import numpy as np

from .base import SourceFileLayout, TypedSource
from .source_data import OpticDiscMeasurements

DV_CONFIG_DIR_NAME = "json"
DV_CONFIG_FILENAME = "DV_params.json"

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

    def optic_disc_measurements(self) -> OpticDiscMeasurements:
        """Read DopplerView's optic-disc measurements without interpreting geometry."""

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
