"""DopplerView source adapter for segmentation and analysis settings."""

from __future__ import annotations

import numpy as np

from .base import SourceFileLayout, TypedSource

DV_CONFIG_DIR_NAME = "config"
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

    layout = DOPPLER_VIEW_LAYOUT

    @classmethod
    def from_context(cls, ctx) -> "DopplerViewSource":
        return cls(ctx.sources.dv, ctx.dv_config)

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

    def optic_disc_center(self) -> np.ndarray | None:
        return self._array(
            "segmentation/OpticDisc/optic_disc_center",
            dtype=np.float32,
            default=None,
        )

    def optic_disc_width(self) -> np.float32 | None:
        return self._array(
            "segmentation/OpticDisc/optic_disc_width",
            dtype=np.float32,
            default=None,
        )

    def optic_disc_height(self) -> np.float32 | None:
        return self._array(
            "segmentation/OpticDisc/optic_disc_height",
            dtype=np.float32,
            default=None,
        )

    def local_background_dist(self) -> int:
        return int(
            self._config_value("VelocityEstimation", "LocalBackgroundDist", 2)
        )

    def peripapillary_inner_radius(self) -> np.float32:
        return np.float32(
            self._config_value("PeripapillaryVascularZone", "InnerRadius", 0.10)
        )

    def peripapillary_outer_radius(self) -> np.float32:
        return np.float32(
            self._config_value("PeripapillaryVascularZone", "OuterRadius", 0.35)
        )

    def peripapillary_ring_count(self) -> np.int32:
        return np.int32(
            self._config_value("PeripapillaryRingAnalysis", "RingsNumber", 10)
        )

    def peripapillary_ring_width(self) -> np.float32:
        return np.float32(
            self._config_value("PeripapillaryRingAnalysis", "RingsWidth", -1.0)
        )
