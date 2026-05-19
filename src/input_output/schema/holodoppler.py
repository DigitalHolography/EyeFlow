"""Holodoppler source adapter for exported HDF5 and config values."""

from __future__ import annotations

import numpy as np

from .base import HolodopplerTiming, SourceFileLayout, TypedSource

HD_CONFIG_DIR_NAME = "json"
HD_CONFIG_FILENAME = "parameters.json"
HD_MOMENT0_PATH = "moment0"
HD_MOMENT2_PATH = "moment2"
HD_OUTPUT_PASSTHROUGH_PATHS = ("registration", "zernike_coefs_radians")
HD_SAMPLING_FREQ_KEY = "sampling_freq"
HD_BATCH_STRIDE_KEY = "batch_stride"

HOLODOPPLER_LAYOUT = SourceFileLayout(
    label="HD",
    companion_suffix="HD",
    h5_folder_name="h5",
    h5_filename_template="{folder}_output.h5",
    config_dir_name=HD_CONFIG_DIR_NAME,
    config_filename=HD_CONFIG_FILENAME,
)


class HolodopplerSource(TypedSource):
    """Typed access to the Holodoppler HDF5 file and sidecar config."""

    layout = HOLODOPPLER_LAYOUT

    @classmethod
    def from_context(cls, ctx) -> "HolodopplerSource":
        return cls(ctx.sources.hd, ctx.hd_config)

    def moment0(self) -> np.ndarray:
        return self._moment(HD_MOMENT0_PATH)

    def moment2(self) -> np.ndarray:
        return self._moment(HD_MOMENT2_PATH)

    def timing(self) -> HolodopplerTiming:
        sampling_freq = self._scalar_h5_or_config(
            HD_SAMPLING_FREQ_KEY,
            HD_SAMPLING_FREQ_KEY,
        )
        batch_stride = self._scalar_h5_or_config(
            HD_BATCH_STRIDE_KEY,
            HD_BATCH_STRIDE_KEY,
        )
        if sampling_freq is None or batch_stride is None:
            raise KeyError(
                "Could not resolve Holodoppler timing from HD HDF5 or config."
            )
        return HolodopplerTiming(float(sampling_freq), float(batch_stride))

    def _moment(self, path: str) -> np.ndarray:
        squeezed = np.squeeze(np.asarray(self._array(path, dtype=np.float32)))
        if squeezed.ndim != 3:
            raise ValueError(
                "Holodoppler moment datasets must become 3-D after squeeze, "
                f"got shape {squeezed.shape}."
            )
        return np.transpose(squeezed, (0, 2, 1))
