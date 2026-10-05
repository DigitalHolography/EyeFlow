"""Holodoppler source adapter for exported HDF5 and config values."""

from __future__ import annotations

import json
from collections.abc import Mapping

import numpy as np
import h5py


from .base import SourceFileLayout, TypedSource
from .source_data import HolodopplerMetadata, HolodopplerTiming, PixelPitch

HD_CONFIG_DIR_NAME = "json"
HD_CONFIG_FILENAME = "parameters.json"
HD_MOMENT0_PATH = "moment0"
HD_MOMENT2_PATH = "moment2"
HD_MOMENT0_PATHS = (HD_MOMENT0_PATH, "M0")
HD_MOMENT2_PATHS = (HD_MOMENT2_PATH, "M2")
HD_MOMENT0_FLAT_FIELD_PATHS = ("moment0ff", "M0FF")
HD_OUTPUT_PASSTHROUGH_PATHS = (
    ("registration", "Meta/registration"),
    ("zernike_coefs_radians", "zernike_coefs_radians"),
)
HD_SAMPLING_FREQ_KEY = "sampling_freq"
HD_BATCH_STRIDE_KEY = "batch_stride"
HD_PARAMETERS_KEY = "HD_parameters"

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

    def moment0_dataset(self):
        return self._moment_dataset(HD_MOMENT0_PATHS)

    def named_root_moment_dataset(self, moment_path: str = HD_MOMENT0_PATH):
        """Resolve a selected root moment, including the legacy M0 alias."""
        from ..writers.h5 import normalize_h5_path

        normalized = normalize_h5_path(moment_path)
        if not normalized or "/" in normalized:
            raise ValueError("The HoloDoppler moment must be a root dataset name.")
        candidates = HD_MOMENT0_PATHS if normalized == HD_MOMENT0_PATH else (normalized,)
        for candidate in candidates:
            found = self._reader.get(candidate)
            if found is None:
                continue
            if not isinstance(found, h5py.Dataset) or found.ndim != 3:
                raise ValueError(
                    f"HoloDoppler moment '{candidate}' must be a 3-D dataset, "
                    f"got {getattr(found, 'shape', None)}."
                )
            return found
        raise KeyError(
            "Missing HoloDoppler root moment dataset. Tried: "
            + ", ".join(repr(candidate) for candidate in candidates)
        )

    def moment2_dataset(self):
        return self._moment_dataset(HD_MOMENT2_PATHS)

    def moment0_flat_field_dataset(self):
        """Return a precomputed flat-field moment, when exported by Holodoppler."""
        return self._optional_moment_dataset(
            HD_MOMENT0_FLAT_FIELD_PATHS, description="flat-field moment"
        )

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
            raise KeyError("Could not resolve Holodoppler timing from HD HDF5 or config.")
        return HolodopplerTiming(float(sampling_freq), float(batch_stride))

    def pixel_pitch(self) -> PixelPitch:
        """Return the native ``(x, y)`` pixel pitch from ``HD_parameters``."""

        raw = self._value(HD_PARAMETERS_KEY, default=None)
        if raw is None:
            raise KeyError("Missing Holodoppler dataset 'HD_parameters'.")
        if isinstance(raw, np.ndarray):
            if raw.size != 1:
                raise ValueError("HD_parameters must be a scalar JSON value.")
            raw = raw.reshape(()).item()
        if isinstance(raw, (bytes, np.bytes_)):
            raw = raw.decode("utf-8")
        if isinstance(raw, str):
            try:
                parameters = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise ValueError("HD_parameters must contain valid JSON.") from exc
        else:
            parameters = raw
        if not isinstance(parameters, Mapping):
            raise TypeError("HD_parameters must decode to a dictionary.")
        if "pixel_pitch" not in parameters:
            raise KeyError("HD_parameters does not contain 'pixel_pitch'.")
        values = np.asarray(parameters["pixel_pitch"], dtype=np.float64).reshape(-1)
        if values.size != 2:
            raise ValueError("HD_parameters['pixel_pitch'] must contain exactly two (x, y) values.")
        return PixelPitch(values[0], values[1])

    def metadata(self) -> HolodopplerMetadata:
        """Return all Holodoppler acquisition metadata used by EyeFlow."""

        return HolodopplerMetadata(
            timing=self.timing(),
            pixel_pitch=self.pixel_pitch(),
        )

    def _moment_dataset(self, paths: tuple[str, ...]):
        dataset = self._optional_moment_dataset(paths)
        if dataset is not None:
            return dataset
        raise KeyError(
            "Missing Holodoppler moment dataset. Tried: "
            + ", ".join(repr(candidate) for candidate in paths)
        )

    def _optional_moment_dataset(
        self, paths: tuple[str, ...], *, description: str = "moment"
    ):
        path = self._first_path(paths)
        if path is None:
            return None
        dataset = self._dataset(path)
        if dataset.ndim != 3:
            raise ValueError(
                f"Holodoppler {description} datasets must be 3-D for lazy processing, "
                f"got shape {dataset.shape}."
            )
        return dataset

    def _first_path(self, paths: tuple[str, ...]) -> str | None:
        return next((candidate for candidate in paths if candidate in self._reader), None)
