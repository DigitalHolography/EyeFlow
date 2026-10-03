"""Holodoppler source adapter for exported HDF5 and config values."""

from __future__ import annotations

import json
from collections.abc import Mapping

import numpy as np

from .base import SourceFileLayout, TypedSource
from .source_data import HolodopplerMetadata, HolodopplerTiming, PixelPitch

HD_CONFIG_DIR_NAME = "json"
HD_CONFIG_FILENAME = "parameters.json"
HD_MOMENT0_PATH = "moment0"
HD_MOMENT2_PATH = "moment2"
HD_BAND_LF_PATH = "band_0_3000_9000"
HD_BAND_HF_PATH = "band_1_9000_18000"
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

    @classmethod
    def from_context(cls, ctx) -> HolodopplerSource:
        return cls(ctx.inputs.hd.h5, ctx.inputs.hd.config)

    def moment0_dataset(self):
        return self._moment_dataset(HD_MOMENT0_PATHS)

    def moment2_dataset(self):
        return self._moment_dataset(HD_MOMENT2_PATHS)

    def optional_moment0_dataset(self):
        """Return the raw zeroth moment when it is present."""

        return self._optional_moment_dataset(HD_MOMENT0_PATHS)

    def optional_moment2_dataset(self):
        """Return the raw second moment when it is present."""

        return self._optional_moment_dataset(HD_MOMENT2_PATHS)

    def frequency_band_datasets(self):
        """Return the exact low/high PSD bands used by the ratio estimator.

        Band discovery is deliberately not heuristic: changing the configured
        frequency limits changes the meaning of the ratio, so this first
        implementation accepts only HoloDoppler's standard two-band export.
        """

        paths = (HD_BAND_LF_PATH, HD_BAND_HF_PATH)
        missing = [f"/{path}" for path in paths if path not in self._reader]
        if missing:
            source = str(self.filename or "<unknown HD source>")
            raise KeyError(
                "velocity_estimation_method='frequency_bands' requires "
                "HoloDoppler datasets "
                f"{', '.join(missing)}; missing from HD source file {source!r}."
            )

        low_frequency = self._frequency_band_dataset(HD_BAND_LF_PATH)
        high_frequency = self._frequency_band_dataset(HD_BAND_HF_PATH)
        if tuple(low_frequency.shape) != tuple(high_frequency.shape):
            raise ValueError(
                "HoloDoppler frequency-band datasets must have identical "
                f"(frame, y, x) shapes; /{HD_BAND_LF_PATH} has shape "
                f"{low_frequency.shape} and /{HD_BAND_HF_PATH} has shape "
                f"{high_frequency.shape}."
            )
        return low_frequency, high_frequency

    def moment0_flat_field_dataset(self):
        """Return a precomputed flat-field moment, when exported by Holodoppler."""
        return self._optional_moment_dataset(HD_MOMENT0_FLAT_FIELD_PATHS)

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
        path = self._first_path(paths)
        if path is None:
            raise KeyError(
                "Missing Holodoppler moment dataset. Tried: "
                + ", ".join(repr(candidate) for candidate in paths)
            )
        dataset = self._dataset(path)
        if dataset.ndim != 3:
            raise ValueError(
                "Holodoppler moment datasets must be 3-D for lazy processing, "
                f"got shape {dataset.shape}."
            )
        return dataset

    def _optional_moment_dataset(self, paths: tuple[str, ...]):
        path = self._first_path(paths)
        if path is None:
            return None
        dataset = self._dataset(path)
        if dataset.ndim != 3:
            raise ValueError(
                "Holodoppler flat-field moment datasets must be 3-D for lazy "
                f"processing, got shape {dataset.shape}."
            )
        return dataset

    def _frequency_band_dataset(self, path: str):
        dataset = self._dataset(path)
        if dataset.ndim != 3:
            raise ValueError(
                f"HoloDoppler frequency-band dataset /{path} must be a "
                f"3-D (frame, y, x) array, got shape {dataset.shape}."
            )
        if not np.issubdtype(dataset.dtype, np.number):
            raise TypeError(
                f"HoloDoppler frequency-band dataset /{path} must be "
                f"numeric, got dtype {dataset.dtype}."
            )
        return dataset

    def _first_path(self, paths: tuple[str, ...]) -> str | None:
        return next((candidate for candidate in paths if candidate in self._reader), None)
