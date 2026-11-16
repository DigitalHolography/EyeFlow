"""Holodoppler source adapter for exported HDF5 and config values."""

from __future__ import annotations

import json
import re
from collections.abc import Mapping

from input_output.readers.nodes import is_dataset
import numpy as np

from .base import SourceFileLayout, TypedSource, scalar_from_value
from retinal_source_data import HolodopplerMetadata, HolodopplerTiming, PixelPitch

HD_CONFIG_DIR_NAME = "json"
HD_CONFIG_FILENAME = "parameters.json"
HD_MOMENT0_PATH = "moment0"
HD_MOMENT2_PATH = "moment2"
HD_BAND_LF_PATH = "band_0_3000_9000"
HD_BAND_HF_PATH = "band_1_9000_18000"
HD_BAND_NAME_PATTERN = re.compile(r"band_[0-9]+_([0-9]+)_([0-9]+)")
HD_SPECTOGRAMS_PATH = "spectrograms"
HD_MOMENT0_PATHS = (HD_MOMENT0_PATH, "M0")
HD_MOMENT2_PATHS = (HD_MOMENT2_PATH, "M2")
HD_MOMENT0_FLAT_FIELD_PATHS = ("moment0ff", "M0FF")
HD_MOMENT_PATHS = {
    HD_MOMENT0_PATH: HD_MOMENT0_PATHS,
    HD_MOMENT2_PATH: HD_MOMENT2_PATHS,
    "moment0ff": HD_MOMENT0_FLAT_FIELD_PATHS,
}
HD_OUTPUT_PASSTHROUGH_PATHS = (
    ("registration", "Meta/registration"),
    ("zernike_coefs_radians", "zernike_coefs_radians"),
)
HD_SAMPLING_FREQ_KEY = "sampling_freq"
HD_BATCH_STRIDE_KEY = "batch_stride"
HD_PARAMETERS_KEY = "HD_parameters"
HD_BATCH_SIZE_KEY = "batch_size"
HD_SPECTROGRAM_FREQUENCY_START_KEY = "spectrogram_frequency_start"
HD_SPECTROGRAM_FREQUENCY_END_KEY = "spectrogram_frequency_end"
HD_SPECTROGRAM_FREQUENCY_BINNING_KEY = "spectrogram_frequency_binning"

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

    def __init__(self, reader, config=None, *, cache_parameters: bool = False) -> None:
        super().__init__(reader, config)
        self._cache_parameters = cache_parameters

    def moment(self, name: str = HD_MOMENT0_PATH, *, required: bool = False):
        """Return a lazy root moment, resolving canonical and legacy names."""
        name = str(name).replace("\\", "/").strip("/")
        if not name or "/" in name:
            raise ValueError("The HoloDoppler moment must be a root dataset name.")
        paths = HD_MOMENT_PATHS.get(name, (name,))
        return self._reader.dataset(paths, required=required, ndim=3)

    def moment0(self, *, required: bool = False):
        return self.moment(HD_MOMENT0_PATH, required=required)

    def moment2(self, *, required: bool = False):
        return self.moment(HD_MOMENT2_PATH, required=required)

    def moment0ff(self, *, required: bool = False):
        return self.moment("moment0ff", required=required)

    def spectrograms(self, *, required: bool = False):
        """Return lazy (time, freq, y, x) power; legacy groups are unavailable."""
        return self._reader.dataset(HD_SPECTOGRAMS_PATH, required=required, ndim=4,
                            kinds="iuf", ignore_groups=True)

    def frequency_bands(self):
        """Return the first/last PSD videos sorted by numeric frequency bounds."""

        bands = []
        for path in self._reader.keys():
            match = HD_BAND_NAME_PATTERN.fullmatch(path)
            if match is None:
                continue
            if not is_dataset(self._reader.get(path)):
                continue
            bands.append((int(match[1]), int(match[2]), path))
        bands.sort()
        if len(bands) < 2:
            source = str(self.filename or "<unknown HD source>")
            raise KeyError(
                "velocity_estimation_method='frequency_bands' requires at least "
                "two HoloDoppler datasets named band_<i>_<n>_<m>; "
                f"found {len(bands)} in HD source file {source!r}."
            )

        low_path, high_path = bands[0][2], bands[-1][2]
        low_frequency = self._reader.dataset(low_path, required=True, ndim=3, kinds="iufc")
        high_frequency = self._reader.dataset(high_path, required=True, ndim=3, kinds="iufc")
        if tuple(low_frequency.shape) != tuple(high_frequency.shape):
            raise ValueError(
                "HoloDoppler frequency-band datasets must have identical "
                f"(frame, y, x) shapes; /{low_path} has shape "
                f"{low_frequency.shape} and /{high_path} has shape "
                f"{high_frequency.shape}."
            )
        return low_frequency, high_frequency

    def sampling_frequency_hz(self) -> float:
        """Resolve positive acquisition sampling frequency without requiring a frame stride."""
        parameters = self._parameters() or {}
        frequency = scalar_from_value(parameters.get(HD_SAMPLING_FREQ_KEY))
        if frequency is None:
            frequency = scalar_from_value(self._config.get(HD_SAMPLING_FREQ_KEY))
        if frequency is None or not np.isfinite(frequency) or frequency <= 0:
            raise ValueError("HoloDoppler sampling_freq must be a positive finite frequency in Hz.")
        return float(frequency)

    def timing(self) -> HolodopplerTiming:
        parameters = self._parameters() or {}
        sampling_freq = scalar_from_value(parameters.get(HD_SAMPLING_FREQ_KEY))
        batch_stride = scalar_from_value(parameters.get(HD_BATCH_STRIDE_KEY))
        if sampling_freq is None:
            sampling_freq = scalar_from_value(self._config.get(HD_SAMPLING_FREQ_KEY))
        if batch_stride is None:
            batch_stride = scalar_from_value(self._config.get(HD_BATCH_STRIDE_KEY))
        if sampling_freq is None or batch_stride is None:
            raise KeyError("Could not resolve Holodoppler timing from HD HDF5 or config.")
        return HolodopplerTiming(float(sampling_freq), float(batch_stride))

    def spectrogram_frequency_parameters(self, bin_count: int) -> dict[str, object]:
        """Resolve the original FFT size and validate the 4-D HD export layout.

        The FFT processes one full batch, independently of registration batches
        and Shack-Hartmann reference windows. Mirror main_default's frequency
        selection before averaging; unresolved HD bounds mean the full spectrum.
        """
        parameters = self._parameters() or {}
        batch_size = scalar_from_value(parameters.get(HD_BATCH_SIZE_KEY))
        if batch_size is None:
            batch_size = scalar_from_value(self._config.get(HD_BATCH_SIZE_KEY))
        if (
            isinstance(batch_size, (bool, np.bool_))
            or not isinstance(batch_size, (int, float, np.integer, np.floating))
            or not np.isfinite(batch_size) or batch_size < 2 or batch_size != int(batch_size)
        ):
            raise ValueError("HoloDoppler spectrogram batch_size must be an integer of at least two.")
        fft_bin_count = int(batch_size)
        start = scalar_from_value(parameters.get(HD_SPECTROGRAM_FREQUENCY_START_KEY))
        end = scalar_from_value(parameters.get(HD_SPECTROGRAM_FREQUENCY_END_KEY))
        bounds_available = all(
            isinstance(value, (int, float, np.integer, np.floating))
            and not isinstance(value, (bool, np.bool_))
            and np.isfinite(value) and 0 <= value <= 1
            for value in (start, end)
        )
        if not bounds_available:
            start, end = 0.0, 0.5
        start, end = float(start), float(end)
        start_index, end_index = int(fft_bin_count * start), int(fft_bin_count * end)
        selected_count = end_index + max(0, fft_bin_count - (start_index + fft_bin_count // 2))
        if bin_count < 1 or selected_count == 0 or selected_count % bin_count:
            raise ValueError(
                f"HD /spectrograms has {bin_count} bins, but batch_size={fft_bin_count} "
                f"and frequency bounds ({start}, {end}) select {selected_count} FFT samples, "
                "which must divide into whole averaged bins."
            )
        binning_factor = selected_count // bin_count
        configured_binning = scalar_from_value(parameters.get(HD_SPECTROGRAM_FREQUENCY_BINNING_KEY))
        if configured_binning is not None and (
            isinstance(configured_binning, (bool, np.bool_)) or configured_binning != binning_factor
        ):
            raise ValueError(
                f"HD /spectrograms requires averaging {binning_factor} samples per bin, "
                f"inconsistent with spectrogram_frequency_binning={configured_binning}."
            )
        return {
            "fft_bin_count": fft_bin_count,
            "frequency_binning_factor": binning_factor,
            "fft_bin_count_parameter": HD_BATCH_SIZE_KEY,
            "fft_bin_count_source_path": (
                f"/{HD_PARAMETERS_KEY}" if parameters.get(HD_BATCH_SIZE_KEY) is not None
                else HD_CONFIG_FILENAME
            ),
            "frequency_selection": "HD_spectrogram_frequency_bounds",
            HD_SPECTROGRAM_FREQUENCY_START_KEY: start,
            HD_SPECTROGRAM_FREQUENCY_END_KEY: end,
            "frequency_selection_source": (
                f"/{HD_PARAMETERS_KEY}" if bounds_available else "full_spectrum_fallback"
            ),
            "selected_fft_sample_count": selected_count,
        }

    def pixel_pitch(self) -> PixelPitch:
        """Return the native ``(x, y)`` pixel pitch from ``HD_parameters``."""

        parameters = self._parameters()
        if parameters is None:
            raise KeyError("Missing Holodoppler dataset 'HD_parameters'.")

        if "pixel_pitch" not in parameters:
            raise KeyError("HD_parameters does not contain 'pixel_pitch'.")
        values = np.asarray(parameters["pixel_pitch"], dtype=np.float64).reshape(-1)
        if values.size != 2:
            raise ValueError("HD_parameters['pixel_pitch'] must contain exactly two (x, y) values.")
        return PixelPitch(values[0], values[1])

    def _parameters(self) -> Mapping[str, object] | None:
        if self._cache_parameters and hasattr(self, "_cached_parameters"):
            return self._cached_parameters
        raw = self._reader.value(HD_PARAMETERS_KEY, default=None)
        if raw is None:
            if self._cache_parameters:
                self._cached_parameters = None
            return None
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
        if self._cache_parameters:
            self._cached_parameters = parameters
        return parameters

    def metadata(self) -> HolodopplerMetadata:
        """Return all Holodoppler acquisition metadata used by EyeFlow."""

        return HolodopplerMetadata(
            timing=self.timing(),
            pixel_pitch=self.pixel_pitch(),
        )

    INPUTS = {
        "metadata": lambda source: source._reader.get(HD_PARAMETERS_KEY),
        "moment0": moment0,
        "moment2": moment2,
        "moment0ff": moment0ff,
        "band_low": lambda source: source.frequency_bands()[0],
        "band_high": lambda source: source.frequency_bands()[1],
        "spectrograms": spectrograms,
    }
