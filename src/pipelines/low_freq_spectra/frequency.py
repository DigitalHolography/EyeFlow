"""Adapt shared frequency calibration to spectral metadata, output, and plot ticks."""

from __future__ import annotations

import numpy as np

from calculations.math.frequency import frequency_bin_spacing_khz, frequency_coordinates_khz
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import (
    HD_SPECTOGRAMS_PATH, HD_SPECTROGRAM_FREQUENCY_START_KEY, HD_SPECTROGRAM_FREQUENCY_END_KEY,
)
from utils.logger import Logger


def read_sampling_frequency_hz(output) -> float | None:
    """Read shared HD calibration, returning None when it is absent or invalid."""
    source = output.h5.get(EyeFlowOutputPaths.active().sampling_frequency_hz)
    if source is None:
        return None
    try:
        value = float(np.asarray(source[()]).reshape(()))
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


def write_frequency_bin_table(output, bin_count: int, *, hd=None) -> None:
    """Persist signed kHz frequencies, indexed by zero-based source bin."""
    output.h5.write_many(pack_frequency_bin_table(output, bin_count, hd=hd))


def pack_frequency_bin_table(output, bin_count: int, *, hd=None) -> dict:
    """Pack calibration for the standard graph handoff/writer contract."""
    sampling_frequency = read_sampling_frequency_hz(output)
    if sampling_frequency is None:
        Logger.log_warning("Skipping FrequencyBinTable: a positive HD sampling frequency is required.")
        return {}
    try:
        if hd is None:
            raise ValueError("HD batch_size and spectrogram export metadata are required.")
        parameters = hd.spectrogram_frequency_parameters(bin_count)
    except (KeyError, TypeError, ValueError) as exc:
        Logger.log_warning(f"Skipping FrequencyBinTable: {exc}")
        return {}
    schema = EyeFlowOutputPaths.active()
    calibration = dict(
        fft_bin_count=parameters["fft_bin_count"],
        binning_factor=parameters["frequency_binning_factor"],
        frequency_start=parameters.get(HD_SPECTROGRAM_FREQUENCY_START_KEY, 0.0),
        frequency_end=parameters.get(HD_SPECTROGRAM_FREQUENCY_END_KEY, 0.5),
    )
    frequencies = frequency_coordinates_khz(bin_count, sampling_frequency, **calibration)
    attrs = dict(
        dimDesc=["freq"],
        unit="kHz",
        index_base=0,
        source_path=f"/{HD_SPECTOGRAMS_PATH}",
        source_bin_count=bin_count,
        sampling_frequency_hz=sampling_frequency,
        sampling_frequency_source_path="/" + schema.sampling_frequency_hz,
        **parameters,
        frequency_bin_spacing_khz=frequency_bin_spacing_khz(bin_count, sampling_frequency, **calibration),
        frequency_calibration="mean_of_contiguous_fft_samples",
        frequency_bin_order="positive_then_negative_averaged_FFT_centers",
        bin_domain_split_index=int(np.count_nonzero(frequencies >= 0)),
        definition="mean signed FFT frequency in kHz of the original samples contributing to each HD plane",
    )
    return {schema.spectrograms_frequency_bin_table: (frequencies, attrs)}


def read_frequency_bin_table(output, bin_count: int | None = None):
    """Return full-spectrum centers reconstructed from HD FFT and binning metadata."""
    source = output.h5.get(EyeFlowOutputPaths.active().spectrograms_frequency_bin_table)
    if (
        source is None or source.ndim != 1 or source.size == 0 or source.dtype.kind not in "iuf"
        or (bin_count is not None and source.size != bin_count)
    ):
        return None
    if source.attrs.get("unit") != "kHz" or not np.all(np.isfinite(source[()])):
        return None
    if source.attrs.get("frequency_calibration") != "mean_of_contiguous_fft_samples":
        return None
    if "fft_bin_count" not in source.attrs or "frequency_binning_factor" not in source.attrs:
        return None
    return source


def frequency_ticks(frequencies_khz: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Label the stored bin positions in signed kHz without reordering arrays."""
    frequencies = np.asarray(frequencies_khz)
    bin_count = frequencies.size
    negative_bins = np.flatnonzero(frequencies < 0)
    split = negative_bins[:1]
    positions = np.unique(np.r_[np.arange(0, bin_count, max(1, int(np.ceil(bin_count / 8)))),
                                split, bin_count - 1])
    return positions, [f"{value:.3g}" for value in frequencies[positions]]
