"""Adapt shared frequency calibration to spectral metadata, output, and plot ticks."""

from __future__ import annotations

import numpy as np

from calculations.math.frequency import frequency_coordinates_khz
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import HD_SPECTOGRAMS_PATH
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


def write_frequency_bin_table(output, bin_count: int) -> None:
    """Persist every source bin and its signed kHz frequency once per run."""
    sampling_frequency = read_sampling_frequency_hz(output)
    if sampling_frequency is None:
        Logger.log_warning("Skipping FrequencyBinTable: a positive HD sampling frequency is required.")
        return
    schema = EyeFlowOutputPaths.active()
    frequencies = frequency_coordinates_khz(bin_count, sampling_frequency)
    output.h5.write(
        schema.spectrograms_frequency_bin_table,
        np.column_stack((np.arange(bin_count), frequencies)),
        dimDesc=["freq", "column"],
        column_names=["bin_index", "frequency_khz"],
        unit=["1", "kHz"],
        index_base=0,
        source_path=f"/{HD_SPECTOGRAMS_PATH}",
        source_bin_count=bin_count,
        sampling_frequency_hz=sampling_frequency,
        sampling_frequency_source_path="/" + schema.sampling_frequency_hz,
        frequency_bin_order="positive_then_negative_DC_zero",
        bin_domain_split_index=bin_count // 2,
        definition="zero-based HD bin index and signed frequency in kHz with full-spectrum fs/N spacing",
    )


def frequency_ticks(bin_count: int, sampling_frequency_hz: float) -> tuple[np.ndarray, list[str]]:
    """Label the stored bin positions in signed kHz without reordering arrays."""
    positions = np.unique(np.r_[np.arange(0, bin_count, max(1, int(np.ceil(bin_count / 8)))),
                                bin_count // 2, bin_count - 1])
    frequencies = frequency_coordinates_khz(bin_count, sampling_frequency_hz)
    return positions, [f"{value:.3g}" for value in frequencies[positions]]
