"""Requested DC-zero frequency calibration for exported spectral bins."""

from __future__ import annotations

import numpy as np

POSITIVE_FREQUENCY_BIN_COUNT = 32


def positive_frequency_coordinates_khz(sampling_frequency_hz: float) -> np.ndarray:
    """Return display bins 0..32 with bin 32 fixed at fs/4, independent of array size.

    Bin 32 is the calibration endpoint; the positive-spectrum figures report
    only bins through 31. These coordinates are used only for plotting.
    """
    if not np.isfinite(sampling_frequency_hz) or sampling_frequency_hz <= 0:
        raise ValueError("Frequency calibration requires a positive finite sampling frequency.")
    return np.linspace(0.0, sampling_frequency_hz / 4000, POSITIVE_FREQUENCY_BIN_COUNT + 1)


def frequency_coordinates_khz(bin_count: int, sampling_frequency_hz: float) -> np.ndarray:
    """Map HD positive-then-negative bin order with spacing fs / (2 * N).

    For 64 bins, indices 0..31 are positive/DC and index 32 is -fs/4.
    The corresponding positive edge +fs/4 is just beyond the positive half.
    """
    if bin_count < 1 or not np.isfinite(sampling_frequency_hz) or sampling_frequency_hz <= 0:
        raise ValueError("Frequency calibration requires bins and a positive finite sampling frequency.")
    bins = np.arange(bin_count, dtype=np.float64)
    if bin_count > 1:
        bins[bin_count // 2:] -= bin_count
    return bins * (sampling_frequency_hz / (2 * bin_count * 1000))


def frequency_ticks(bin_count: int, sampling_frequency_hz: float) -> tuple[np.ndarray, list[str]]:
    """Label the stored bin positions in signed kHz without reordering arrays."""
    positions = np.unique(np.r_[np.arange(0, bin_count, max(1, int(np.ceil(bin_count / 8)))),
                                bin_count // 2, bin_count - 1])
    frequencies = frequency_coordinates_khz(bin_count, sampling_frequency_hz)
    return positions, [f"{value:.3g}" for value in frequencies[positions]]
