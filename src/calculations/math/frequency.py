"""Full-spectrum frequency calibration independent of pipeline and output policy."""

from __future__ import annotations

import numpy as np


def frequency_bin_spacing_khz(bin_count: int, sampling_frequency_hz: float) -> float:
    """Return the positive full-spectrum bin spacing in kHz."""
    if bin_count < 1 or not np.isfinite(sampling_frequency_hz) or sampling_frequency_hz <= 0:
        raise ValueError("Frequency calibration requires bins and a positive finite sampling frequency.")
    return float(sampling_frequency_hz / (bin_count * 1000))


def frequency_coordinates_khz(bin_count: int, sampling_frequency_hz: float) -> np.ndarray:
    """Map positive-then-negative bins with spacing fs/N and a floor(N/2) split.

    For 64 bins, indices 0..31 are positive/DC and index 32 is -fs/2.
    A single-bin spectrum is DC. Preserve source order, including odd domains.
    """
    spacing = frequency_bin_spacing_khz(bin_count, sampling_frequency_hz)
    bins = np.arange(bin_count, dtype=np.float64)
    if bin_count > 1:
        bins[bin_count // 2:] -= bin_count
    return bins * spacing


def positive_frequency_coordinates_khz(bin_count: int, sampling_frequency_hz: float) -> np.ndarray:
    """Return the actual positive/DC half of the full-spectrum frequency grid."""
    return frequency_coordinates_khz(bin_count, sampling_frequency_hz)[:bin_count // 2]
