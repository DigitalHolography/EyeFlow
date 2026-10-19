"""DC-zero frequency calibration for full-spectrum HD exports."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths


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


def frequency_bin_spacing_khz(bin_count: int, sampling_frequency_hz: float) -> float:
    """Return the positive full-spectrum bin spacing in kHz."""
    if bin_count < 1 or not np.isfinite(sampling_frequency_hz) or sampling_frequency_hz <= 0:
        raise ValueError("Frequency calibration requires bins and a positive finite sampling frequency.")
    return float(sampling_frequency_hz / (bin_count * 1000))


def positive_frequency_coordinates_khz(bin_count: int, sampling_frequency_hz: float) -> np.ndarray:
    """Return the actual positive/DC half of the full-spectrum frequency grid."""
    return frequency_coordinates_khz(bin_count, sampling_frequency_hz)[:bin_count // 2]


def frequency_coordinates_khz(bin_count: int, sampling_frequency_hz: float) -> np.ndarray:
    """Map HD full-spectrum positive-then-negative bins with spacing fs / N.

    For 64 bins, indices 0..31 are positive/DC and index 32 is -fs/2.
    The corresponding positive Nyquist edge +fs/2 is just beyond bin 31.
    """
    spacing = frequency_bin_spacing_khz(bin_count, sampling_frequency_hz)
    bins = np.arange(bin_count, dtype=np.float64)
    if bin_count > 1:
        bins[bin_count // 2:] -= bin_count
    return bins * spacing


def frequency_ticks(bin_count: int, sampling_frequency_hz: float) -> tuple[np.ndarray, list[str]]:
    """Label the stored bin positions in signed kHz without reordering arrays."""
    positions = np.unique(np.r_[np.arange(0, bin_count, max(1, int(np.ceil(bin_count / 8)))),
                                bin_count // 2, bin_count - 1])
    frequencies = frequency_coordinates_khz(bin_count, sampling_frequency_hz)
    return positions, [f"{value:.3g}" for value in frequencies[positions]]
