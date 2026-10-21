"""Frequency calibration for contiguous averages of signed FFT samples."""

from __future__ import annotations

import numpy as np


def frequency_bin_spacing_khz(
    bin_count: int, sampling_frequency_hz: float, *,
    fft_bin_count: int | None = None, binning_factor: int = 1,
    frequency_start: float = 0.0, frequency_end: float = 0.5,
) -> float:
    """Return center spacing; an unspecified FFT size describes unbinned input."""
    if bin_count < 1 or not np.isfinite(sampling_frequency_hz) or sampling_frequency_hz <= 0:
        raise ValueError("Frequency calibration requires bins and a positive finite sampling frequency.")
    fft_bin_count = bin_count if fft_bin_count is None else fft_bin_count
    for count in (bin_count, fft_bin_count, binning_factor):
        if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 1:
            raise ValueError("Frequency bin counts and binning factor must be positive integers.")
    selected = _selected_fft_indices(fft_bin_count, frequency_start, frequency_end)
    if bin_count * binning_factor != selected.size:
        raise ValueError("Averaged frequency bins must cover the selected original FFT samples.")
    return float(sampling_frequency_hz * binning_factor / (fft_bin_count * 1000))


def frequency_coordinates_khz(
    bin_count: int, sampling_frequency_hz: float, *,
    fft_bin_count: int | None = None, binning_factor: int = 1,
    frequency_start: float = 0.0, frequency_end: float = 0.5,
) -> np.ndarray:
    """Use NumPy's FFT frequencies, then apply source selection and averaging.

    Select samples [:int(FFT_size*end)] and [int(FFT_size*start)+FFT_size//2:].
    Bounds (0, 0.5) retain the full FFT. Average signed frequencies exactly as
    their associated powers are averaged, including groups crossing a sign edge.
    With factor one this also describes an ordinary unbinned signed spectrum.
    """
    frequency_bin_spacing_khz(
        bin_count, sampling_frequency_hz,
        fft_bin_count=fft_bin_count, binning_factor=binning_factor,
        frequency_start=frequency_start, frequency_end=frequency_end,
    )
    fft_bin_count = bin_count if fft_bin_count is None else fft_bin_count
    raw = np.fft.fftfreq(fft_bin_count, d=1 / sampling_frequency_hz) / 1000
    indices = _selected_fft_indices(fft_bin_count, frequency_start, frequency_end)
    return raw[indices].reshape(bin_count, binning_factor).mean(axis=1)


def positive_frequency_coordinates_khz(
    bin_count: int, sampling_frequency_hz: float, *,
    fft_bin_count: int | None = None, binning_factor: int = 1,
    frequency_start: float = 0.0, frequency_end: float = 0.5,
) -> np.ndarray:
    """Return centers in the positive half, including the average containing DC."""
    frequencies = frequency_coordinates_khz(
        bin_count, sampling_frequency_hz,
        fft_bin_count=fft_bin_count, binning_factor=binning_factor,
        frequency_start=frequency_start, frequency_end=frequency_end,
    )
    return frequencies[frequencies >= 0]


def _selected_fft_indices(fft_bin_count: int, start: float, end: float) -> np.ndarray:
    if not np.all(np.isfinite([start, end])) or not (0 <= start <= 1 and 0 <= end <= 1):
        raise ValueError("Frequency selection fractions must be finite and within [0, 1].")
    bins = np.arange(fft_bin_count)
    return np.r_[bins[:int(fft_bin_count * end)], bins[int(fft_bin_count * start) + fft_bin_count // 2:]]
