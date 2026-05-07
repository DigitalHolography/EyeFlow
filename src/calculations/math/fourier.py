"""Fourier transform helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np


def interpft_real(signal, target_length: int) -> np.ndarray:
    source = np.asarray(signal, dtype=np.float32).reshape(-1)
    source_length = int(source.size)
    if source_length == 0:
        raise ValueError("interpft requires a non-empty signal.")
    if target_length <= 0:
        raise ValueError("interpft target_length must be positive.")
    if target_length == source_length:
        return source.copy()

    spectrum = np.fft.fft(source)
    resized = np.zeros(int(target_length), dtype=np.complex64)
    _copy_resized_spectrum(spectrum, resized, source_length)
    interpolated = np.fft.ifft(resized) * (float(target_length) / source_length)
    return interpolated.real.astype(np.float32, copy=False)


def band_limited_ifft_abs(
    spectrum,
    output_length: int,
    harmonic_count: int,
) -> np.ndarray:
    fft_values = np.asarray(spectrum, dtype=np.complex64).reshape(-1)
    if output_length <= 0:
        raise ValueError("output_length must be positive.")
    if harmonic_count < 1:
        raise ValueError("harmonic_count must be positive.")
    count = min(int(harmonic_count), int(output_length), fft_values.size)
    band_limited = fft_values[:count].copy() * 2.0
    band_limited[0] = fft_values[0]
    result = np.abs(np.fft.ifft(band_limited, n=int(output_length)))
    return result.astype(np.float32, copy=False)


def _copy_resized_spectrum(
    spectrum: np.ndarray,
    resized: np.ndarray,
    source_length: int,
) -> None:
    if source_length % 2 == 0:
        _copy_even_spectrum(spectrum, resized, source_length)
        return
    pivot = source_length // 2 + 1
    resized[:pivot] = spectrum[:pivot]
    resized[-(source_length // 2) :] = spectrum[pivot:]


def _copy_even_spectrum(
    spectrum: np.ndarray,
    resized: np.ndarray,
    source_length: int,
) -> None:
    half = source_length // 2
    resized[:half] = spectrum[:half]
    resized[-(source_length - half - 1) :] = spectrum[half + 1 :]
    resized[half] = spectrum[half] / 2.0
    resized[resized.size - half] = spectrum[half] / 2.0
