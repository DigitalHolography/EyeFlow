"""Fourier transform helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np


def rfft_normalized(
    signal: np.ndarray,
    axis: int = -1,
    dtype=np.float32,
) -> np.ndarray:
    """Return the real FFT normalized by the number of samples."""
    values = np.asarray(signal, dtype=dtype)
    sample_count = values.shape[axis]
    if sample_count <= 0:
        raise ValueError("The Fourier transform axis must contain samples")
    return np.fft.rfft(values, axis=axis) / float(sample_count)


def _axis_slice(ndim: int, axis: int, start: int, stop: int) -> tuple[slice, ...]:
    selection = [slice(None)] * ndim
    selection[axis] = slice(start, stop)
    return tuple(selection)


def truncate_harmonics(
    coefficients: np.ndarray,
    max_harmonic: int,
    axis: int = -1,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Keep harmonics ``0..H`` and return retained/full-truncated arrays."""
    values = np.asarray(coefficients)
    harmonic_count = values.shape[axis]
    if harmonic_count == 0:
        raise ValueError("The Fourier coefficient axis must not be empty")
    if max_harmonic < 0:
        raise ValueError("max_harmonic must be non-negative")

    highest = min(int(max_harmonic), harmonic_count - 1)
    kept = np.array(
        values[_axis_slice(values.ndim, axis, 0, highest + 1)],
        copy=True,
    )
    truncated = np.zeros_like(values)
    truncated[_axis_slice(values.ndim, axis, 0, highest + 1)] = kept
    return kept, truncated, highest


def irfft_normalized(
    coefficients: np.ndarray,
    sample_count: int,
    axis: int = -1,
) -> np.ndarray:
    """Reconstruct a signal from coefficients normalized by sample count."""
    if sample_count <= 0:
        raise ValueError("sample_count must be positive")
    values = np.asarray(coefficients)
    return np.fft.irfft(values * float(sample_count), n=sample_count, axis=axis)


def harmonic_pack(
    signal: np.ndarray,
    max_harmonic: int,
    axis: int = -1,
    dtype=np.float32,
) -> dict[str, np.ndarray | int | None]:
    """Return full/truncated normalized Fourier representations."""
    values = np.asarray(signal, dtype=dtype)
    sample_count = values.shape[axis]
    if sample_count < 2:
        return {"V": None, "H": 0, "vb": None, "Vfull": None}

    full = rfft_normalized(values, axis=axis, dtype=dtype)
    kept, truncated, highest = truncate_harmonics(full, max_harmonic, axis=axis)
    reconstructed = irfft_normalized(truncated, sample_count, axis=axis)
    return {"V": kept, "H": highest, "vb": reconstructed, "Vfull": full}


def next_power_of_two(value: int) -> int:
    if value < 1:
        raise ValueError("next_power_of_two expects a strictly positive integer.")
    return 1 << (value - 1).bit_length()


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
