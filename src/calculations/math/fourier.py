"""Fourier transform helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np
from scipy.signal import resample


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


def interpft_axis0(values, target_length: int, *, array_module=np):
    """Fourier-resample an array's first axis while preserving empty pixels.

    Columns containing no finite samples remain NaN. This is the vectorized
    counterpart of :func:`interpft_real` used for image and profile stacks.
    Passing CuPy as ``array_module`` keeps both computation and results on the
    device, with the same real-spectrum Nyquist and missing-pixel behavior.
    """

    xp = array_module
    source = xp.asarray(values, dtype=xp.float32)
    if source.ndim == 0:
        raise ValueError("interpft expects an array with a sample axis.")
    source_length = int(source.shape[0])
    if source_length == 0:
        raise ValueError("interpft requires a non-empty sample axis.")
    if target_length <= 0:
        raise ValueError("interpft target_length must be positive.")
    if target_length == source_length:
        return source.copy()

    trailing_shape = source.shape[1:]
    flattened = source.reshape(source_length, -1)
    active = xp.any(xp.isfinite(flattened), axis=0)
    if xp is not np:
        # Resample all columns in one batch. Gathering active columns on CUDA
        # forces synchronization and adds expensive gather/scatter operations.
        if flattened.shape[1] == 0:
            return xp.empty((int(target_length), *trailing_shape), dtype=xp.float32)
        finite_columns = xp.all(xp.isfinite(flattened), axis=0)
        # cuFFT can pair real transforms internally: a NaN spectrum may poison
        # an adjacent column during irfft. Transform finite columns alongside
        # zero-filled placeholders, then restore each missing column's NaNs.
        safe_values = xp.where(finite_columns[None, :], flattened, xp.float32(0.0))
        result = _fourier_resample_real_axis0(safe_values, int(target_length), xp)
        xp.copyto(result, xp.nan, where=~finite_columns[None, :])
        # Preserve SciPy's special infinity behavior (notably a length-one
        # result), without passing non-finite spectra into batched cuFFT.
        infinite_columns = (
            xp.any(xp.isinf(flattened), axis=0)
            & active
            & ~xp.any(xp.isnan(flattened), axis=0)
        )
        if bool(xp.any(infinite_columns)):
            result[:, infinite_columns] = xp.asarray(
                interpft_axis0(xp.asnumpy(flattened[:, infinite_columns]), int(target_length)),
                dtype=xp.float32,
            )
        return result.astype(xp.float32, copy=False).reshape(int(target_length), *trailing_shape)
    interpolated = xp.full(
        (int(target_length), flattened.shape[1]),
        xp.nan,
        dtype=xp.float32,
    )
    if bool(xp.any(active)):
        result = resample(flattened[:, active], int(target_length), axis=0)
        interpolated[:, active] = result.astype(xp.float32, copy=False)
    return interpolated.reshape(int(target_length), *trailing_shape)


def _fourier_resample_real_axis0(values, target_length: int, xp):
    """Real Fourier resampling with SciPy's even-length Nyquist handling."""

    source_length = int(values.shape[0])
    spectrum = xp.fft.rfft(values, axis=0)
    resized = xp.zeros((target_length // 2 + 1, *spectrum.shape[1:]), dtype=spectrum.dtype)
    common_length = min(source_length, target_length)
    nyquist_stop = common_length // 2 + 1
    resized[:nyquist_stop] = spectrum[:nyquist_stop]
    if common_length % 2 == 0:
        resized[common_length // 2] *= xp.float32(2.0 if target_length < source_length else 0.5)
    result = xp.fft.irfft(resized, n=target_length, axis=0)
    result *= xp.float32(target_length / source_length)
    return result


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
