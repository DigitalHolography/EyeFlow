"""Periodic Fourier bandlimiting with an explicit non-DC harmonic count.

``harmonics=12`` always means DC plus harmonics 1 through 12. This is
independent of beat duration and is not a cutoff in Hz.
"""

from __future__ import annotations

import operator
from typing import Literal

import numpy as np

DEFAULT_BAND_LIMITED_HARMONICS = 12
Reconstruction = Literal["matlab_magnitude", "signed"]


def validate_harmonics(harmonics: int) -> int:
    """Validate the calculation API's count, accepting NumPy integer scalars."""
    if isinstance(harmonics, (bool, np.bool_)):
        # Keep the numerical API's ValueError contract for every invalid count.
        raise ValueError("harmonics must be a non-negative integer, excluding DC.")  # noqa: TRY004
    try:
        count = operator.index(harmonics)
    except TypeError as exc:
        raise ValueError("harmonics must be a non-negative integer, excluding DC.") from exc
    if count < 0:
        raise ValueError("harmonics must be a non-negative integer, excluding DC.")
    return count


def band_limit(
    signal,
    *,
    harmonics: int = DEFAULT_BAND_LIMITED_HARMONICS,
    reconstruction: Reconstruction = "matlab_magnitude",
    axis: int = -1,
) -> np.ndarray:
    """Retain DC and the first ``harmonics`` of each periodic signal.

    ``matlab_magnitude`` matches EyeFlowMatlab: double positive-frequency
    coefficients except DC, zero-pad, and take ``abs(ifft(...))``. Magnitude
    can change the mean and sign; it is not a signed Fourier-series sum.

    ``signed`` restores conjugate frequency pairs, preserving DC and signs.
    Signed reconstruction retains Nyquist once; MATLAB magnitude doubles all
    retained non-DC coefficients, including Nyquist. Counts above Nyquist clamp.
    Nonfinite samples propagate within their own series; callers that need
    incomplete series unchanged must select fully finite series first.
    """
    values = np.asarray(signal, dtype=np.float32)
    if values.ndim == 0 or values.shape[axis] == 0:
        raise ValueError("band_limit requires a nonempty sample axis.")
    return reconstruct_band_limited(
        np.fft.fft(values, axis=axis),
        harmonics=harmonics,
        reconstruction=reconstruction,
        axis=axis,
    )


def reconstruct_band_limited(
    spectrum,
    *,
    harmonics: int = DEFAULT_BAND_LIMITED_HARMONICS,
    reconstruction: Reconstruction = "matlab_magnitude",
    axis: int = -1,
) -> np.ndarray:
    """Reconstruct from a full, unnormalized FFT using :func:`band_limit` semantics."""
    harmonics = validate_harmonics(harmonics)
    if reconstruction not in ("matlab_magnitude", "signed"):
        raise ValueError("reconstruction must be 'matlab_magnitude' or 'signed'.")
    values = np.moveaxis(np.asarray(spectrum), axis, -1)
    if values.ndim == 0 or values.shape[-1] == 0:
        raise ValueError("spectrum must have a nonempty sample axis.")
    count = values.shape[-1]
    highest = min(harmonics, count // 2)
    retained = np.zeros_like(values)
    retained[..., : highest + 1] = values[..., : highest + 1]
    if reconstruction == "signed":
        if highest:
            retained[..., -highest:] = values[..., -highest:]
        result = np.fft.ifft(retained, axis=-1).real
    else:
        retained[..., 1 : highest + 1] *= 2
        result = np.abs(np.fft.ifft(retained, axis=-1))
    return np.moveaxis(result.astype(np.float32), -1, axis)
