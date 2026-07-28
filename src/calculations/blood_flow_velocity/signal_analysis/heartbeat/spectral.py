"""Spectral heartbeat analysis translated from SpectralWaveformAnalysis.m."""

from __future__ import annotations

import numpy as np
from scipy import signal as scipy_signal

from calculations.math.arrays import nan_to_mean

from .models import SpectralHeartbeatResult


MATLAB_PADDING_FACTOR = 16
MATLAB_MAX_HARMONICS = 6


def spectral_heartbeat_analysis(
    values,
    dt_seconds: float,
    systole_count: int,
    *,
    max_harmonics: int = MATLAB_MAX_HARMONICS,
) -> SpectralHeartbeatResult:
    """Reproduce the calculations in MATLAB SpectralWaveformAnalysis.m."""
    clean = nan_to_mean(values).reshape(-1)
    _validate_inputs(clean, dt_seconds, systole_count, max_harmonics)

    sampling_frequency_hz = 1.0 / float(dt_seconds)
    duration_seconds = clean.size * float(dt_seconds)
    estimated_fundamental_hz = float(systole_count) / duration_seconds

    window = np.hamming(clean.size)
    windowed = clean * window
    pad_width = clean.size * MATLAB_PADDING_FACTOR
    padded = np.pad(windowed, (pad_width, pad_width))
    fft_coefficients = np.fft.fft(padded)
    frequencies_hz = np.arange(
        padded.size // 2 + 1,
        dtype=np.float64,
    ) * (sampling_frequency_hz / padded.size)

    coherent_gain = max(float(np.mean(window)), np.finfo(np.float32).eps)
    magnitude = np.abs(fft_coefficients[: frequencies_hz.size]) / coherent_gain
    maximum = float(np.max(magnitude))
    if maximum > 0:
        magnitude = magnitude / maximum
    else:
        magnitude = np.zeros_like(magnitude)

    peak_indexes = _prominent_peak_indexes(
        magnitude,
        estimated_fundamental_hz,
    )
    fundamental_hz, valid_harmonics_hz = _valid_harmonics(
        magnitude,
        frequencies_hz,
        peak_indexes,
        sampling_frequency_hz,
        max_harmonics,
    )
    heart_rate_hz, heart_rate_ste_hz = _harmonic_heart_rate(
        frequencies_hz,
        peak_indexes,
        fundamental_hz,
    )
    period_seconds = (
        1.0 / heart_rate_hz
        if np.isfinite(heart_rate_hz) and heart_rate_hz > 0
        else np.nan
    )

    positive_fft = fft_coefficients[: frequencies_hz.size]
    return SpectralHeartbeatResult(
        fft_coefficients=fft_coefficients.astype(np.complex64),
        frequencies_hz=frequencies_hz.astype(np.float32),
        magnitude=magnitude.astype(np.float32),
        phase_rad=np.angle(positive_fft).astype(np.float32),
        peak_indexes=peak_indexes.astype(np.int32, copy=False),
        fundamental_hz=fundamental_hz,
        valid_harmonics_hz=valid_harmonics_hz.astype(np.float32, copy=False),
        heart_rate_hz=heart_rate_hz,
        heart_rate_bpm=heart_rate_hz * 60.0,
        heart_rate_ste_hz=heart_rate_ste_hz,
        heart_rate_ste_bpm=heart_rate_ste_hz * 60.0,
        period_seconds=period_seconds,
        estimated_fundamental_hz=estimated_fundamental_hz,
    )


def _validate_inputs(
    values: np.ndarray,
    dt_seconds: float,
    systole_count: int,
    max_harmonics: int,
) -> None:
    if values.size < 2:
        raise ValueError("signal must contain at least two samples.")
    if not np.isfinite(dt_seconds) or dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive.")
    if int(systole_count) < 0:
        raise ValueError("systole_count must be nonnegative.")
    if int(max_harmonics) < 1:
        raise ValueError("max_harmonics must be positive.")


def _prominent_peak_indexes(
    magnitude: np.ndarray,
    estimated_fundamental_hz: float,
) -> np.ndarray:
    first_distance = _matlab_peak_distance(estimated_fundamental_hz * 0.6)
    candidates, _ = scipy_signal.find_peaks(
        magnitude,
        distance=first_distance,
    )
    if candidates.size == 0:
        return np.asarray([], dtype=np.int32)

    strongest = int(candidates[np.argmax(magnitude[candidates])])
    minimum_prominence = float(magnitude[strongest]) * 0.2
    peaks, _ = scipy_signal.find_peaks(
        magnitude,
        distance=_matlab_peak_distance(estimated_fundamental_hz * 0.8),
        prominence=minimum_prominence,
    )
    return np.sort(peaks).astype(np.int32, copy=False)


def _matlab_peak_distance(value: float) -> int:
    # MATLAB receives findpeaks(Y, ...) without an X vector, so this value is
    # interpreted in spectrum-index units even though it is derived in Hz.
    if not np.isfinite(value) or value <= 1:
        return 1
    return int(np.ceil(value))


def _valid_harmonics(
    magnitude: np.ndarray,
    frequencies_hz: np.ndarray,
    peak_indexes: np.ndarray,
    sampling_frequency_hz: float,
    max_harmonics: int,
) -> tuple[float, np.ndarray]:
    if peak_indexes.size == 0:
        return np.nan, np.asarray([], dtype=np.float32)

    fundamental_hz = float(frequencies_hz[peak_indexes[0]])
    harmonics = [fundamental_hz]
    if peak_indexes.size >= 2:
        local_maxima, _ = scipy_signal.find_peaks(magnitude)
        for harmonic_number in range(2, int(max_harmonics) + 1):
            expected_hz = fundamental_hz * harmonic_number
            half_window_hz = fundamental_hz * 0.48
            in_window = local_maxima[
                (frequencies_hz[local_maxima] >= expected_hz - half_window_hz)
                & (frequencies_hz[local_maxima] <= expected_hz + half_window_hz)
            ]
            if in_window.size == 0:
                continue
            closest = int(
                in_window[
                    np.argmin(np.abs(frequencies_hz[in_window] - expected_hz))
                ]
            )
            harmonics.append(float(frequencies_hz[closest]))

    valid = np.asarray(harmonics, dtype=np.float32)
    return fundamental_hz, valid[valid <= sampling_frequency_hz / 2.0]


def _harmonic_heart_rate(
    frequencies_hz: np.ndarray,
    peak_indexes: np.ndarray,
    fundamental_hz: float,
) -> tuple[float, float]:
    if (
        peak_indexes.size == 0
        or not np.isfinite(fundamental_hz)
        or fundamental_hz <= 0
    ):
        return np.nan, np.nan

    locations = frequencies_hz[peak_indexes]
    harmonic_numbers = np.floor(locations / fundamental_hz + 0.5)
    denominator = float(np.dot(harmonic_numbers, harmonic_numbers))
    if denominator <= 0:
        return np.nan, np.nan
    heart_rate_hz = float(np.dot(harmonic_numbers, locations) / denominator)
    residuals = locations - harmonic_numbers * heart_rate_hz
    rmse_hz = float(np.sqrt(np.mean(residuals**2)))
    standard_error_hz = rmse_hz / np.sqrt(locations.size)
    return heart_rate_hz, standard_error_hz
