"""One-signal spectral analysis for waveform-shape velocity figures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from calculations.math.arrays import nan_to_mean


@dataclass(frozen=True)
class SpectrumData:
    frequencies: np.ndarray
    magnitude: np.ndarray
    phase: np.ndarray
    peak_indexes: np.ndarray
    fundamental_hz: float
    heart_rate_bpm: float
    heart_rate_se_bpm: float


def spectrum_signal_analysis(values: np.ndarray, dt_seconds: float) -> SpectrumData:
    clean = nan_to_mean(values)
    window = np.hamming(clean.size)
    windowed = clean * window
    padded = np.pad(windowed, (0, clean.size * 2))
    fft = np.fft.rfft(padded)
    freq = np.fft.rfftfreq(padded.size, dt_seconds)
    mag = np.abs(fft) / max(float(np.mean(window)), np.finfo(np.float32).eps)
    mag = mag / max(float(np.nanmax(mag)), np.finfo(np.float32).eps)
    peaks, _ = signal.find_peaks(
        mag,
        distance=max(1, int(0.25 / dt_seconds)),
        prominence=max(float(np.nanmax(mag)) * 0.05, np.finfo(np.float32).eps),
    )
    if peaks.size == 0:
        peaks, _ = signal.find_peaks(
            mag,
            distance=max(1, int(0.25 / dt_seconds)),
        )
    peaks = peaks[np.argsort(freq[peaks])]
    fundamental, heart_rate, heart_rate_se = spectral_heart_rate(freq, peaks)
    return SpectrumData(
        freq.astype(np.float32),
        mag.astype(np.float32),
        np.angle(fft).astype(np.float32),
        peaks.astype(np.int32),
        fundamental,
        heart_rate,
        heart_rate_se,
    )


def spectral_heart_rate(
    frequencies: np.ndarray,
    peak_indexes: np.ndarray,
) -> tuple[float, float, float]:
    if peak_indexes.size == 0:
        return np.nan, np.nan, np.nan
    fundamental = float(frequencies[peak_indexes[0]])
    if not np.isfinite(fundamental) or fundamental <= 0:
        return fundamental, np.nan, np.nan
    peak_freqs = frequencies[peak_indexes]
    harmonic_numbers = np.rint(peak_freqs / fundamental)
    valid = harmonic_numbers > 0
    if not np.any(valid):
        return fundamental, np.nan, np.nan
    numbers = harmonic_numbers[valid]
    locations = peak_freqs[valid]
    heart_rate_hz = float(np.dot(numbers, locations) / np.dot(numbers, numbers))
    residuals = locations - numbers * heart_rate_hz
    heart_rate_se_hz = float(np.sqrt(np.mean(residuals**2)) / np.sqrt(locations.size))
    return fundamental, heart_rate_hz * 60.0, heart_rate_se_hz * 60.0
