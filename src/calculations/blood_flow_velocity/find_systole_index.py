"""Systole detection from BloodFlowVelocity/find_systole_index.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SystoleDetectionResult:
    systole_indexes: np.ndarray
    artery_signal_filtered: np.ndarray
    derivative_signal: np.ndarray
    min_peak_distance: int
    min_peak_height: np.float32


def find_systole_index(
    pulse_artery,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
    min_duration_seconds: float = 0.5,
    validation_distance: int = 10,
) -> SystoleDetectionResult:
    butter, filtfilt, find_peaks = _scipy_signal_dependencies()
    pulse = np.asarray(pulse_artery, dtype=np.float32).reshape(-1)
    filtered = _lowpass_filter(pulse, dt_seconds, lowpass_freq_hz, butter, filtfilt)
    derivative = np.gradient(filtered).astype(np.float32)
    min_peak_height = np.float32(np.percentile(derivative, 95))
    min_peak_distance = _min_peak_distance(dt_seconds, min_duration_seconds)
    peaks, _ = find_peaks(
        derivative,
        height=min_peak_height,
        distance=min_peak_distance,
    )
    indexes = _validate_peaks(peaks.astype(np.int32), validation_distance)
    if indexes.size == 0:
        raise ValueError("No systole peaks detected. Check signal quality or parameters.")
    return SystoleDetectionResult(
        systole_indexes=indexes,
        artery_signal_filtered=filtered,
        derivative_signal=derivative,
        min_peak_distance=min_peak_distance,
        min_peak_height=min_peak_height,
    )


def _lowpass_filter(
    pulse: np.ndarray,
    dt_seconds: float,
    lowpass_freq_hz: float,
    butter,
    filtfilt,
) -> np.ndarray:
    normalized_cutoff = _normalized_cutoff(dt_seconds, lowpass_freq_hz)
    b, a = butter(4, normalized_cutoff, btype="low")
    return filtfilt(b, a, pulse).astype(np.float32)


def _normalized_cutoff(dt_seconds: float, lowpass_freq_hz: float) -> float:
    if dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive for systole detection.")
    nyquist_hz = 0.5 / float(dt_seconds)
    if nyquist_hz <= 0:
        raise ValueError("Nyquist frequency must be positive for systole detection.")
    return float(np.clip(float(lowpass_freq_hz) / nyquist_hz, 1e-6, 0.99))


def _min_peak_distance(dt_seconds: float, min_duration_seconds: float) -> int:
    if dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive for systole detection.")
    return max(1, int(np.floor(float(min_duration_seconds) / float(dt_seconds))))


def _validate_peaks(peaks: np.ndarray, min_distance: int) -> np.ndarray:
    if peaks.size == 0:
        return peaks.astype(np.int32, copy=False)
    validated = [int(peaks[0])]
    for peak in peaks[1:]:
        if int(peak) - validated[-1] >= int(min_distance):
            validated.append(int(peak))
    return np.asarray(validated, dtype=np.int32)


def _scipy_signal_dependencies():
    try:
        from scipy.signal import butter, filtfilt, find_peaks
    except ModuleNotFoundError as exc:
        raise ImportError("Systole detection requires scipy.") from exc
    return butter, filtfilt, find_peaks
