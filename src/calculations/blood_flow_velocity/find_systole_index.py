"""Systole detection based on low-pass filtering and peak detection of the derivative of the signal."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.math import butter_lowpass_filtfilt


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
    dt: np.float32,
    lowpass_freq_hz: np.float32 = np.float32(15.0),
    min_duration_seconds: np.float32 = np.float32(0.5),
    validation_distance: int = 10,
) -> SystoleDetectionResult:
    find_peaks = _scipy_signal_dependencies()

    pulse = np.asarray(pulse_artery, dtype=np.float32).reshape(-1)
    filtered_pulse = butter_lowpass_filtfilt(
        pulse,
        dt_seconds=np.float32(dt),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
        order=4,
    )
    derivative = np.gradient(filtered_pulse).astype(np.float32)
    min_peak_height = np.float32(np.percentile(derivative, 95))
    min_peak_distance = _min_peak_distance(dt, min_duration_seconds)

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
        artery_signal_filtered=filtered_pulse,
        derivative_signal=derivative,
        min_peak_distance=min_peak_distance,
        min_peak_height=min_peak_height,
    )

def _min_peak_distance(dt_seconds: np.float32, min_duration_seconds: np.float32) -> int:
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
        from scipy.signal import find_peaks
    except ModuleNotFoundError as exc:
        raise ImportError("Systole detection requires scipy.") from exc
    return find_peaks
