"""Velocity signal filtering from BloodFlowVelocity/pulseAnalysis.m."""

from __future__ import annotations

import numpy as np


def lowpass_velocity_signal(
    signal,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
) -> np.ndarray:
    butter, filtfilt = _scipy_signal_dependencies()
    values = np.asarray(signal, dtype=np.float32).reshape(-1)
    normalized_cutoff = _normalized_cutoff(dt_seconds, lowpass_freq_hz)
    b, a = butter(4, normalized_cutoff, btype="low")
    if values.size <= _filtfilt_pad_length(b, a):
        return values.copy()
    return filtfilt(b, a, values).astype(np.float32)


def _normalized_cutoff(dt_seconds: float, lowpass_freq_hz: float) -> float:
    if dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive for velocity filtering.")
    nyquist_hz = 0.5 / float(dt_seconds)
    if nyquist_hz <= 0:
        raise ValueError("Nyquist frequency must be positive for velocity filtering.")
    return float(np.clip(float(lowpass_freq_hz) / nyquist_hz, 1e-6, 0.99))


def _scipy_signal_dependencies():
    try:
        from scipy.signal import butter, filtfilt
    except ModuleNotFoundError as exc:
        raise ImportError("Velocity signal filtering requires scipy.") from exc
    return butter, filtfilt


def _filtfilt_pad_length(b, a) -> int:
    return 3 * max(len(a), len(b))
