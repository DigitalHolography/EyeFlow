"""Filtering helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np


def butter_lowpass_filtfilt(
    signal,
    *,
    dt_seconds: np.float32,
    lowpass_freq_hz: np.float32,
    order: int = 4,
) -> np.ndarray:
    butter, filtfilt = _scipy_signal_dependencies()
    values = np.asarray(signal, dtype=np.float32).reshape(-1)
    cutoff = normalized_lowpass_cutoff(dt_seconds, lowpass_freq_hz)
    b, a = butter(int(order), cutoff, btype="low")
    if values.size <= filtfilt_pad_length(b, a):
        return values.copy()
    return filtfilt(b, a, values).astype(np.float32)


def normalized_lowpass_cutoff(dt_seconds: np.float32, lowpass_freq_hz: np.float32) -> np.float32:
    if dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive for low-pass filtering.")
    nyquist_hz = np.float32(0.5) / dt_seconds
    if nyquist_hz <= 0:
        raise ValueError("Nyquist frequency must be positive for low-pass filtering.")
    return np.float32(np.clip(float(lowpass_freq_hz) / float(nyquist_hz), 1e-6, 0.99))


def filtfilt_pad_length(b, a) -> int:
    return 3 * max(len(a), len(b))


def _scipy_signal_dependencies():
    try:
        from scipy.signal import butter, filtfilt
    except ModuleNotFoundError as exc:
        raise ImportError("Low-pass filtering requires scipy.") from exc
    return butter, filtfilt
