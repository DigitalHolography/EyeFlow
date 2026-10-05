"""Filtering operations shared by waveform-velocity products."""

from __future__ import annotations

import numpy as np

from calculations.math import butter_lowpass_filtfilt
from pipelines.retinal_velocity.signal_processing import (
    DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
)


def lowpass_velocity_signals(
    values,
    *,
    dt_seconds: float,
    axis: int = -1,
) -> np.ndarray:
    """Low-pass each velocity signal independently along ``axis``."""

    signals = np.asarray(values, dtype=np.float32)
    if signals.ndim == 0:
        raise ValueError("Velocity signals must include a sample axis.")
    normalized_axis = int(axis)
    if normalized_axis < 0:
        normalized_axis += signals.ndim
    if not 0 <= normalized_axis < signals.ndim:
        raise ValueError(f"axis {axis} is out of bounds for {signals.ndim} dimensions.")
    moved = np.moveaxis(signals, normalized_axis, -1)
    flattened = moved.reshape(-1, moved.shape[-1])
    filtered = np.full(flattened.shape, np.nan, dtype=np.float32)
    for index, signal in enumerate(flattened):
        if np.any(np.isfinite(signal)):
            filtered[index] = butter_lowpass_filtfilt(
                signal,
                dt_seconds=np.float32(dt_seconds),
                lowpass_freq_hz=np.float32(DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ),
                order=4,
            )
    restored = filtered.reshape(moved.shape)
    return np.moveaxis(restored, -1, normalized_axis)


__all__ = ["lowpass_velocity_signals"]
