"""Cycle extraction helpers for blood-flow velocity waveforms."""

from __future__ import annotations

import numpy as np

from calculations.math.arrays import as_nonnegative_int_indexes


def average_cycle(values: np.ndarray, peaks: np.ndarray, samples: int) -> np.ndarray | None:
    """Return a fixed-length mean cycle, matching MATLAB ``interpSignal``.

    Each interval between consecutive systolic boundaries is interpolated onto
    the same normalized time grid before the cycles are averaged.  Taking the
    extrema after this function is therefore different from taking extrema on
    the concatenated recording (or averaging per-cycle extrema).
    """
    values = np.asarray(values, dtype=np.float32).reshape(-1)
    peaks = as_nonnegative_int_indexes(peaks)
    samples = int(samples)
    if samples < 1:
        raise ValueError("samples must be positive.")
    if values.size == 0:
        return None

    if peaks.size < 2:
        return np.interp(
            np.linspace(0, values.size - 1, samples, dtype=np.float32),
            np.arange(values.size, dtype=np.float32),
            values,
        ).astype(np.float32)

    cycles: list[np.ndarray] = []
    for start, stop in zip(peaks[:-1], peaks[1:]):
        start = int(start)
        stop = int(stop)
        if start < 0 or stop > values.size or stop <= start + 1:
            continue
        source_x = np.arange(start, stop, dtype=np.float32)
        target_x = np.linspace(start, stop - 1, samples, dtype=np.float32)
        cycles.append(np.interp(target_x, source_x, values[start:stop]))
    if not cycles:
        return None

    return np.mean(np.asarray(cycles, dtype=np.float32), axis=0).astype(
        np.float32,
        copy=False,
    )


def cycle_extrema(values: np.ndarray, peaks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    maxima: list[int] = []
    minima: list[int] = []
    if peaks.size == 0:
        return np.asarray(maxima, dtype=np.int32), np.asarray(minima, dtype=np.int32)
    for first, second in zip(peaks[:-1], peaks[1:]):
        if second <= first + 2:
            continue
        midpoint = first + int(np.floor((second - first) / 2.0 + 0.5))
        maxima.append(first + int(np.nanargmax(values[first : midpoint + 1])))
        minima.append(midpoint + int(np.nanargmin(values[midpoint:second])))
    if peaks[0] > 0:
        minima.insert(0, int(np.nanargmin(values[: peaks[0] + 1])))
    if peaks[-1] < values.size:
        maxima.append(peaks[-1] + int(np.nanargmax(values[peaks[-1]:])))
    return np.asarray(maxima, dtype=np.int32), np.asarray(minima, dtype=np.int32)


def cycle_time(cycle: np.ndarray, period_seconds: float) -> np.ndarray:
    return np.linspace(0, period_seconds, cycle.size, dtype=np.float32)


def mean_period_seconds(
    beat_indexes: np.ndarray,
    dt_seconds: float,
    *,
    default_samples: int = 128,
) -> float:
    peaks = as_nonnegative_int_indexes(beat_indexes)
    if peaks.size > 1:
        return float(np.nanmean(np.diff(peaks)) * dt_seconds)
    return float(dt_seconds * default_samples)


def padded_cycle(
    pulse_time: np.ndarray,
    values: np.ndarray,
    period_seconds: float,
) -> tuple[np.ndarray, np.ndarray]:
    half = values.size // 2
    if half == 0:
        return pulse_time, values
    left_time = np.linspace(-period_seconds / 2.0, 0.0, half, endpoint=False)
    right_count = values.size - half
    right_time = np.linspace(
        period_seconds,
        period_seconds * 1.5,
        right_count,
        endpoint=False,
    )
    return (
        np.concatenate((left_time, pulse_time, right_time)).astype(np.float32),
        np.concatenate((values[-half:], values, values[:right_count])).astype(np.float32),
    )
