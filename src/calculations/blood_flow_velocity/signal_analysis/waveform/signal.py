"""Waveform signal calculations used by blood-flow velocity analyses."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal


@dataclass(frozen=True)
class PulseMetricData:
    name: str
    value: float
    maximum: float
    minimum: float
    mean: float


@dataclass(frozen=True)
class ArterialWaveformAnalysis:
    pulse_time: np.ndarray
    period_seconds: float
    gradient: np.ndarray
    padded_time: np.ndarray
    padded_signal: np.ndarray
    padded_gradient: np.ndarray
    peak_indexes: np.ndarray
    end_min_index: int
    notch_index: int | None


@dataclass(frozen=True)
class VenousWaveformAnalysis:
    pulse_time: np.ndarray
    period_seconds: float
    padded_time: np.ndarray
    padded_signal: np.ndarray
    peak_index: int
    trough_index: int


def average_cycle(values: np.ndarray, peaks: np.ndarray, samples: int) -> np.ndarray | None:
    if peaks.size < 2:
        return None
    cycles = []
    base_x = np.linspace(0, 1, samples, dtype=np.float32)
    for start, stop in zip(peaks[:-1], peaks[1:]):
        if stop <= start + 1:
            continue
        beat = values[start:stop]
        cycles.append(np.interp(base_x, np.linspace(0, 1, beat.size), beat))
    if not cycles:
        return None
    return np.nanmean(np.asarray(cycles, dtype=np.float32), axis=0).astype(np.float32)


def cycle_extrema(values: np.ndarray, peaks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    maxima: list[int] = []
    minima: list[int] = []
    if peaks.size == 0:
        return np.asarray(maxima, dtype=np.int32), np.asarray(minima, dtype=np.int32)
    for first, second in zip(peaks[:-1], peaks[1:]):
        if second <= first + 2:
            continue
        midpoint = first + int(round((second - first) / 2))
        maxima.append(first + int(np.nanargmax(values[first:midpoint])))
        minima.append(midpoint + int(np.nanargmin(values[midpoint:second])))
    if peaks[0] > 0:
        minima.insert(0, int(np.nanargmin(values[: peaks[0]])))
    if peaks[-1] < values.size:
        maxima.append(peaks[-1] + int(np.nanargmax(values[peaks[-1]:])))
    return np.asarray(maxima, dtype=np.int32), np.asarray(minima, dtype=np.int32)


def pulse_metric(cycle: np.ndarray, metric_name: str) -> PulseMetricData:
    maximum = float(np.nanmax(cycle))
    minimum = float(np.nanmin(cycle))
    mean = float(np.nanmean(cycle))
    if metric_name == "RI" and maximum:
        value = (maximum - minimum) / maximum
    else:
        value = (maximum - minimum) / mean if mean else np.nan
    return PulseMetricData(metric_name, float(value), maximum, minimum, mean)


def arterial_waveform_analysis(
    cycle: np.ndarray,
    beat_indexes: np.ndarray,
    dt_seconds: float,
) -> ArterialWaveformAnalysis:
    period = mean_period_seconds(beat_indexes, dt_seconds, default_samples=cycle.size)
    pulse_time = cycle_time(cycle, period)
    gradient = np.gradient(cycle)
    peak_indexes, properties = signal.find_peaks(
        cycle,
        height=float(np.nanmax(cycle)) * 0.3,
        distance=max(1, cycle.size // 4),
    )
    if peak_indexes.size == 0:
        peak_indexes = np.asarray([int(np.nanargmax(cycle))])
    else:
        order = np.argsort(properties["peak_heights"])[::-1][:2]
        peak_indexes = peak_indexes[order]
    end_min = int(np.nanargmin(cycle[int(0.75 * cycle.size):]) + int(0.75 * cycle.size))
    padded_time, padded_signal = padded_cycle(pulse_time, cycle, period)
    _, padded_gradient = padded_cycle(pulse_time, gradient, period)
    return ArterialWaveformAnalysis(
        pulse_time=pulse_time,
        period_seconds=period,
        gradient=gradient,
        padded_time=padded_time,
        padded_signal=padded_signal,
        padded_gradient=padded_gradient,
        peak_indexes=peak_indexes.astype(np.int32),
        end_min_index=end_min,
        notch_index=dicrotic_notch_index(cycle, peak_indexes),
    )


def venous_waveform_analysis(
    cycle: np.ndarray,
    beat_indexes: np.ndarray,
    dt_seconds: float,
) -> VenousWaveformAnalysis:
    period = mean_period_seconds(beat_indexes, dt_seconds, default_samples=cycle.size)
    pulse_time = cycle_time(cycle, period)
    padded_time, padded_signal = padded_cycle(pulse_time, cycle, period)
    return VenousWaveformAnalysis(
        pulse_time=pulse_time,
        period_seconds=period,
        padded_time=padded_time,
        padded_signal=padded_signal,
        peak_index=int(np.nanargmax(cycle)),
        trough_index=int(np.nanargmin(cycle)),
    )


def dicrotic_notch_index(cycle: np.ndarray, peak_indexes: np.ndarray) -> int | None:
    if peak_indexes.size < 2:
        return None
    first, second = sorted(int(index) for index in peak_indexes[:2])
    if second <= first:
        return None
    notch = first + int(np.nanargmin(cycle[first : second + 1]))
    return notch if second - notch > cycle.size * 0.05 else None


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


def cycle_time(cycle: np.ndarray, period_seconds: float) -> np.ndarray:
    return np.linspace(0, period_seconds, cycle.size, dtype=np.float32)


def mean_period_seconds(
    beat_indexes: np.ndarray,
    dt_seconds: float,
    *,
    default_samples: int = 128,
) -> float:
    peaks = safe_indexes(beat_indexes)
    if peaks.size > 1:
        return float(np.nanmean(np.diff(peaks)) * dt_seconds)
    return float(dt_seconds * default_samples)


def safe_indexes(values) -> np.ndarray:
    if values is None:
        return np.asarray([], dtype=np.int32)
    indexes = np.asarray(values, dtype=np.int32).reshape(-1)
    return indexes[indexes >= 0]


def vector(values) -> np.ndarray:
    if values is None:
        return np.asarray([], dtype=np.float32)
    return np.asarray(values, dtype=np.float32).reshape(-1)


def finite_image(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(data)
    if not np.any(finite):
        return np.zeros_like(data)
    clean = data.copy()
    clean[~finite] = float(np.nanmin(data[finite]))
    return clean


def nan_to_mean(values: np.ndarray) -> np.ndarray:
    clean = np.asarray(values, dtype=np.float32).copy()
    finite = np.isfinite(clean)
    if not np.any(finite):
        return np.zeros_like(clean)
    clean[~finite] = float(np.nanmean(clean[finite]))
    return clean


def standardize(values: np.ndarray) -> np.ndarray:
    clean = nan_to_mean(values)
    std = float(np.nanstd(clean))
    if std <= 0:
        return np.zeros_like(clean)
    return ((clean - float(np.nanmean(clean))) / std).astype(np.float32)


def rescale(values: np.ndarray) -> np.ndarray:
    data = finite_image(values).astype(np.float32)
    vmin = float(np.nanmin(data))
    vmax = float(np.nanmax(data))
    if vmax <= vmin:
        return np.zeros_like(data)
    return (data - vmin) / (vmax - vmin)
