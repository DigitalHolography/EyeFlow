"""Pulse metrics for blood-flow velocity waveforms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cycles import average_cycle


@dataclass(frozen=True)
class PulseMetricData:
    name: str
    value: float
    maximum: float
    minimum: float
    mean: float


def pulse_metric(cycle: np.ndarray, metric_name: str) -> PulseMetricData:
    values = np.asarray(cycle, dtype=np.float32).reshape(-1)
    return _pulse_metric(values, metric_name)


def pulse_metric_from_cycles(
    cycles: np.ndarray,
    metric_name: str,
    *,
    beat_axis: int = 0,
) -> PulseMetricData:
    """Calculate a pulse metric from the mean of phase-aligned beats."""
    values = np.asarray(cycles, dtype=float)
    metric_name = _metric_name(metric_name)
    if values.size == 0:
        return PulseMetricData(metric_name, np.nan, np.nan, np.nan, np.nan)
    cycle = np.mean(values, axis=beat_axis)
    if cycle.size == 0 or not np.all(np.isfinite(cycle)):
        return PulseMetricData(metric_name, np.nan, np.nan, np.nan, np.nan)
    return _pulse_metric(np.asarray(cycle).reshape(-1), metric_name)


def _pulse_metric(values: np.ndarray, metric_name: str) -> PulseMetricData:
    metric_name = _metric_name(metric_name)
    values = np.asarray(values).reshape(-1)
    if values.size == 0 or not np.any(np.isfinite(values)):
        return PulseMetricData(metric_name, np.nan, np.nan, np.nan, np.nan)

    maximum = float(np.max(values))
    minimum = float(np.min(values))
    mean = float(np.mean(values))
    excursion = maximum - minimum

    with np.errstate(divide="ignore", invalid="ignore"):
        if metric_name == "RI":
            value = float(np.divide(np.float64(excursion), np.float64(maximum)))
            value = float(np.clip(value, 0.0, 1.0))
        elif metric_name == "PI":
            value = float(np.divide(np.float64(excursion), np.float64(mean)))
            value = max(value, 0.0)
        else:
            value = float(np.divide(np.float64(maximum), np.float64(minimum)))

    return PulseMetricData(metric_name, value, maximum, minimum, mean)


def _metric_name(metric_name: str) -> str:
    metric_name = str(metric_name).upper()
    if metric_name not in {"RI", "PI", "PR"}:
        raise ValueError(f"Unsupported pulse metric '{metric_name}'.")
    return metric_name


def pulse_metric_from_signal(
    values: np.ndarray,
    peaks: np.ndarray,
    metric_name: str,
    *,
    samples: int = 60,
) -> PulseMetricData | None:
    """Calculate a pulse metric from the interpolated mean of repeating cycles."""
    cycle = average_cycle(values, peaks, samples)
    if cycle is None:
        return None
    return pulse_metric(cycle, metric_name)
