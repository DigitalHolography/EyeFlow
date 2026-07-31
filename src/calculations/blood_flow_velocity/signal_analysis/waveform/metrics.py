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
    metric_name = str(metric_name).upper()
    if metric_name not in {"RI", "PI", "PR"}:
        raise ValueError(f"Unsupported pulse metric '{metric_name}'.")
    if values.size == 0 or not np.any(np.isfinite(values)):
        return PulseMetricData(metric_name, np.nan, np.nan, np.nan, np.nan)

    maximum = float(np.max(values))
    minimum = float(np.min(values))
    mean = float(np.mean(values))
    excursion = maximum - minimum

    with np.errstate(divide="ignore", invalid="ignore"):
        if metric_name == "RI":
            value = float(np.divide(np.float64(excursion), np.float64(maximum)))
            if value < 0:
                value = 0.0
            if value > 1:
                value = 1.0
        elif metric_name == "PI":
            value = float(np.divide(np.float64(excursion), np.float64(mean)))
            if value < 0:
                value = 0.0
        else:
            value = float(np.divide(np.float64(maximum), np.float64(minimum)))

    return PulseMetricData(metric_name, value, maximum, minimum, mean)


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
