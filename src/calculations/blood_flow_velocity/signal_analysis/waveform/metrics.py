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

    maximum = float(np.nanmax(values))
    minimum = float(np.nanmin(values))
    mean = float(np.nanmean(values))
    excursion = maximum - minimum

    if metric_name == "RI":
        value = excursion / maximum if maximum else np.nan
        value = float(np.clip(value, 0.0, 1.0)) if np.isfinite(value) else np.nan
    elif metric_name == "PI":
        value = excursion / mean if mean else np.nan
        value = max(0.0, float(value)) if np.isfinite(value) else np.nan
    else:
        value = maximum / minimum if minimum else np.nan

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
