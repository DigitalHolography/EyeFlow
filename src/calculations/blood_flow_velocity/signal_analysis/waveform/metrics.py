"""Pulse metrics for blood-flow velocity waveforms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PulseMetricData:
    name: str
    value: float
    maximum: float
    minimum: float
    mean: float


def pulse_metric(cycle: np.ndarray, metric_name: str) -> PulseMetricData:
    maximum = float(np.nanmax(cycle))
    minimum = float(np.nanmin(cycle))
    mean = float(np.nanmean(cycle))
    if metric_name == "RI" and maximum:
        value = (maximum - minimum) / maximum
    else:
        value = (maximum - minimum) / mean if mean else np.nan
    return PulseMetricData(metric_name, float(value), maximum, minimum, mean)
