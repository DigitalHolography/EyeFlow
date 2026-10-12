"""Arterial and venous morphology analysis for average velocity cycles."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from .cycles import cycle_time, padded_cycle


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


def arterial_waveform_analysis(
    cycle: np.ndarray,
    period_seconds: float,
) -> ArterialWaveformAnalysis | None:
    """Return no morphology when required extrema have no finite support."""
    cycle = np.asarray(cycle)
    finite = np.isfinite(cycle)
    end_start = int(0.75 * cycle.size)
    end_finite = np.flatnonzero(finite[end_start:])
    if cycle.size < 2 or not np.any(finite) or not end_finite.size:
        return None
    extrema_cycle = np.where(finite, cycle, np.nan)
    period = float(period_seconds)
    pulse_time = cycle_time(cycle, period)
    gradient = np.gradient(cycle)
    peak_indexes, properties = signal.find_peaks(
        extrema_cycle,
        height=float(np.max(cycle[finite])) * 0.3,
        distance=max(1, cycle.size // 4),
    )
    if peak_indexes.size == 0:
        peak_indexes = np.asarray([int(np.nanargmax(extrema_cycle))])
    else:
        order = np.argsort(properties["peak_heights"])[::-1][:2]
        peak_indexes = peak_indexes[order]
    end_min = end_start + int(end_finite[np.argmin(cycle[end_start:][end_finite])])
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
    period_seconds: float,
) -> VenousWaveformAnalysis | None:
    cycle = np.asarray(cycle)
    finite = np.isfinite(cycle)
    if not np.any(finite):
        return None
    extrema_cycle = np.where(finite, cycle, np.nan)
    period = float(period_seconds)
    pulse_time = cycle_time(cycle, period)
    padded_time, padded_signal = padded_cycle(pulse_time, cycle, period)
    return VenousWaveformAnalysis(
        pulse_time=pulse_time,
        period_seconds=period,
        padded_time=padded_time,
        padded_signal=padded_signal,
        peak_index=int(np.nanargmax(extrema_cycle)),
        trough_index=int(np.nanargmin(extrema_cycle)),
    )


def dicrotic_notch_index(cycle: np.ndarray, peak_indexes: np.ndarray) -> int | None:
    if peak_indexes.size < 2:
        return None
    first, second = sorted(int(index) for index in peak_indexes[:2])
    if second <= first:
        return None
    section = cycle[first : second + 1]
    finite = np.flatnonzero(np.isfinite(section))
    if not finite.size:
        return None
    notch = first + int(finite[np.argmin(section[finite])])
    return notch if second - notch > cycle.size * 0.05 else None
