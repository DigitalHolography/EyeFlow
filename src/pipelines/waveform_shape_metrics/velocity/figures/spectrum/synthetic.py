"""Synthetic repeated-cycle spectrum calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from calculations.blood_flow_velocity.signal_analysis.waveform.cycles import average_cycle


@dataclass(frozen=True)
class SyntheticSpectrumData:
    frequencies: np.ndarray
    magnitude: np.ndarray
    phase: np.ndarray
    peak_indexes: np.ndarray


def synthetic_spectrum_analysis(
    cycle: np.ndarray,
    dt_seconds: float,
    beat_indexes: np.ndarray,
) -> SyntheticSpectrumData:
    repeated = np.tile(cycle, 512)
    if beat_indexes.size > 1:
        period = float(np.nanmean(np.diff(beat_indexes)) * dt_seconds)
    else:
        period = cycle.size * dt_seconds
    fs = cycle.size / max(period, np.finfo(np.float32).eps)
    frequencies = np.fft.rfftfreq(repeated.size, 1.0 / fs)
    fft = np.fft.rfft(repeated)
    magnitude = np.abs(fft)
    magnitude = (
        magnitude
        / max(float(magnitude[0]), np.finfo(np.float32).eps)
        * float(np.nanmean(repeated))
    )
    phase = np.angle(fft)
    peaks, _ = signal.find_peaks(magnitude, height=0.001)
    return SyntheticSpectrumData(
        frequencies.astype(np.float32),
        magnitude.astype(np.float32),
        phase.astype(np.float32),
        peaks.astype(np.int32),
    )


def synthetic_spectrum_from_signals(
    first_values: np.ndarray,
    fallback_values: np.ndarray,
    beat_indexes: np.ndarray,
    dt_seconds: float,
    *,
    samples: int = 128,
) -> SyntheticSpectrumData | None:
    cycle = average_cycle(first_values, beat_indexes, samples)
    if cycle is None:
        cycle = average_cycle(fallback_values, beat_indexes, samples)
    if cycle is None:
        return None
    return synthetic_spectrum_analysis(cycle, dt_seconds, beat_indexes)
