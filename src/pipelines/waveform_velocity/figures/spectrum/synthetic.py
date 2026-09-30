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
    period_seconds: float,
) -> SyntheticSpectrumData:
    repeated = np.tile(cycle, 512)
    period = float(period_seconds)
    fs = cycle.size / max(period, np.finfo(np.float32).eps)
    positive_count = repeated.size // 2
    frequencies = (
        np.arange(positive_count, dtype=np.float64) * (fs / repeated.size)
    )
    fft = np.fft.fft(repeated)[:positive_count]
    fft[1 : positive_count // 2] *= 2.0
    magnitude = np.abs(fft)
    magnitude = (
        magnitude
        / max(float(magnitude[0]), np.finfo(np.float32).eps)
        * float(np.mean(repeated))
    )
    phase = np.angle(fft)
    peaks, _ = signal.find_peaks(magnitude, height=0.001)
    peaks = np.concatenate((np.asarray([0], dtype=np.int64), peaks))
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
    period_seconds: float,
    *,
    samples: int = 128,
) -> SyntheticSpectrumData | None:
    cycle = average_cycle(first_values, beat_indexes, samples)
    if cycle is None:
        cycle = average_cycle(fallback_values, beat_indexes, samples)
    if cycle is None:
        return None
    return synthetic_spectrum_analysis(cycle, period_seconds)
