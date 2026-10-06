"""Typed results for cardiac-cycle timing analyses."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class SuspectedMissedBeatGap:
    """A retained systole interval close to a multiple of the typical period."""

    start_index: int
    stop_index: int
    interval_samples: int
    estimated_period_samples: float
    estimated_multiple: int


@dataclass(frozen=True)
class SystoleDetectionResult:
    systole_indexes: np.ndarray
    signal_filtered: np.ndarray
    derivative_signal: np.ndarray
    min_peak_distance: int
    min_peak_height: np.float32
    min_peak_prominence: np.float32 = field(default_factory=lambda: np.float32(np.nan))
    estimated_period_samples: np.float32 = field(default_factory=lambda: np.float32(np.nan))
    suspected_missed_beat_gaps: tuple[SuspectedMissedBeatGap, ...] = ()


@dataclass(frozen=True)
class SpectralCardiacCycleAnalysis:
    fft_coefficients: np.ndarray
    frequencies_hz: np.ndarray
    magnitude: np.ndarray
    phase_rad: np.ndarray
    peak_indexes: np.ndarray
    fundamental_hz: float
    valid_harmonics_hz: np.ndarray
    heart_rate_hz: float
    heart_rate_bpm: float
    heart_rate_ste_hz: float
    heart_rate_ste_bpm: float
    period_seconds: float
    estimated_fundamental_hz: float

    @property
    def frequencies(self) -> np.ndarray:
        """Figure-facing alias retained for existing spectrum consumers."""
        return self.frequencies_hz

    @property
    def phase(self) -> np.ndarray:
        """Figure-facing alias retained for existing spectrum consumers."""
        return self.phase_rad


@dataclass(frozen=True)
class CardiacCycleAnalysis:
    systole: SystoleDetectionResult
    spectral: SpectralCardiacCycleAnalysis
