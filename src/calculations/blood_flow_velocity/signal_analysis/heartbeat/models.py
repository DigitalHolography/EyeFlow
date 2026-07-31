"""Typed results for cardiac timing analyses."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SystoleDetectionResult:
    systole_indexes: np.ndarray
    artery_signal_filtered: np.ndarray
    derivative_signal: np.ndarray
    min_peak_distance: int
    min_peak_height: np.float32


@dataclass(frozen=True)
class SpectralHeartbeatResult:
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
class HeartbeatAnalysisResult:
    systole: SystoleDetectionResult
    spectral: SpectralHeartbeatResult
