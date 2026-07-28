"""Orchestration for systolic and spectral heartbeat analyses."""

from __future__ import annotations

import numpy as np

from .models import HeartbeatAnalysisResult
from .spectral import MATLAB_MAX_HARMONICS, spectral_heartbeat_analysis
from .systole import find_systole_index


def run_heartbeat_analysis(
    arterial_signal,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
    max_harmonics: int = MATLAB_MAX_HARMONICS,
) -> HeartbeatAnalysisResult:
    systole = find_systole_index(
        arterial_signal,
        dt=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
    )
    spectral = spectral_heartbeat_analysis(
        arterial_signal,
        dt_seconds,
        systole.systole_indexes.size,
        max_harmonics=max_harmonics,
    )
    return HeartbeatAnalysisResult(systole=systole, spectral=spectral)
