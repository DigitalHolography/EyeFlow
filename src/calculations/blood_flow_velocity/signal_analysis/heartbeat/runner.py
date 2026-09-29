"""Orchestration for systolic and spectral heartbeat analyses."""

from __future__ import annotations

import numpy as np

from .models import HeartbeatAnalysisResult, SystoleDetectionResult
from .spectral import spectral_heartbeat_analysis
from .systole import find_systole_index


def run_heartbeat_analysis(
    arterial_signal,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
) -> HeartbeatAnalysisResult:
    systole = find_systole_index(
        arterial_signal,
        dt=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
    )
    spectral = spectral_heartbeat_analysis(
        systole.artery_signal_filtered,
        dt_seconds,
        systole.systole_indexes.size,
    )
    return HeartbeatAnalysisResult(systole=systole, spectral=spectral)


def heartbeat_from_available_vessel(
    artery_signal,
    vein_signal,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
) -> tuple[HeartbeatAnalysisResult, str]:
    """Detect shared beat boundaries from artery, then vein, then the record."""

    artery = np.asarray(artery_signal, dtype=np.float32).reshape(-1)
    vein = np.asarray(vein_signal, dtype=np.float32).reshape(-1)
    if artery.shape != vein.shape:
        raise ValueError("artery and vein signals must have the same length.")
    for source_name, signal in (("artery", artery), ("vein", vein)):
        if not np.any(np.isfinite(signal)):
            continue
        try:
            return (
                run_heartbeat_analysis(
                    signal,
                    dt_seconds=dt_seconds,
                    lowpass_freq_hz=lowpass_freq_hz,
                ),
                source_name,
            )
        except ValueError as exc:
            if "No systole peaks detected" not in str(exc):
                raise
    return missing_vessel_heartbeat(artery.size, dt_seconds), "none"


def missing_vessel_heartbeat(
    signal_length: int,
    dt_seconds: float,
) -> HeartbeatAnalysisResult:
    """Return one full-record cycle when no vessel can provide timing."""

    if signal_length < 2:
        raise ValueError("At least two frames are required for per-beat analysis.")
    if not np.isfinite(dt_seconds) or dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive.")
    boundaries = np.asarray([0, signal_length - 1], dtype=np.int32)
    missing = np.full(signal_length, np.nan, dtype=np.float32)
    return HeartbeatAnalysisResult(
        systole=SystoleDetectionResult(
            systole_indexes=boundaries,
            artery_signal_filtered=missing,
            derivative_signal=missing.copy(),
            min_peak_distance=max(1, int(np.floor(0.5 / dt_seconds))),
            min_peak_height=np.float32(np.nan),
        ),
        spectral=spectral_heartbeat_analysis(
            missing,
            dt_seconds,
            systole_count=0,
        ),
    )
