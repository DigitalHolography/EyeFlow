"""Orchestration for systolic and spectral cardiac-cycle analyses."""

from __future__ import annotations

import numpy as np

from .models import CardiacCycleAnalysis, SystoleDetectionResult
from .spectral import spectral_cardiac_cycle_analysis
from .systole import (
    DEFAULT_MIN_PERIOD_SECONDS,
    SystoleDetectionError,
    find_systole_index,
)


def analyze_cardiac_cycles(
    arterial_signal,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
) -> CardiacCycleAnalysis:
    systole = find_systole_index(
        arterial_signal,
        dt=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
    )
    spectral = spectral_cardiac_cycle_analysis(
        systole.signal_filtered,
        dt_seconds,
        systole.systole_indexes.size,
    )
    return CardiacCycleAnalysis(systole=systole, spectral=spectral)


def cardiac_cycles_from_available_vessel(
    artery_signal,
    vein_signal,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = 15.0,
) -> tuple[CardiacCycleAnalysis, str]:
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
                analyze_cardiac_cycles(
                    signal,
                    dt_seconds=dt_seconds,
                    lowpass_freq_hz=lowpass_freq_hz,
                ),
                source_name,
            )
        except SystoleDetectionError:
            continue
    return missing_vessel_cardiac_cycles(artery.size, dt_seconds), "none"


def missing_vessel_cardiac_cycles(
    signal_length: int,
    dt_seconds: float,
) -> CardiacCycleAnalysis:
    """Return one full-record cycle when no vessel can provide timing."""

    if signal_length < 2:
        raise ValueError("At least two frames are required for per-beat analysis.")
    if not np.isfinite(dt_seconds) or dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive.")
    boundaries = np.asarray([0, signal_length - 1], dtype=np.int32)
    missing = np.full(signal_length, np.nan, dtype=np.float32)
    return CardiacCycleAnalysis(
        systole=SystoleDetectionResult(
            systole_indexes=boundaries,
            signal_filtered=missing,
            derivative_signal=missing.copy(),
            min_peak_distance=max(
                1,
                int(np.floor(float(DEFAULT_MIN_PERIOD_SECONDS) / dt_seconds)),
            ),
            min_peak_height=np.float32(np.nan),
        ),
        spectral=spectral_cardiac_cycle_analysis(
            missing,
            dt_seconds,
            systole_count=0,
        ),
    )
