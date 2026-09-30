"""Cardiac-cycle detection from global retinal-velocity signals."""

from __future__ import annotations

from collections.abc import Mapping

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
    cardiac_cycles_from_available_vessel,
)

from .signal_processing import DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ


def detect_cardiac_cycles(
    estimation: Mapping[str, object],
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
) -> tuple[CardiacCycleAnalysis, str]:
    """Detect systolic boundaries from artery, vein, or a fallback record."""

    return cardiac_cycles_from_available_vessel(
        estimation["retinal_artery_velocity_signal"],
        estimation["retinal_vein_velocity_signal"],
        dt_seconds=dt_seconds,
        lowpass_freq_hz=lowpass_freq_hz,
    )


__all__ = ["detect_cardiac_cycles"]
