"""Cardiac-cycle detection from global retinal-velocity signals."""

from __future__ import annotations

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
    cardiac_cycles_from_available_vessel,
)

from .models import RetinalVelocityData
from .signal_processing import DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ


def detect_cardiac_cycles(
    data: RetinalVelocityData,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float = DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
) -> tuple[CardiacCycleAnalysis, str]:
    """Detect systolic boundaries from artery, vein, or a fallback record."""

    return cardiac_cycles_from_available_vessel(
        data.artery.velocity,
        data.vein.velocity,
        dt_seconds=dt_seconds,
        lowpass_freq_hz=lowpass_freq_hz,
    )


__all__ = ["detect_cardiac_cycles"]
