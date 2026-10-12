"""Workflow-specific cardiac-cycle detection from physical retinal velocity."""

from __future__ import annotations

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
    cardiac_cycles_from_available_vessel,
)

from .models import RetinalVelocityData
from .signal_processing import DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ


def detect_velocity_cardiac_cycles(
    velocity_data: RetinalVelocityData,
    *,
    dt_seconds: float,
) -> tuple[CardiacCycleAnalysis, str]:
    """Detect cycles from signed, background-corrected velocity in mm/s.

    The detector owns low-pass filtering; pass the estimator's unfiltered
    whole-vessel signals to avoid filtering them twice.
    """
    return cardiac_cycles_from_available_vessel(
        velocity_data.artery.velocity,
        velocity_data.vein.velocity,
        dt_seconds=float(dt_seconds),
        lowpass_freq_hz=DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
    )


__all__ = ["detect_velocity_cardiac_cycles"]
