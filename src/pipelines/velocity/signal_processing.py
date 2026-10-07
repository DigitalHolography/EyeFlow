"""Signal preparation for the velocity core pipeline."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from calculations.math import butter_lowpass_filtfilt

from .models import RetinalVelocity, RetinalVelocityData, VesselVelocity

DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ = 15.0


def build_velocity(
    data: RetinalVelocityData,
    cardiac_cycle: CardiacCycleAnalysis,
    cardiac_cycle_source: str,
    *,
    dt_seconds: float,
    index_base: int = 0,
    lowpass_freq_hz: float = DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
) -> RetinalVelocity:
    """Build the canonical typed result from an estimator output."""

    artery_raw = np.asarray(data.artery.velocity, dtype=np.float32)
    vein_raw = np.asarray(data.vein.velocity, dtype=np.float32)
    # Workflow cycles are detected from frequency; filter physical velocity here.
    artery_filtered = _filter(artery_raw, dt_seconds, lowpass_freq_hz)
    vein_filtered = _filter(vein_raw, dt_seconds, lowpass_freq_hz)
    return RetinalVelocity(
        maps=data.maps,
        artery=VesselVelocity(
            signals=data.artery,
            velocity_filtered=artery_filtered,
        ),
        vein=VesselVelocity(
            signals=data.vein,
            velocity_filtered=vein_filtered,
        ),
        vessel_frms_background=data.vessel_frms_background,
        cardiac_cycle=cardiac_cycle,
        cardiac_cycle_source=cardiac_cycle_source,
        dt_seconds=float(dt_seconds),
        provenance=data.provenance,
        index_base=int(index_base),
    )


def _filter(signal, dt_seconds: float, lowpass_freq_hz: float) -> np.ndarray:
    return butter_lowpass_filtfilt(
        signal,
        dt_seconds=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
        order=4,
    )


__all__ = [
    "DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ",
    "build_velocity",
]
