"""Signal preparation for the retinal-velocity core pipeline."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from calculations.math import butter_lowpass_filtfilt

from .models import RetinalVelocity, RetinalVelocityData, VesselVelocity

DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ = 15.0


def build_retinal_velocity(
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
    artery_filtered = _filtered_vessel_signal(
        artery_raw,
        "artery",
        cardiac_cycle,
        cardiac_cycle_source,
        dt_seconds,
        lowpass_freq_hz,
    )
    vein_filtered = _filtered_vessel_signal(
        vein_raw,
        "vein",
        cardiac_cycle,
        cardiac_cycle_source,
        dt_seconds,
        lowpass_freq_hz,
    )
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


def _filtered_vessel_signal(
    signal: np.ndarray,
    vessel: str,
    cardiac_cycle: CardiacCycleAnalysis,
    cardiac_cycle_source: str,
    dt_seconds: float,
    lowpass_freq_hz: float,
) -> np.ndarray:
    if cardiac_cycle_source == vessel:
        return np.asarray(cardiac_cycle.systole.signal_filtered, dtype=np.float32)
    return _filter(signal, dt_seconds, lowpass_freq_hz)


__all__ = [
    "DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ",
    "build_retinal_velocity",
]
