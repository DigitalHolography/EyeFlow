"""Figure adapter for the shared spectral heartbeat calculation."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.heartbeat import (
    SpectralHeartbeatResult,
    spectral_heartbeat_analysis,
)


SpectrumData = SpectralHeartbeatResult


def spectrum_signal_analysis(
    values: np.ndarray,
    dt_seconds: float,
    systole_count: int = 0,
) -> SpectrumData:
    return spectral_heartbeat_analysis(
        values,
        dt_seconds,
        systole_count,
    )
