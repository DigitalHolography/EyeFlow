"""Figure adapter for the shared spectral cardiac-cycle calculation."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    SpectralCardiacCycleAnalysis,
    spectral_cardiac_cycle_analysis,
)


SpectrumData = SpectralCardiacCycleAnalysis


def spectrum_signal_analysis(
    values: np.ndarray,
    dt_seconds: float,
    systole_count: int = 0,
) -> SpectrumData:
    return spectral_cardiac_cycle_analysis(
        values,
        dt_seconds,
        systole_count,
    )
