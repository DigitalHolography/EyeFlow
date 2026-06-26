"""Figure-facing waveform analysis orchestration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signal import average_cycle


@dataclass(frozen=True)
class VesselCycleAnalysis:
    artery: np.ndarray | None
    vein: np.ndarray | None


def vessel_cycle_analysis(
    artery_values: np.ndarray,
    vein_values: np.ndarray,
    beat_indexes: np.ndarray,
    samples: int,
) -> VesselCycleAnalysis:
    return VesselCycleAnalysis(
        artery=average_cycle(artery_values, beat_indexes, samples),
        vein=average_cycle(vein_values, beat_indexes, samples),
    )
