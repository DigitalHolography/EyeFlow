"""Paired artery/vein cycle analysis helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cycles import average_cycle


@dataclass(frozen=True)
class PairedVesselCycles:
    artery: np.ndarray | None
    vein: np.ndarray | None


def paired_vessel_cycles(
    artery_values: np.ndarray,
    vein_values: np.ndarray,
    beat_indexes: np.ndarray,
    samples: int,
) -> PairedVesselCycles:
    return PairedVesselCycles(
        artery=average_cycle(artery_values, beat_indexes, samples),
        vein=average_cycle(vein_values, beat_indexes, samples),
    )
