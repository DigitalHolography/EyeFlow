"""Shared context and helpers for waveform velocity PNG diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from calculations.math.arrays import (
    as_float32_vector as _vector,
    as_nonnegative_int_indexes as _safe_indexes,
    finite_image as _finite_image,
)
from utils.logger import Logger

from .signal_inputs import display_frequency, display_velocity

__all__ = [
    "PulseFigureContext",
    "_finite_image",
    "_log",
    "_matplotlib",
    "_output_stem",
    "_plt",
    "_safe_indexes",
    "_vector",
    "display_frequency",
    "display_velocity",
]

if TYPE_CHECKING:
    from calculations.blood_flow_velocity import PerBeatAnalysisResult
    from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
        SpectralCardiacCycleAnalysis,
    )
    from pipelines.retinal_velocity.models import RetinalVelocity


@dataclass(frozen=True)
class PulseFigureContext:
    output: object
    stem: str
    time: np.ndarray
    dt_seconds: float
    moment0_average: np.ndarray
    artery_mask: np.ndarray
    vein_mask: np.ndarray
    section_mask: np.ndarray
    retinal_velocity: RetinalVelocity
    per_beat_result: PerBeatAnalysisResult

    @property
    def artery_section_mask(self) -> np.ndarray:
        return self.artery_mask & self.section_mask

    @property
    def vein_section_mask(self) -> np.ndarray:
        return self.vein_mask & self.section_mask

    @property
    def vessel_section_mask(self) -> np.ndarray:
        return (self.artery_mask | self.vein_mask) & self.section_mask

    @property
    def cycle_boundary_indexes(self) -> np.ndarray:
        """Systoles used only to align and delimit individual cycles."""
        return _safe_indexes(self.per_beat_result.cycle_boundary_indexes)

    @property
    def cardiac_cycle(self) -> SpectralCardiacCycleAnalysis:
        """Shared spectral cardiac-cycle analysis used for figure timing."""
        cardiac_cycle = getattr(self.per_beat_result, "cardiac_cycle", None)
        if cardiac_cycle is None:
            raise ValueError(
                "Pulse figures require the shared spectral cardiac-cycle analysis."
            )
        return cardiac_cycle

    @property
    def cardiac_cycle_period_seconds(self) -> float:
        period = float(self.cardiac_cycle.period_seconds)
        if not np.isfinite(period) or period <= 0:
            raise ValueError(
                "Pulse figures require a positive spectral cardiac-cycle period."
            )
        return period


def _output_stem(output) -> str:
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    return str(stem or "eyeflow")


def _log(ctx: PulseFigureContext, message: str) -> None:
    Logger.log(message)


def _matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=True)
    return matplotlib


def _plt():
    _matplotlib()
    import matplotlib.pyplot as plt

    return plt
