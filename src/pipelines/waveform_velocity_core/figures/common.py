"""Shared context and helpers for waveform velocity PNG diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .signal_inputs import (
    array_or_none as _array_or_none,
    display_frequency,
    display_velocity,
    section_mask as _section_mask,
)
from calculations.math.arrays import (
    as_float32_vector as _vector,
    as_nonnegative_int_indexes as _safe_indexes,
    finite_image as _finite_image,
)
from utils.logger import Logger
from pipelines.retinal_velocity.semantics import (
    VelocitySemantics,
    resolve_velocity_semantics,
)

__all__ = [
    "PulseFigureContext",
    "_array_or_none",
    "_finite_image",
    "_log",
    "_matplotlib",
    "_output_stem",
    "_plt",
    "_safe_indexes",
    "_section_mask",
    "_vector",
    "display_frequency",
    "display_velocity",
]

if TYPE_CHECKING:
    from calculations.blood_flow_velocity import PerBeatAnalysisResult
    from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
        SpectralCardiacCycleAnalysis,
    )


@dataclass(frozen=True)
class PulseFigureContext:
    output: object
    stem: str
    time: np.ndarray
    dt_seconds: float
    moment0_avg: np.ndarray
    artery_mask: np.ndarray
    vein_mask: np.ndarray
    section_mask: np.ndarray
    velocity_analysis: dict[str, object]
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

    @property
    def velocity_semantics(self) -> VelocitySemantics:
        return resolve_velocity_semantics(self.velocity_analysis)

    @property
    def velocity_axis_label(self) -> str:
        return self.velocity_semantics.axis_label

    @property
    def velocity_colorbar_label(self) -> str:
        return self.velocity_semantics.unit

    @property
    def rms_quantity_label(self) -> str:
        return "RMS frequency (kHz)"

    @property
    def delta_rms_quantity_label(self) -> str:
        return "Delta Doppler RMS frequency (kHz)"

    def display_rms(self, values) -> np.ndarray:
        return display_frequency(values)


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
