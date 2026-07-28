"""Shared context and helpers for waveform-shape velocity PNG figures."""

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

if TYPE_CHECKING:
    from calculations.blood_flow_velocity import PerBeatAnalysisResult
    from calculations.blood_flow_velocity.signal_analysis.heartbeat import (
        SpectralHeartbeatResult,
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
    analysis: dict[str, object]
    per_beat_result: PerBeatAnalysisResult
    log: object | None = None

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
    def heartbeat(self) -> SpectralHeartbeatResult:
        """MATLAB-compatible spectral heartbeat used for all figure timing."""
        heartbeat = getattr(self.per_beat_result, "heartbeat", None)
        if heartbeat is None:
            raise ValueError("Pulse figures require the shared spectral heartbeat result.")
        return heartbeat

    @property
    def heartbeat_period_seconds(self) -> float:
        period = float(self.heartbeat.period_seconds)
        if not np.isfinite(period) or period <= 0:
            raise ValueError("Pulse figures require a positive spectral heartbeat period.")
        return period


def _output_stem(output) -> str:
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    return str(stem or "eyeflow")
def _log(ctx: PulseFigureContext, message: str) -> None:
    if callable(ctx.log):
        ctx.log(message)
def _matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=True)
    return matplotlib
def _plt():
    _matplotlib()
    import matplotlib.pyplot as plt

    return plt
