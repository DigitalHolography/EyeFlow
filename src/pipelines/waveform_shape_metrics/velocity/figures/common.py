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
