"""Typed inputs for absolute waveform metrics."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class VesselAbsoluteWaveformInputs:
    """Raw and band-limited global/segment waveforms for one vessel."""

    raw_global: np.ndarray | None = None
    bandlimited_global: np.ndarray | None = None
    raw_segments: np.ndarray | None = None
    bandlimited_segments: np.ndarray | None = None


@dataclass(frozen=True)
class AbsoluteWaveformMetricInputs:
    """All values needed by the absolute waveform metric calculator."""

    beat_period_seconds: np.ndarray
    artery: VesselAbsoluteWaveformInputs
    vein: VesselAbsoluteWaveformInputs
