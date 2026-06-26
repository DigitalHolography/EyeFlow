"""Typed inputs for the AngioEye waveform-shape calculations."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class VesselWaveformInputs:
    """Optional global and segment waveforms for one vessel type."""

    raw_global: np.ndarray | None = None
    bandlimited_global: np.ndarray | None = None
    raw_segments: np.ndarray | None = None
    bandlimited_segments: np.ndarray | None = None


@dataclass(frozen=True)
class WaveformShapeMetricInputs:
    """All values needed by the AngioEye waveform-shape calculator."""

    beat_period_seconds: np.ndarray
    artery: VesselWaveformInputs
    vein: VesselWaveformInputs
