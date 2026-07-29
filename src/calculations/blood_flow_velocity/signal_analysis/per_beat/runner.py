"""Per-beat blood-flow velocity calculations from BloodFlowVelocity/perBeatAnalysis.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.heartbeat import (
    SpectralHeartbeatResult,
)

from .segments import (
    PerBeatSegmentAnalysisResult,
    per_beat_segment_analysis,
)
from .signal import PerBeatSignalAnalysisResult, per_beat_signal_analysis
from ._signal_utils import normalize_cycle_boundaries


@dataclass(frozen=True)
class VesselPerBeatAnalysisResult:
    signal: PerBeatSignalAnalysisResult
    vmax_band_limited: np.ndarray
    vmin_band_limited: np.ndarray
    vti_per_beat: np.ndarray
    segments: PerBeatSegmentAnalysisResult | None = None


@dataclass(frozen=True)
class PerBeatAnalysisInput:
    arterial_velocity_signal: np.ndarray
    venous_velocity_signal: np.ndarray
    cycle_boundary_indexes: np.ndarray
    band_limited_signal_harmonic_count: int
    heartbeat: SpectralHeartbeatResult
    dt_seconds: float
    arterial_velocity_segments: np.ndarray | None = None
    venous_velocity_segments: np.ndarray | None = None
    index_base: int | None = None


@dataclass(frozen=True)
class PerBeatAnalysisResult:
    beat_period_seconds: np.ndarray
    heartbeat: SpectralHeartbeatResult
    cycle_boundary_indexes: np.ndarray
    artery: VesselPerBeatAnalysisResult
    vein: VesselPerBeatAnalysisResult


def run_per_beat_analysis(inputs: PerBeatAnalysisInput) -> PerBeatAnalysisResult:
    cycle_boundaries = normalize_cycle_boundaries(
        inputs.cycle_boundary_indexes,
        np.asarray(inputs.arterial_velocity_signal).size,
        index_base=inputs.index_base,
    )
    return PerBeatAnalysisResult(
        beat_period_seconds=_beat_period_seconds(
            cycle_boundaries,
            inputs.dt_seconds,
        ),
        heartbeat=inputs.heartbeat,
        cycle_boundary_indexes=cycle_boundaries,
        vein=_run_vessel(
            inputs.venous_velocity_signal,
            inputs.venous_velocity_segments,
            cycle_boundaries,
            inputs,
        ),
        artery=_run_vessel(
            inputs.arterial_velocity_signal,
            inputs.arterial_velocity_segments,
            cycle_boundaries,
            inputs,
        ),
    )


def _beat_period_seconds(
    cycle_boundaries: np.ndarray,
    dt_seconds: float,
) -> np.ndarray:
    return (
        np.diff(cycle_boundaries).astype(np.float32)
        * np.float32(dt_seconds)
    ).astype(np.float32, copy=False)


def _run_vessel(
    velocity_signal,
    velocity_segments,
    cycle_boundaries: np.ndarray,
    inputs: PerBeatAnalysisInput,
) -> VesselPerBeatAnalysisResult:
    signal = per_beat_signal_analysis(
        velocity_signal,
        cycle_boundaries,
        inputs.band_limited_signal_harmonic_count,
        index_base=0,
    )
    return VesselPerBeatAnalysisResult(
        signal=signal,
        vmax_band_limited=np.max(signal.velocity_signal_per_beat_band_limited, axis=1),
        vmin_band_limited=np.min(signal.velocity_signal_per_beat_band_limited, axis=1),
        vti_per_beat=(
            np.sum(signal.velocity_signal_per_beat, axis=1)
            * np.float32(inputs.dt_seconds)
        ).astype(np.float32, copy=False),
        segments=_run_segments(velocity_segments, cycle_boundaries, inputs),
    )


def _run_segments(
    velocity_segments,
    cycle_boundaries: np.ndarray,
    inputs: PerBeatAnalysisInput,
) -> PerBeatSegmentAnalysisResult | None:
    if velocity_segments is None:
        return None
    return per_beat_segment_analysis(
        velocity_segments,
        cycle_boundaries,
        inputs.band_limited_signal_harmonic_count,
        index_base=0,
    )
