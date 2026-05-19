"""Per-beat blood-flow velocity calculations from BloodFlowVelocity/perBeatAnalysis.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..segments.per_beat_segments import (
    PerBeatSegmentAnalysisResult,
    per_beat_segment_analysis,
)
from .per_beat_signal import PerBeatSignalAnalysisResult, per_beat_signal_analysis


@dataclass(frozen=True)
class VesselPerBeatAnalysisResult:
    signal: PerBeatSignalAnalysisResult
    vmax_band_limited: np.ndarray
    vmin_band_limited: np.ndarray
    segments: PerBeatSegmentAnalysisResult | None = None


@dataclass(frozen=True)
class PerBeatAnalysisInput:
    arterial_velocity_signal: np.ndarray
    venous_velocity_signal: np.ndarray
    systolic_acceleration_peak_indexes: np.ndarray
    band_limited_signal_harmonic_count: int
    dt_seconds: float
    arterial_velocity_segments: np.ndarray | None = None
    venous_velocity_segments: np.ndarray | None = None
    beat_period_seconds: np.ndarray | None = None
    index_base: int | None = None


@dataclass(frozen=True)
class PerBeatAnalysisResult:
    beat_period_idx: np.ndarray
    beat_period_seconds: np.ndarray
    artery: VesselPerBeatAnalysisResult
    vein: VesselPerBeatAnalysisResult


def run_per_beat_analysis(inputs: PerBeatAnalysisInput) -> PerBeatAnalysisResult:
    sys_idx_list = np.asarray(
        inputs.systolic_acceleration_peak_indexes,
        dtype=np.int32,
    ).reshape(-1)
    beat_period_idx = np.diff(sys_idx_list).astype(np.int32, copy=False)
    return PerBeatAnalysisResult(
        beat_period_idx=beat_period_idx,
        beat_period_seconds=_beat_period_seconds(inputs, beat_period_idx),
        vein=_run_vessel(
            inputs.venous_velocity_signal,
            inputs.venous_velocity_segments,
            sys_idx_list,
            inputs,
        ),
        artery=_run_vessel(
            inputs.arterial_velocity_signal,
            inputs.arterial_velocity_segments,
            sys_idx_list,
            inputs,
        ),
    )


def _beat_period_seconds(
    inputs: PerBeatAnalysisInput,
    beat_period_idx: np.ndarray,
) -> np.ndarray:
    if inputs.beat_period_seconds is not None:
        periods = np.asarray(inputs.beat_period_seconds, dtype=np.float32).reshape(-1)
        if periods.size == beat_period_idx.size:
            return periods
    periods = beat_period_idx.astype(np.float32, copy=False) * np.float32(
        inputs.dt_seconds
    )
    return periods.astype(np.float32, copy=False)


def _run_vessel(
    velocity_signal,
    velocity_segments,
    sys_idx_list: np.ndarray,
    inputs: PerBeatAnalysisInput,
) -> VesselPerBeatAnalysisResult:
    signal = per_beat_signal_analysis(
        velocity_signal,
        sys_idx_list,
        inputs.band_limited_signal_harmonic_count,
        index_base=inputs.index_base,
    )
    return VesselPerBeatAnalysisResult(
        signal=signal,
        vmax_band_limited=np.max(signal.velocity_signal_per_beat_band_limited, axis=1),
        vmin_band_limited=np.min(signal.velocity_signal_per_beat_band_limited, axis=1),
        segments=_run_segments(velocity_segments, sys_idx_list, inputs),
    )


def _run_segments(
    velocity_segments,
    sys_idx_list: np.ndarray,
    inputs: PerBeatAnalysisInput,
) -> PerBeatSegmentAnalysisResult | None:
    if velocity_segments is None:
        return None
    return per_beat_segment_analysis(
        velocity_segments,
        sys_idx_list,
        inputs.band_limited_signal_harmonic_count,
        index_base=inputs.index_base,
    )

