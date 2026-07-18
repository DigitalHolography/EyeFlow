"""Per-beat blood-flow velocity calculations from BloodFlowVelocity/perBeatAnalysis.m."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .segments import (
    PerBeatSegmentAnalysisResult,
    per_beat_segment_analysis,
)
from .quality import BeatQualityResult, BeatQualitySettings, assess_beat_quality
from .signal import PerBeatSignalAnalysisResult, per_beat_signal_analysis


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
    interval_period_ratio: np.ndarray | None = None
    recovered_boundary_indexes: np.ndarray | None = None
    nominal_period_samples: float | None = None
    quality_settings: BeatQualitySettings | None = None


@dataclass(frozen=True)
class PerBeatAnalysisResult:
    beat_period_idx: np.ndarray
    beat_period_seconds: np.ndarray
    artery: VesselPerBeatAnalysisResult
    vein: VesselPerBeatAnalysisResult
    candidate_start_indexes: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    candidate_stop_indexes: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    candidate_beat_period_idx: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    candidate_beat_period_seconds: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float32)
    )
    recovered_boundary_mask: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.bool_)
    )
    nominal_period_samples: np.float32 = np.float32(np.nan)
    index_base: int | None = None
    quality: BeatQualityResult | None = None


def run_per_beat_analysis(inputs: PerBeatAnalysisInput) -> PerBeatAnalysisResult:
    sys_idx_list = np.asarray(
        inputs.systolic_acceleration_peak_indexes,
        dtype=np.int32,
    ).reshape(-1)
    candidate_period_idx = np.diff(sys_idx_list).astype(np.int32, copy=False)
    candidate_period_seconds = _beat_period_seconds(inputs, candidate_period_idx)
    quality = _quality_result(inputs, sys_idx_list)
    accepted_mask = (
        quality.accepted_mask
        if quality is not None
        else np.ones(candidate_period_idx.size, dtype=np.bool_)
    )
    if not np.any(accepted_mask):
        raise ValueError("No beat intervals passed waveform quality control.")

    return PerBeatAnalysisResult(
        beat_period_idx=candidate_period_idx[accepted_mask],
        beat_period_seconds=candidate_period_seconds[accepted_mask],
        vein=_run_vessel(
            inputs.venous_velocity_signal,
            inputs.venous_velocity_segments,
            sys_idx_list,
            inputs,
            interval_mask=accepted_mask,
        ),
        artery=_run_vessel(
            inputs.arterial_velocity_signal,
            inputs.arterial_velocity_segments,
            sys_idx_list,
            inputs,
            interval_mask=accepted_mask,
        ),
        candidate_start_indexes=sys_idx_list[:-1],
        candidate_stop_indexes=sys_idx_list[1:],
        candidate_beat_period_idx=candidate_period_idx,
        candidate_beat_period_seconds=candidate_period_seconds,
        recovered_boundary_mask=_recovered_boundary_mask(inputs, sys_idx_list),
        nominal_period_samples=np.float32(
            np.nan
            if inputs.nominal_period_samples is None
            else inputs.nominal_period_samples
        ),
        index_base=inputs.index_base,
        quality=quality,
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
    *,
    interval_mask: np.ndarray | None = None,
) -> VesselPerBeatAnalysisResult:
    signal = per_beat_signal_analysis(
        velocity_signal,
        sys_idx_list,
        inputs.band_limited_signal_harmonic_count,
        index_base=inputs.index_base,
        interval_mask=interval_mask,
    )
    return VesselPerBeatAnalysisResult(
        signal=signal,
        vmax_band_limited=np.max(signal.velocity_signal_per_beat_band_limited, axis=1),
        vmin_band_limited=np.min(signal.velocity_signal_per_beat_band_limited, axis=1),
        segments=_run_segments(
            velocity_segments,
            sys_idx_list,
            inputs,
            interval_mask=interval_mask,
        ),
    )


def _run_segments(
    velocity_segments,
    sys_idx_list: np.ndarray,
    inputs: PerBeatAnalysisInput,
    *,
    interval_mask: np.ndarray | None = None,
) -> PerBeatSegmentAnalysisResult | None:
    if velocity_segments is None:
        return None
    return per_beat_segment_analysis(
        velocity_segments,
        sys_idx_list,
        inputs.band_limited_signal_harmonic_count,
        index_base=inputs.index_base,
        interval_mask=interval_mask,
    )


def _quality_result(
    inputs: PerBeatAnalysisInput,
    sys_idx_list: np.ndarray,
) -> BeatQualityResult | None:
    if inputs.quality_settings is None:
        return None
    artery = per_beat_signal_analysis(
        inputs.arterial_velocity_signal,
        sys_idx_list,
        inputs.band_limited_signal_harmonic_count,
        index_base=inputs.index_base,
    )
    vein = per_beat_signal_analysis(
        inputs.venous_velocity_signal,
        sys_idx_list,
        inputs.band_limited_signal_harmonic_count,
        index_base=inputs.index_base,
    )
    period_ratio = inputs.interval_period_ratio
    if period_ratio is None:
        period_ratio = _period_ratios(inputs, sys_idx_list)
    return assess_beat_quality(
        artery,
        vein,
        period_ratio=period_ratio,
        settings=inputs.quality_settings,
    )


def _period_ratios(
    inputs: PerBeatAnalysisInput,
    sys_idx_list: np.ndarray,
) -> np.ndarray:
    periods = np.diff(sys_idx_list).astype(np.float32, copy=False)
    nominal = inputs.nominal_period_samples
    if nominal is None or not np.isfinite(nominal) or nominal <= 0:
        return np.ones(periods.size, dtype=np.float32)
    return (periods / np.float32(nominal)).astype(np.float32, copy=False)


def _recovered_boundary_mask(
    inputs: PerBeatAnalysisInput,
    sys_idx_list: np.ndarray,
) -> np.ndarray:
    if inputs.recovered_boundary_indexes is None:
        return np.zeros(sys_idx_list.size, dtype=np.bool_)
    recovered = np.asarray(inputs.recovered_boundary_indexes, dtype=np.int32).reshape(-1)
    return np.isin(sys_idx_list, recovered)

