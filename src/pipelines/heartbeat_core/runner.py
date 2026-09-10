"""Build run-scoped heartbeat state shared by independent pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.heartbeat import (
    HeartbeatAnalysisResult,
    run_heartbeat_analysis,
)
from calculations.retinal_velocity import (
    VelocityEstimatorCacheKey,
    run_chunked_velocity_estimator,
    velocity_estimator_cache_key,
)
from utils.logger import Logger

from .scratch import heartbeat_scratch_h5
from .sources import load_heartbeat_inputs

HEARTBEAT_RESULT_STATE = "heartbeat.result"
_ARTERIAL_SIGNAL_CACHE_STATE = "heartbeat.cache.arterial_velocity_signal"
_ANALYSIS_CACHE_STATE = "heartbeat.cache.analysis"
_VELOCITY_ESTIMATION_CACHE_STATE = "heartbeat.cache.velocity_estimation"
_DEFAULT_LOWPASS_HZ = 15.0


@dataclass(frozen=True)
class HeartbeatResult:
    """Frame indexes delimiting complete heartbeat cycles."""

    cycle_boundary_indexes: np.ndarray
    index_base: int


@dataclass(frozen=True)
class _CachedVelocityEstimation:
    key: VelocityEstimatorCacheKey
    analysis: dict[str, object]


def run_heartbeat_core(ctx) -> HeartbeatResult:
    """Compute heartbeat boundaries once and place them in run-scoped state."""

    started = perf_counter()
    Logger.log("Starting shared heartbeat analysis (scratch=RAM)...")
    inputs = load_heartbeat_inputs(ctx)
    retain_velocity_video = _pipeline_scheduled(ctx, "waveform_velocity_core")
    velocity_video = (
        np.empty(tuple(int(size) for size in inputs.moment0.shape), dtype=np.float32)
        if retain_velocity_video
        else None
    )
    with heartbeat_scratch_h5(ctx) as scratch_h5:
        velocity = run_chunked_velocity_estimator(
            moment0=inputs.moment0,
            moment2=inputs.moment2,
            artery_mask=inputs.artery_mask,
            vein_mask=inputs.vein_mask,
            optic_disc_center=inputs.optic_disc_center,
            optic_disc_width=inputs.optic_disc_width,
            optic_disc_height=inputs.optic_disc_height,
            local_background_dist=inputs.local_background_dist,
            scratch_h5=scratch_h5,
            retain_velocity_video=retain_velocity_video,
            velocity_video_output=velocity_video,
        )
    if retain_velocity_video:
        ctx.state.set(
            _VELOCITY_ESTIMATION_CACHE_STATE,
            _CachedVelocityEstimation(
                key=_velocity_estimator_key(inputs),
                analysis=velocity,
            ),
        )
    arterial_signal = np.asarray(
        velocity["retinal_artery_velocity_signal"],
        dtype=np.float32,
    )
    analysis = run_heartbeat_analysis(
        arterial_signal,
        dt_seconds=float(inputs.timing.dt_seconds),
        lowpass_freq_hz=_DEFAULT_LOWPASS_HZ,
    )
    result = HeartbeatResult(
        cycle_boundary_indexes=np.asarray(
            analysis.systole.systole_indexes,
            dtype=np.int32,
        ),
        index_base=inputs.index_base,
    )
    ctx.state.set(HEARTBEAT_RESULT_STATE, result)
    ctx.state.set(_ARTERIAL_SIGNAL_CACHE_STATE, arterial_signal)
    ctx.state.set(_ANALYSIS_CACHE_STATE, analysis)
    Logger.log(
        f"Completed shared heartbeat analysis in {perf_counter() - started:.1f}s."
    )
    return result


def heartbeat_result(ctx) -> HeartbeatResult:
    """Return the beat boundaries produced by the declared upstream pipeline."""

    value = ctx.state.get(HEARTBEAT_RESULT_STATE)
    if not isinstance(value, HeartbeatResult):
        raise RuntimeError(
            "Shared heartbeat state is unavailable; check the pipeline DAG dependency."
        )
    return value


def cached_heartbeat_analysis(ctx) -> HeartbeatAnalysisResult:
    """Return private heartbeat details used by declared downstream consumers."""

    value = ctx.state.get(_ANALYSIS_CACHE_STATE)
    if not isinstance(value, HeartbeatAnalysisResult):
        raise RuntimeError("Cached heartbeat analysis is unavailable.")
    return value


def cached_velocity_estimation(ctx, source_data) -> dict[str, object] | None:
    """Return the heartbeat estimator result when every input still matches."""

    value = ctx.state.get(_VELOCITY_ESTIMATION_CACHE_STATE)
    if not isinstance(value, _CachedVelocityEstimation):
        return None
    if value.key != _velocity_estimator_key(source_data):
        return None
    velocity_map = value.analysis.get("velocity_map")
    if velocity_map is None:
        return None
    return dict(value.analysis)


def _velocity_estimator_key(source) -> VelocityEstimatorCacheKey:
    return velocity_estimator_cache_key(
        moment0=source.moment0,
        moment2=source.moment2,
        artery_mask=source.artery_mask
        if hasattr(source, "artery_mask")
        else source.retinal_artery_mask,
        vein_mask=source.vein_mask
        if hasattr(source, "vein_mask")
        else source.retinal_vein_mask,
        optic_disc_center=source.optic_disc_center,
        optic_disc_width=source.optic_disc_width,
        optic_disc_height=source.optic_disc_height,
        local_background_dist=source.local_background_dist,
    )


def _pipeline_scheduled(ctx, name: str) -> bool:
    predicate = getattr(ctx, "pipeline_scheduled", None)
    return bool(callable(predicate) and predicate(name))


__all__ = [
    "HEARTBEAT_RESULT_STATE",
    "HeartbeatResult",
    "cached_heartbeat_analysis",
    "cached_velocity_estimation",
    "heartbeat_result",
    "run_heartbeat_core",
]
