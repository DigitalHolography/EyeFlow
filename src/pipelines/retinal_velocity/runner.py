"""Orchestrate shared retinal-velocity core processing."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from utils.logger import Logger

from .cardiac_cycle import detect_cardiac_cycles
from .estimation import estimate_retinal_velocity
from .models import RetinalVelocity
from .outputs import pack_retinal_velocity_outputs
from .scratch import retinal_velocity_scratch_h5
from .signal_processing import build_retinal_velocity
from .sources import load_retinal_velocity_inputs

RETINAL_VELOCITY_STATE = "retinal_velocity"


def run_retinal_velocity(ctx) -> tuple[RetinalVelocity, dict[str, object]]:
    """Compute retinal velocity once and publish reusable typed state."""

    started = perf_counter()
    Logger.log("Starting retinal velocity core processing (scratch=RAM)...")
    source = load_retinal_velocity_inputs(ctx)
    images = source.image_maps
    segmentation = source.segmentation
    timing = source.holodoppler.timing
    retain_velocity_map = _pipeline_scheduled(ctx, "waveform_velocity")
    velocity_map_output = (
        np.empty(tuple(int(size) for size in images.moment0.shape), dtype=np.float32)
        if retain_velocity_map
        else None
    )
    with retinal_velocity_scratch_h5(ctx) as scratch_h5:
        estimation = estimate_retinal_velocity(
            moment0=images.moment0,
            moment2=images.moment2,
            artery_mask=segmentation.vessels.artery,
            vein_mask=segmentation.vessels.vein,
            background_mask=segmentation.vessels.velocity_background,
            optic_disc_center=segmentation.optic_disc.center,
            local_background_dist=source.doppler_view.local_background_dist,
            scratch_h5=scratch_h5,
            retain_velocity_video=retain_velocity_map,
            velocity_video_output=velocity_map_output,
        )
    cycle_analysis, cycle_source = detect_cardiac_cycles(
        estimation,
        dt_seconds=float(timing.dt_seconds),
    )
    velocity = build_retinal_velocity(
        estimation,
        cycle_analysis,
        cycle_source,
        dt_seconds=float(timing.dt_seconds),
    )
    ctx.state.set(RETINAL_VELOCITY_STATE, velocity)
    Logger.log(
        f"Completed retinal velocity core processing in {perf_counter() - started:.1f}s."
    )
    return velocity, pack_retinal_velocity_outputs(velocity)


def retinal_velocity(ctx) -> RetinalVelocity:
    """Return the canonical result produced by the DAG dependency."""

    value = ctx.state.get(RETINAL_VELOCITY_STATE)
    if not isinstance(value, RetinalVelocity):
        raise RuntimeError(
            "Retinal velocity state is unavailable; check the pipeline DAG dependency."
        )
    return value


def cardiac_cycles(ctx) -> CardiacCycleAnalysis:
    """Return reusable cardiac-cycle timing without exposing pipeline internals."""

    return retinal_velocity(ctx).cardiac_cycle


def cardiac_cycle_indexes(ctx) -> np.ndarray:
    """Return canonical zero-based frame indexes delimiting cardiac cycles."""

    return retinal_velocity(ctx).cycle_boundary_indexes


def _pipeline_scheduled(ctx, name: str) -> bool:
    predicate = getattr(ctx, "pipeline_scheduled", None)
    return bool(callable(predicate) and predicate(name))


__all__ = [
    "RETINAL_VELOCITY_STATE",
    "cardiac_cycle_indexes",
    "cardiac_cycles",
    "retinal_velocity",
    "run_retinal_velocity",
]
