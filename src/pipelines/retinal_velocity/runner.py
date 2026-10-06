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
    velocity_source = _velocity_source_volume(source)
    velocity_map_output = (
        np.empty(tuple(int(size) for size in velocity_source.shape), dtype=np.float32)
        if retain_velocity_map
        else None
    )
    with retinal_velocity_scratch_h5(ctx) as scratch_h5:
        velocity_data = estimate_retinal_velocity(
            moment0=images.moment0,
            moment2=images.moment2,
            band_lf=images.band_lf,
            band_hf=images.band_hf,
            velocity_estimation_method=source.velocity_estimation_method,
            band_ratio_frequency_scale_hz=(
                source.band_ratio_frequency_scale_hz
            ),
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
        velocity_data,
        dt_seconds=float(timing.dt_seconds),
    )
    _log_cardiac_cycle_warnings(
        cycle_analysis,
        cycle_source,
        dt_seconds=float(timing.dt_seconds),
    )
    velocity = build_retinal_velocity(
        velocity_data,
        cycle_analysis,
        cycle_source,
        dt_seconds=float(timing.dt_seconds),
    )
    ctx.state.set(RETINAL_VELOCITY_STATE, velocity)
    Logger.log(
        f"Completed retinal velocity core processing in {perf_counter() - started:.1f}s."
    )
    return velocity, pack_retinal_velocity_outputs(velocity)


def _velocity_source_volume(source):
    """Return the primary volume for the configured estimator."""

    images = source.image_maps
    sources = {
        "doppler_moments": images.moment0,
        "frequency_bands": images.band_lf,
    }
    try:
        volume = sources[source.velocity_estimation_method]
    except KeyError as exc:
        raise ValueError(
            "Unsupported velocity_estimation_method "
            f"{source.velocity_estimation_method!r}."
        ) from exc
    if volume is None:
        raise ValueError(
            f"The {source.velocity_estimation_method!r} estimator has no primary input."
        )
    return volume


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


def _log_cardiac_cycle_warnings(
    analysis: CardiacCycleAnalysis,
    source: str,
    *,
    dt_seconds: float,
) -> None:
    """Record retained gaps that look like two or three cardiac periods."""

    if source == "none":
        Logger.log_warning(
            "No usable systole sequence was detected in the artery or vein; "
            "the full recording is retained as one fallback cardiac cycle."
        )
        return

    for gap in getattr(analysis.systole, "suspected_missed_beat_gaps", ()):
        Logger.log_warning(
            "Possible missed systole in the "
            f"{source} signal: frames {gap.start_index} to {gap.stop_index} "
            f"span {gap.estimated_multiple}x the estimated cardiac period "
            f"({gap.interval_samples * dt_seconds:.3f}s versus "
            f"{gap.estimated_period_samples * dt_seconds:.3f}s). "
            "The interval was retained for per-beat analysis."
        )


__all__ = [
    "RETINAL_VELOCITY_STATE",
    "cardiac_cycle_indexes",
    "cardiac_cycles",
    "retinal_velocity",
    "run_retinal_velocity",
]
