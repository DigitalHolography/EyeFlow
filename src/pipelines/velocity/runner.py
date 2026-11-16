"""Orchestrate dual velocity estimation and workflow-specific cardiac cycles."""

from __future__ import annotations

from copy import copy
from dataclasses import replace

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
    cardiac_cycles_from_available_vessel,
)
from data_manager import InvocationResource
from input_output.schema.eyeflow_output import (
    VELOCITY_WORKFLOW_FOLDERS,
    VELOCITY_WORKFLOW_ROOTS,
    processing_path,
)
from input_output.writers.h5_output import PipelineH5Output
from pipeline_engine import ExecutionVariant, ProductSpec, StageGraph, StageSpec
from pipeline_engine.stage_graph import InputAccess, InputUse
from pipeline_engine.context import require_managed_runtime
from pipelines.vessel_inputs import load_vessel_topology_inputs
from utils.logger import Logger
from velocity_calibration import physical_velocity_provenance

from .estimation import estimate_retinal_velocity
from .models import RetinalVelocity, RetinalVelocityData, VelocityWorkflows
from .outputs import pack_velocity_outputs
from .signal_processing import build_velocity, DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ

VELOCITY_STATE = "velocity"
VELOCITY_WORKFLOWS_STATE = "velocity_workflows"


def _inputs(ctx):
    Logger.log("Starting velocity core processing...")
    sources = {}
    failures = {}
    for method in getattr(ctx, "velocity_workflows", VELOCITY_WORKFLOW_ROOTS):
        method_ctx = copy(ctx)
        method_ctx.velocity_estimation_method = method
        try:
            sources[method] = load_vessel_topology_inputs(method_ctx)
        except Exception as exc:  # noqa: BLE001 - isolate a failed scientific workflow
            failures[method] = str(exc)
            Logger.log_warning(f"Skipping {method}: {exc}")

    if not sources:
        raise RuntimeError(f"No selected velocity workflow has usable inputs: {failures}")
    reference_shape = _source_shape(next(iter(sources.values())))
    ctx.dag.handoff("inputs", (
        {method: InvocationResource(source) for method, source in sources.items()},
        failures, reference_shape,
    ))


def _estimates(ctx):
    sources, initial_failures, reference_shape = ctx.dag.require("inputs")
    failures = dict(initial_failures)
    estimates = {method: _estimate_available(ctx, source, reference_shape, failures)
                 for method, source in sources.items()}
    ctx.dag.handoff("estimates", (estimates, failures))


def _cycles(ctx):
    sources, _, _ = ctx.dag.require("inputs")
    estimates, initial_failures = ctx.dag.require("estimates")
    failures = dict(initial_failures)
    variants = []
    for method, velocity_data in estimates.items():
        if velocity_data is None:
            continue
        source = sources[method]
        try:
            cycle_analysis, cycle_source = detect_velocity_cardiac_cycles(
                velocity_data,
                dt_seconds=float(source.holodoppler.timing.dt_seconds),
            )
            _log_cardiac_cycle_warnings(
                cycle_analysis,
                cycle_source,
                dt_seconds=float(source.holodoppler.timing.dt_seconds),
            )
            provenance = {
                **physical_velocity_provenance(
                    velocity_estimation_method=method,
                    band_ratio_frequency_scale_hz=source.band_ratio_frequency_scale_hz,
                ),
                **velocity_data.provenance,
                "cardiac_cycle_detection_method": method,
                "cardiac_cycle_detection_quantity": "physical_velocity",
                "cardiac_cycle_detection_unit": "mm/s",
                "cardiac_cycle_detection_source": cycle_source,
            }
            result = build_velocity(
                replace(velocity_data, provenance=provenance),
                cycle_analysis,
                cycle_source,
                dt_seconds=float(source.holodoppler.timing.dt_seconds),
            )
            root = VELOCITY_WORKFLOW_ROOTS[method]
            variants.append(
                ExecutionVariant(
                    name=method,
                    state={VELOCITY_STATE: result, "velocity_source": InvocationResource(source)},
                    output_namespace=root,
                    artifact_namespace=VELOCITY_WORKFLOW_FOLDERS[method],
                    provenance=provenance,
                )
            )
        except Exception as exc:
            failures[method] = str(exc)
            Logger.log_warning(f"Skipping {method}: {exc}")
    ctx.dag.handoff("workflows", (variants, failures))


def _publish(ctx):
    candidates, initial_failures = ctx.dag.require("workflows")
    failures, variants = dict(initial_failures), []
    for variant in candidates:
        variant = replace(variant, state=dict(variant.state))
        try:
            packed = pack_velocity_outputs(variant.state[VELOCITY_STATE])
            qualified = {processing_path(path, variant.output_namespace): value for path, value in packed.items()}
            output = PipelineH5Output(ctx.runtime.work_h5, processing_root=variant.output_namespace,
                                      provenance=variant.provenance)
            output.set_attrs(variant.provenance)
            ctx.dag.handoff(variant.name + "_outputs", qualified)
        except Exception as exc:
            failures[variant.name] = str(exc)
            Logger.log_warning(f"Skipping {variant.name}: {exc}")
        else:
            variants.append(variant)
    if not variants:
        raise RuntimeError(f"No selected velocity workflow completed: {failures}")
    results = {variant.name: variant.state[VELOCITY_STATE] for variant in variants}
    ctx.state.set(VELOCITY_WORKFLOWS_STATE, VelocityWorkflows(
        moments=results.get("doppler_moments"),
        bandratio=results.get("frequency_bands"),
    ))
    ctx.publish_execution_variants(variants, failures=failures)
    canonical = variants[0]
    for key, value in canonical.state.items():
        ctx.state.set(key, value)


def velocity_graph(*, options=()):
    return StageGraph((
        StageSpec("inputs", _inputs, products=(ProductSpec("inputs"),), inputs=(
            InputUse("hd.metadata", InputAccess.SMALL),
            InputUse("dv.vessels", InputAccess.SMALL),
            InputUse("dv.optional_segmentation", InputAccess.SMALL, optional=True),
        )),
        StageSpec("estimates", _estimates, requires=("inputs",), products=(ProductSpec("estimates"),), inputs=(
            InputUse("hd.moment0", InputAccess.SEQUENTIAL, when="doppler_moments"),
            InputUse("hd.moment2", InputAccess.SEQUENTIAL, when="doppler_moments"),
            InputUse("hd.band_low", InputAccess.SEQUENTIAL, when="frequency_bands"),
            InputUse("hd.band_high", InputAccess.SEQUENTIAL, when="frequency_bands"),
        )),
        StageSpec("cardiac_cycles", _cycles, requires=("inputs", "estimates"),
                  products=(ProductSpec("workflows"),)),
        StageSpec("publish", _publish, requires=("workflows",), products=tuple(
            ProductSpec(method + "_outputs", pack=lambda value, ctx: value, optional=True)
            for method in VELOCITY_WORKFLOW_ROOTS)),
    ))


def run_velocity(ctx):
    require_managed_runtime(ctx, "velocity")
    return velocity_graph(options=ctx.options_for("velocity"))


def _source_shape(source):
    images = source.image_maps
    reference = (
        images.moment0 if source.velocity_estimation_method == "doppler_moments" else images.band_lf
    )
    return tuple(reference.shape)


def _estimate_available(ctx, source, reference_shape, failures):
    if source is None:
        return None
    try:
        if _source_shape(source) != reference_shape:
            raise ValueError("Both workflows must share (frame, y, x) dimensions.")
        return _estimate(ctx, source)
    except Exception as exc:  # noqa: BLE001 - isolate a failed scientific workflow
        failures[source.velocity_estimation_method] = str(exc)
        Logger.log_warning(f"Skipping {source.velocity_estimation_method}: {exc}")
        return None


def _estimate(ctx, source):
    images = source.image_maps
    segmentation = source.segmentation
    retain_velocity_map = _pipeline_scheduled(ctx, "velocity_analysis")
    scratch = getattr(getattr(ctx, "runtime", None), "scratch", None)
    return estimate_retinal_velocity(
        image_maps=images,
        velocity_estimation_method=source.velocity_estimation_method,
        band_ratio_frequency_scale_hz=source.band_ratio_frequency_scale_hz,
        artery_mask=segmentation.vessels.artery,
        vein_mask=segmentation.vessels.vein,
        background_mask=segmentation.vessels.velocity_background,
        optic_disc_center=segmentation.optic_disc.center,
        local_background_dist=source.doppler_view.local_background_dist,
        retain_velocity_video=retain_velocity_map,
        velocity_video_allocator=(
            scratch.empty_float32 if retain_velocity_map and scratch is not None else None
        ),
    )


def velocity(ctx) -> RetinalVelocity:
    """Return the canonical result produced by the DAG dependency."""

    value = ctx.state.get(VELOCITY_STATE)
    if not isinstance(value, RetinalVelocity):
        raise RuntimeError("Velocity state is unavailable; check the pipeline DAG dependency.")
    return value


def velocity_workflows(ctx) -> VelocityWorkflows:
    """Return successful workflow results without selecting a canonical method."""
    value = ctx.state.get(VELOCITY_WORKFLOWS_STATE)
    if value is None:
        raise RuntimeError("Velocity workflow state is unavailable; check the pipeline DAG dependency.")
    if not isinstance(value, VelocityWorkflows):
        raise TypeError("Velocity workflow state must be a VelocityWorkflows result.")
    return value


def cardiac_cycle_indexes(ctx) -> np.ndarray:
    """Return canonical zero-based frame indexes delimiting cardiac cycles."""

    return velocity(ctx).cycle_boundary_indexes


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
    "VELOCITY_STATE",
    "cardiac_cycle_indexes",
    "cardiac_cycles",
    "run_velocity",
    "velocity",
]


def detect_velocity_cardiac_cycles(
    velocity_data: RetinalVelocityData,
    *,
    dt_seconds: float,
) -> tuple[CardiacCycleAnalysis, str]:
    """Detect cycles from signed, background-corrected velocity in mm/s.

    The detector owns low-pass filtering; pass the estimator's unfiltered
    whole-vessel signals to avoid filtering them twice.
    """
    return cardiac_cycles_from_available_vessel(
        velocity_data.artery.velocity,
        velocity_data.vein.velocity,
        dt_seconds=float(dt_seconds),
        lowpass_freq_hz=DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
    )
