"""Orchestrate dual velocity estimation and workflow-specific cardiac cycles."""

from __future__ import annotations

from copy import copy
from dataclasses import replace
from time import perf_counter

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from input_output.h5_access import PipelineH5Output
from input_output.schema.eyeflow_output import (
    VELOCITY_WORKFLOW_FOLDERS,
    VELOCITY_WORKFLOW_ROOTS,
    processing_path,
)
from pipeline_engine.imports import ExecutionVariant
from utils.logger import Logger
from velocity_calibration import physical_velocity_provenance

from .cardiac_cycle import detect_source_cardiac_cycles
from .estimation import estimate_retinal_velocity
from .models import RetinalVelocity
from .outputs import pack_velocity_outputs
from .signal_processing import build_velocity
from .sources import load_velocity_inputs

VELOCITY_STATE = "velocity"
VELOCITY_FAILURES_STATE = "velocity_workflow_failures"


def run_velocity(ctx) -> tuple[RetinalVelocity, dict[str, object]]:
    """Compute both methods with independent cycles and isolated workflow state."""

    started = perf_counter()
    Logger.log("Starting velocity core processing...")
    sources = {}
    failures = {}
    ctx.state.set(VELOCITY_FAILURES_STATE, failures)
    for method in VELOCITY_WORKFLOW_ROOTS:
        method_ctx = copy(ctx)
        method_ctx.velocity_estimation_method = method
        try:
            sources[method] = load_velocity_inputs(method_ctx)
        except Exception as exc:  # noqa: BLE001 - isolate a failed scientific workflow
            failures[method] = str(exc)
            Logger.log_warning(f"Skipping {method}: {exc}")

    cycles = {}
    reference_shape = None
    for method, source in tuple(sources.items()):
        try:
            cycle_analysis, cycle_source = detect_source_cardiac_cycles(source)
            _log_cardiac_cycle_warnings(
                cycle_analysis,
                cycle_source,
                dt_seconds=float(source.holodoppler.timing.dt_seconds),
            )
            if reference_shape is None:
                reference_shape = _source_shape(source)
            cycles[method] = (cycle_analysis, cycle_source)
        except Exception as exc:  # noqa: BLE001 - isolate a failed scientific workflow
            failures[method] = str(exc)
            del sources[method]
            Logger.log_warning(f"Skipping {method}: {exc}")
    if not cycles:
        _write_failures(ctx, failures)
        raise RuntimeError(f"Neither velocity workflow has usable inputs: {failures}")

    variants: list[ExecutionVariant] = []
    outputs = {}
    velocity_data_moments = _estimate_available(
        ctx,
        sources.get("doppler_moments"),
        reference_shape,
        failures,
    )
    velocity_data_bandratio = _estimate_available(
        ctx,
        sources.get("frequency_bands"),
        reference_shape,
        failures,
    )
    for method, velocity_data in (
        ("doppler_moments", velocity_data_moments),
        ("frequency_bands", velocity_data_bandratio),
    ):
        if velocity_data is None:
            continue
        source = sources[method]
        try:
            cycle_analysis, cycle_source = cycles[method]
            provenance = {
                **physical_velocity_provenance(
                    velocity_estimation_method=method,
                    band_ratio_frequency_scale_hz=source.band_ratio_frequency_scale_hz,
                ),
                **velocity_data.provenance,
                "cardiac_cycle_detection_method": method,
                "cardiac_cycle_detection_quantity": "raw_rms_frequency",
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
                    state={VELOCITY_STATE: result, "velocity_source": source},
                    output_namespace=root,
                    artifact_namespace=VELOCITY_WORKFLOW_FOLDERS[method],
                    provenance=provenance,
                )
            )
            packed = pack_velocity_outputs(result)
            outputs.update(
                {
                    processing_path(path, root): value
                    for path, value in packed.items()
                }
            )
            if hasattr(ctx, "output"):
                output = PipelineH5Output(
                    ctx.runtime.work_h5, processing_root=root, provenance=provenance
                )
                output.set_attrs(provenance)
                output.write_many(packed)
        except Exception as exc:  # noqa: BLE001 - isolate a failed scientific workflow
            failures[method] = str(exc)
            variants = [variant for variant in variants if variant.name != method]
            root = VELOCITY_WORKFLOW_ROOTS[method]
            outputs = {
                path: value for path, value in outputs.items() if not path.startswith(root + "/")
            }
            if hasattr(ctx, "output") and root in ctx.runtime.work_h5:
                del ctx.runtime.work_h5[root]
            Logger.log_warning(f"Skipping {method}: {exc}")
    ctx.state.set(VELOCITY_FAILURES_STATE, failures)
    _write_failures(ctx, failures)
    if not variants:
        raise RuntimeError(f"Neither velocity workflow completed: {failures}")
    publish_variants = getattr(ctx, "publish_execution_variants", None)
    if callable(publish_variants):
        publish_variants(variants, failures=failures)
    canonical = variants[0]
    for key, value in canonical.state.items():
        ctx.state.set(key, value)
    Logger.log(f"Completed velocity core processing in {perf_counter() - started:.1f}s.")
    return canonical.state[VELOCITY_STATE], outputs


def _write_failures(ctx, failures):
    import json

    if hasattr(ctx, "output"):
        ctx.output.h5.set_attr("velocity_workflow_failures", json.dumps(failures, sort_keys=True))


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
    )


def velocity(ctx) -> RetinalVelocity:
    """Return the canonical result produced by the DAG dependency."""

    value = ctx.state.get(VELOCITY_STATE)
    if not isinstance(value, RetinalVelocity):
        raise RuntimeError("Velocity state is unavailable; check the pipeline DAG dependency.")
    return value


def cardiac_cycles(ctx) -> CardiacCycleAnalysis:
    """Return reusable cardiac-cycle timing without exposing pipeline internals."""

    return velocity(ctx).cardiac_cycle


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
