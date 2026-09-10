"""Orchestrate selectable velocity output products."""

from time import perf_counter

import numpy as np

from input_output import EyeFlowOutputPaths
from pipelines.waveform_velocity_core.per_beat import run_velocity_per_beat_metrics
from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
    VELOCITY_PER_BEAT_RESULT_STATE,
    WAVEFORM_CONTEXT_STATE,
)
from utils.logger import Logger

from .continuous import (
    pack_continuous_velocity_outputs,
    pack_segment_velocity_outputs,
)
from .profiles import (
    pack_blood_volume_rate_outputs,
    pack_cross_section_displacement_profile_outputs,
    pack_cross_section_profile_outputs,
    pack_displacement_magnitude_outputs,
    pack_displacement_profile_outputs,
    pack_velocity_profile_fft_outputs,
)
from .quadrants import pack_quadrant_velocity_outputs
from .segment_maps import (
    pack_displacement_segment_map_outputs,
    pack_segment_map_outputs,
    prepare_segment_velocity_maps_per_beat,
)
from .segment_velocity_map_avi import export_segment_velocity_map_avis
from .spatial_gradient_profiles import (
    cleanup_spatial_gradient_artifacts,
    extract_spatial_gradient_segments,
    pack_spatial_gradient_profile_outputs,
)


def run_waveform_velocity(ctx) -> dict[str, object]:
    """Publish base velocity plus the selected derived velocity products."""
    context = _required_state(ctx, WAVEFORM_CONTEXT_STATE)
    selected = ctx.options_for("waveform_velocity")
    velocity_analysis = getattr(
        context,
        "velocity_analysis",
        getattr(context, "dopplerview_analysis", None),
    )
    metrics = pack_continuous_velocity_outputs(velocity_analysis)
    segments_selected = "segments" in selected
    maps_selected = "segment_velocity_maps" in selected
    profiles_selected = "velocity_profiles" in selected
    artery_velocity_maps_per_beat = None
    vein_velocity_maps_per_beat = None
    if maps_selected:
        map_started = perf_counter()
        Logger.log("Starting shared per-beat segment velocity-map interpolation...")
        artery_velocity_maps_per_beat, vein_velocity_maps_per_beat = (
            prepare_segment_velocity_maps_per_beat(
                context.artery_segment_result,
                context.vein_segment_result,
                context.per_beat_analysis.cycle_boundary_indexes,
                index_base=int(context.source_data.provenance["beat_index_base"]),
            )
        )
        Logger.log(
            "Completed shared per-beat segment velocity-map interpolation in "
            f"{perf_counter() - map_started:.1f}s."
        )
    if segments_selected:
        metrics.update(
            pack_segment_velocity_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                source_data=context.source_data,
            )
        )
    if maps_selected:
        segment_map_outputs = pack_segment_map_outputs(
            context.artery_segment_result,
            context.vein_segment_result,
            artery_velocity_maps_per_beat,
            vein_velocity_maps_per_beat,
        )
        metrics.update(segment_map_outputs)
        metrics.update(
            pack_displacement_segment_map_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                context.per_beat_analysis.cycle_boundary_indexes,
                index_base=int(
                    context.source_data.provenance["beat_index_base"]
                ),
            )
        )
        output = getattr(ctx, "output", None)
        if getattr(output, "available", False):
            avi_started = perf_counter()
            Logger.log("Starting segment velocity-map AVI export...")
            export_segment_velocity_map_avis(
                output,
                context.artery_segment_result,
                context.vein_segment_result,
                segment_map_outputs,
            )
            Logger.log(
                "Completed segment velocity-map AVI export in "
                f"{perf_counter() - avi_started:.1f}s."
            )

    per_beat_result = ctx.state.get(VELOCITY_PER_BEAT_RESULT_STATE)
    velocity_outputs = ctx.state.get(VELOCITY_PER_BEAT_OUTPUTS_STATE, {})
    if "per_beat" in selected or ctx.pipeline_scheduled("pdf_report"):
        if per_beat_result is None:
            per_beat_result, velocity_outputs = run_velocity_per_beat_metrics(context)
            ctx.state.set(VELOCITY_PER_BEAT_RESULT_STATE, per_beat_result)
            ctx.state.set(VELOCITY_PER_BEAT_OUTPUTS_STATE, velocity_outputs)
        if segments_selected:
            metrics.update(velocity_outputs)
        else:
            schema = EyeFlowOutputPaths.active()
            segment_paths = {
                schema.artery_per_beat.segment_velocity_signal,
                schema.artery_per_beat.segment_velocity_signal_band_limited,
                schema.vein_per_beat.segment_velocity_signal,
                schema.vein_per_beat.segment_velocity_signal_band_limited,
                schema.artery_per_beat_safe.velocity_signal,
                schema.artery_per_beat_safe.velocity_signal_band_limited,
                schema.vein_per_beat_safe.velocity_signal,
                schema.vein_per_beat_safe.velocity_signal_band_limited,
            }
            metrics.update(
                {
                    key: value
                    for key, value in velocity_outputs.items()
                    if key not in segment_paths
                }
            )

    if profiles_selected:
        cycle_boundaries = (
            per_beat_result.cycle_boundary_indexes
            if per_beat_result is not None
            else context.per_beat_analysis.cycle_boundary_indexes
        )
        index_base = (
            0
            if per_beat_result is not None
            else int(context.source_data.provenance["beat_index_base"])
        )
        gradient_artery_segments, gradient_vein_segments = (
            extract_spatial_gradient_segments(ctx, context)
        )
        spatial_gradient_outputs = pack_spatial_gradient_profile_outputs(
            gradient_artery_segments,
            gradient_vein_segments,
            cycle_boundaries,
            index_base=index_base,
        )
        metrics.update(spatial_gradient_outputs)

    if profile_products_required:
        velocity_profile_outputs = pack_cross_section_profile_outputs(
            context.artery_segment_result,
            context.vein_segment_result,
            cycle_boundaries,
            index_base=index_base,
        )
        metrics.update(velocity_profile_outputs)
        _validate_profile_segment_alignment(
            "Artery",
            context.artery_segment_result,
            gradient_artery_segments,
        )
        _validate_profile_segment_alignment(
            "Vein",
            context.vein_segment_result,
            gradient_vein_segments,
        )
        metrics.update(
            pack_blood_volume_rate_outputs(
                velocity_profile_outputs,
                spatial_gradient_outputs,
            )
        )
        metrics.update(
            pack_velocity_profile_fft_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
            )
        )
        # Displacement profile metrics are temporarily disabled.
        # metrics.update(
        #     pack_displacement_profile_outputs(
        #         context.artery_segment_result,
        #         context.vein_segment_result,
        #         cycle_boundaries,
        #         index_base=index_base,
        #     )
        # )
        metrics.update(
            pack_displacement_magnitude_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                cycle_boundaries,
                index_base=index_base,
            )
        )
        metrics.update(
            pack_cross_section_displacement_profile_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                cycle_boundaries,
                index_base=index_base,
            )
        )

    if "quadrants" in selected:
        metrics.update(
            pack_quadrant_velocity_outputs(
                velocity_outputs,
                context.source_data,
                context.artery_segment_result,
                context.vein_segment_result,
            )
        )

    cleanup_spatial_gradient_artifacts(ctx)
    return metrics


def _required_state(ctx, key: str):
    value = ctx.state.get(key)
    if value is None:
        raise RuntimeError(
            f"Required pipeline state '{key}' is unavailable; "
            "check the pipeline DAG dependencies."
        )
    return value


def _validate_profile_segment_alignment(
    vessel_name: str,
    velocity_segments,
    gradient_segments,
) -> None:
    """Ensure velocity and gradient profiles share branch/radius identities."""

    if velocity_segments is None or gradient_segments is None:
        raise RuntimeError(
            f"{vessel_name} velocity and gradient segments are required for alignment."
        )

    for field in ("labels", "branch_ids"):
        velocity_value = np.asarray(getattr(velocity_segments, field))
        gradient_value = np.asarray(getattr(gradient_segments, field))
        if not np.array_equal(velocity_value, gradient_value):
            raise RuntimeError(
                f"{vessel_name} velocity and gradient segment {field} do not match."
            )

    velocity_centers = np.asarray(velocity_segments.segment_center_xy)
    gradient_centers = np.asarray(gradient_segments.segment_center_xy)
    if velocity_centers.shape != gradient_centers.shape or not np.allclose(
        velocity_centers,
        gradient_centers,
        equal_nan=True,
    ):
        raise RuntimeError(
            f"{vessel_name} velocity and gradient segment centers do not match."
        )

    velocity_profile_shape = tuple(velocity_segments.velocity_profiles.shape[:2])
    gradient_profile_shape = tuple(gradient_segments.velocity_profiles.shape[:2])
    if velocity_profile_shape != gradient_profile_shape:
        raise RuntimeError(
            f"{vessel_name} velocity and gradient (radius, branch) dimensions "
            f"do not match: {velocity_profile_shape} != {gradient_profile_shape}."
        )
