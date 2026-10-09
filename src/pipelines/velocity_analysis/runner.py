"""Orchestrate selectable velocity output products."""

from time import perf_counter

from input_output import EyeFlowOutputPaths
from utils.logger import Logger

from .analysis.segment_maps import prepare_segment_velocity_maps_per_beat
from .artifacts import export_segment_velocity_map_avis, export_velocity_signals
from .builder import (
    VELOCITY_ANALYSIS_STATE,
    build_velocity_analysis,
)
from .outputs import (
    pack_continuous_velocity_outputs,
    pack_cross_section_profile_outputs,
    pack_quadrant_velocity_outputs,
    pack_segment_map_outputs,
    pack_segment_velocity_outputs,
    pack_velocity_per_beat_outputs,
    pack_velocity_profile_analysis_outputs,
    pack_velocity_profile_fft_outputs,
)


def run_velocity_analysis(ctx) -> dict[str, object]:
    """Publish base velocity plus the selected derived velocity products."""
    analysis = ctx.state.get(VELOCITY_ANALYSIS_STATE)
    if analysis is None:
        analysis = build_velocity_analysis(ctx)
    selected = ctx.options_for("velocity_analysis")
    velocity = analysis.velocity
    velocity_semantics_kwargs = {"velocity": velocity}
    metrics = pack_continuous_velocity_outputs(velocity)
    segments_available = (
        analysis.artery_segments is not None or analysis.vein_segments is not None
    )
    maps_selected = "segment_velocity_maps" in selected
    profiles_selected = bool(
        {
            "velocity_profiles",
            "velocity_profile_analysis",
            "velocity_profile_fft",
        }
        & selected
    )
    profile_analysis_selected = "velocity_profile_analysis" in selected
    profile_fft_selected = "velocity_profile_fft" in selected
    artery_velocity_maps_per_beat = None
    vein_velocity_maps_per_beat = None
    if maps_selected:
        map_started = perf_counter()
        Logger.log("Starting shared per-beat segment velocity-map interpolation...")
        artery_velocity_maps_per_beat, vein_velocity_maps_per_beat = (
            prepare_segment_velocity_maps_per_beat(
                analysis.artery_segments,
                analysis.vein_segments,
                analysis.cycle_boundary_indexes,
                index_base=0,
            )
        )
        Logger.log(
            "Completed shared per-beat segment velocity-map interpolation in "
            f"{perf_counter() - map_started:.1f}s."
        )
    if segments_available:
        metrics.update(
            pack_segment_velocity_outputs(
                analysis.artery_segments,
                analysis.vein_segments,
                source_data=analysis.source_data,
                **velocity_semantics_kwargs,
            )
        )
    if maps_selected:
        segment_map_outputs = pack_segment_map_outputs(
            analysis.artery_segments,
            analysis.vein_segments,
            artery_velocity_maps_per_beat,
            vein_velocity_maps_per_beat,
            **velocity_semantics_kwargs,
        )
        metrics.update(segment_map_outputs)
        output = getattr(ctx, "output", None)
        if getattr(output, "available", False):
            avi_started = perf_counter()
            Logger.log("Starting segment velocity-map AVI export...")
            export_segment_velocity_map_avis(
                output,
                analysis.artery_segments,
                analysis.vein_segments,
                segment_map_outputs,
            )
            Logger.log(
                "Completed segment velocity-map AVI export in "
                f"{perf_counter() - avi_started:.1f}s."
            )

    velocity_outputs = pack_velocity_per_beat_outputs(
        analysis.per_beat_result,
        **velocity_semantics_kwargs,
    )
    metrics.update(velocity_outputs)

    output = getattr(ctx, "output", None)
    if getattr(output, "available", False):
        schema = EyeFlowOutputPaths.active()
        artery_path = schema.artery_per_beat_safe.velocity_signal
        vein_path = schema.vein_per_beat_safe.velocity_signal
        if artery_path in velocity_outputs and vein_path in velocity_outputs:
            export_velocity_signals(
                output,
                velocity_outputs[artery_path],
                velocity_outputs[vein_path],
            )

    if profiles_selected:
        velocity_profile_outputs = pack_cross_section_profile_outputs(
            analysis.artery_segments,
            analysis.vein_segments,
            analysis.cycle_boundary_indexes,
            index_base=0,
            **velocity_semantics_kwargs,
        )
        metrics.update(velocity_profile_outputs)
        if profile_analysis_selected:
            metrics.update(
                pack_velocity_profile_analysis_outputs(velocity_profile_outputs)
            )
        if profile_fft_selected:
            metrics.update(
                pack_velocity_profile_fft_outputs(
                    analysis.artery_segments,
                    analysis.vein_segments,
                )
            )
    if "quadrants" in selected:
        metrics.update(
            pack_quadrant_velocity_outputs(
                velocity_outputs,
                analysis.source_data,
                analysis.artery_segments,
                analysis.vein_segments,
                **velocity_semantics_kwargs,
            )
        )

    return metrics
