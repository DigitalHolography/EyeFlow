"""Orchestrate selectable velocity output products."""

from time import perf_counter

from input_output import EyeFlowOutputPaths
from utils.logger import Logger

from .continuous import (
    pack_continuous_velocity_outputs,
    pack_segment_velocity_outputs,
)
from .outputs import export_velocity_signals
from .per_beat_outputs import pack_velocity_per_beat_outputs
from .profiles import (
    pack_cross_section_profile_outputs,
    pack_velocity_profile_fft_outputs,
)
from .quadrants import pack_quadrant_velocity_outputs
from .segment_maps import (
    pack_segment_map_outputs,
    prepare_segment_velocity_maps_per_beat,
)
from .segment_velocity_map_avi import export_segment_velocity_map_avis
from .workflow import (
    WAVEFORM_VELOCITY_STATE,
    build_waveform_velocity,
)


def run_waveform_velocity(ctx) -> dict[str, object]:
    """Publish base velocity plus the selected derived velocity products."""
    waveform = ctx.state.get(WAVEFORM_VELOCITY_STATE)
    if waveform is None:
        waveform = build_waveform_velocity(ctx)
    selected = ctx.options_for("waveform_velocity")
    metrics = pack_continuous_velocity_outputs(waveform.retinal_velocity)
    segments_selected = "segments" in selected
    maps_selected = "segment_velocity_maps" in selected
    profiles_selected = bool(
        {"velocity_profiles", "velocity_profile_fft"} & selected
    )
    profile_fft_selected = "velocity_profile_fft" in selected
    profile_analysis_scheduled = ctx.pipeline_scheduled(
        "velocity_profile_analysis"
    )
    artery_velocity_maps_per_beat = None
    vein_velocity_maps_per_beat = None
    if maps_selected:
        map_started = perf_counter()
        Logger.log("Starting shared per-beat segment velocity-map interpolation...")
        artery_velocity_maps_per_beat, vein_velocity_maps_per_beat = (
            prepare_segment_velocity_maps_per_beat(
                waveform.artery_segments,
                waveform.vein_segments,
                waveform.cycle_boundary_indexes,
                index_base=0,
            )
        )
        Logger.log(
            "Completed shared per-beat segment velocity-map interpolation in "
            f"{perf_counter() - map_started:.1f}s."
        )
    if segments_selected:
        metrics.update(
            pack_segment_velocity_outputs(
                waveform.artery_segments,
                waveform.vein_segments,
                source_data=waveform.source_data,
            )
        )
    if maps_selected:
        segment_map_outputs = pack_segment_map_outputs(
            waveform.artery_segments,
            waveform.vein_segments,
            artery_velocity_maps_per_beat,
            vein_velocity_maps_per_beat,
        )
        metrics.update(segment_map_outputs)
        output = getattr(ctx, "output", None)
        if getattr(output, "available", False):
            avi_started = perf_counter()
            Logger.log("Starting segment velocity-map AVI export...")
            export_segment_velocity_map_avis(
                output,
                waveform.artery_segments,
                waveform.vein_segments,
                segment_map_outputs,
            )
            Logger.log(
                "Completed segment velocity-map AVI export in "
                f"{perf_counter() - avi_started:.1f}s."
            )

    per_beat_result = waveform.per_beat_result
    velocity_outputs = (
        pack_velocity_per_beat_outputs(per_beat_result)
        if per_beat_result is not None
        else {}
    )
    if "per_beat" in selected or ctx.pipeline_scheduled("pdf_report"):
        if per_beat_result is None:
            raise RuntimeError(
                "Per-beat waveform products were selected but were not computed."
            )
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

    profile_products_required = profiles_selected or profile_analysis_scheduled
    if profile_products_required:
        velocity_profile_outputs = pack_cross_section_profile_outputs(
            waveform.artery_segments,
            waveform.vein_segments,
            waveform.cycle_boundary_indexes,
            index_base=0,
        )
        metrics.update(velocity_profile_outputs)
        if profile_fft_selected:
            metrics.update(
                pack_velocity_profile_fft_outputs(
                    waveform.artery_segments,
                    waveform.vein_segments,
                )
            )
    if "quadrants" in selected:
        metrics.update(
            pack_quadrant_velocity_outputs(
                velocity_outputs,
                waveform.source_data,
                waveform.artery_segments,
                waveform.vein_segments,
            )
        )

    return metrics
