"""Executable velocity analysis and its shared scientific input adapters."""

from contextlib import contextmanager
from dataclasses import replace
from functools import partial
from time import perf_counter

import numpy as np

from calculations.blood_flow_velocity import (
    PerBeatAnalysisInput,
    PerBeatAnalysisResult,
    run_per_beat_analysis,
)
from calculations.topology import AnnulusGeometry
from data_manager import InvocationResource
from input_output import EyeFlowOutputPaths
from pipeline_engine import ProductSpec, StageGraph, StageSpec, read_int_setting
from pipeline_engine.context import require_managed_runtime
from pipelines.signal_settings import band_limited_harmonics
from pipelines.topology_core.runner import prepared_topologies
from pipelines.velocity.models import RetinalVelocity
from pipelines.velocity.runner import velocity
from pipelines.velocity.signal_processing import DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ
from utils.logger import Logger

from .analysis.filtering import lowpass_velocity_signals
from .analysis.segment_maps import prepare_segment_velocity_maps_per_beat
from .analysis.segments import analyze_velocity_segment_profiles
from .artifacts import (
    export_segment_velocity_map_avis,
    export_velocity_signals,
    export_branch_identity_stage_pngs,
    export_pulse_pngs,
    export_rotated_mean_pngs,
)
from .constants import NUMBER_OF_RADII_IN_FOV, SEGMENT_INNER_RADIUS_FRAC, SEGMENT_OUTER_RADIUS_FRAC
from .models import VelocityAnalysis, VelocitySegmentResult
from .outputs import (
    pack_velocity_per_beat_inputs,
    pack_continuous_velocity_outputs,
    pack_cross_section_profile_outputs,
    pack_quadrant_velocity_outputs,
    pack_segment_map_outputs,
    pack_segment_velocity_outputs,
    pack_velocity_profile_analysis_outputs,
    pack_velocity_profile_fft_outputs,
)
from .sources import VelocityAnalysisSourceData, load_velocity_analysis_source_data


VELOCITY_ANALYSIS_STATE = "velocity_analysis"


def velocity_analysis_graph(*, options=(), selection=None, segments_required=None, pulse_figures=True):
    if selection is not None:
        segments_required = _segments_required(selection)
        pulse_figures = _pulse_pngs_required(selection)
    options = set(options)
    has_segments = bool(options) if segments_required is None else segments_required
    packed = lambda value, ctx: value
    stages = [StageSpec("inputs", _inputs, products=(ProductSpec("inputs"),))]
    if has_segments:
        stages.append(StageSpec("segments", _segments, requires=("inputs",), products=(ProductSpec("segments"),)))
    stages.extend((
        StageSpec("per_beat", partial(_per_beat, has_segments), requires=("inputs",) + (("segments",) if has_segments else ()),
                  products=(ProductSpec("analysis"),)),
        StageSpec("continuous", _continuous, requires=("analysis",), products=(ProductSpec("continuous", pack=packed),)),
    ))
    if has_segments:
        stages.append(StageSpec("segment_signals", _segment_signals, requires=("analysis",),
                                products=(ProductSpec("segment_signals", pack=packed),)))
    if "segment_velocity_maps" in options:
        stages.append(StageSpec("maps", _maps, requires=("analysis",), products=(ProductSpec("maps", pack=packed),)))
    stages.append(StageSpec("waveform_inputs", _per_beat_outputs, requires=("analysis",),
                           products=(ProductSpec("waveform_inputs", pack=_pack_per_beat),)))
    if {"velocity_profiles", "velocity_profile_analysis", "velocity_profile_fft"} & options:
        stages.append(StageSpec("profiles", _profiles, requires=("analysis",), products=(ProductSpec("profiles", pack=packed),)))
    if "velocity_profile_analysis" in options:
        stages.append(StageSpec("profile_fits", _profile_fits, requires=("profiles",), products=(ProductSpec("profile_fits", pack=packed),)))
    if "velocity_profile_fft" in options:
        stages.append(StageSpec("profile_fft", _profile_fft, requires=("analysis",), products=(ProductSpec("profile_fft", pack=packed),)))
    if "quadrants" in options:
        stages.append(StageSpec("quadrants", _quadrants, requires=("analysis", "waveform_inputs"), products=(ProductSpec("quadrants", pack=packed),)))
    stages.append(StageSpec("publish", _publish_analysis, requires=("analysis", "waveform_inputs")))
    if has_segments:
        stages.append(StageSpec("segment_figures", _segment_figures, requires=("inputs", "segments"), artifact=True))
        stages.append(StageSpec("signal_figures", _signal_figures, requires=("waveform_inputs",), artifact=True))
    if pulse_figures:
        stages.append(StageSpec("pulse_figures", _pulse_figures, requires=("analysis", "waveform_inputs", "continuous"), artifact=True))
    if "segment_velocity_maps" in options:
        stages.append(StageSpec("map_videos", _map_videos, requires=("analysis", "maps"), artifact=True))
    return StageGraph(tuple(stages))


def run_velocity_analysis(ctx):
    require_managed_runtime(ctx, "velocity_analysis")
    return velocity_analysis_graph(options=ctx.options_for("velocity_analysis"),
                                   selection=ctx)


def _inputs(ctx):
    ctx.require_inputs("hd", "dv")
    retinal = velocity(ctx)
    source = load_velocity_analysis_source_data(ctx)
    harmonic_count = band_limited_harmonics(ctx)
    attrs = _context_attrs(source, retinal, harmonic_count, "eyeflow_velocity_analysis",
                           _number_of_radii_in_fov(ctx), retinal.cardiac_cycle_source)
    ctx.dag.handoff("inputs", (retinal, InvocationResource(source), harmonic_count, attrs))


def _segments(ctx):
    retinal, source, _, _ = ctx.dag.require("inputs")
    artery, vein = _build_segments(ctx, retinal, source, segments_required=True)
    dt = float(source.source.holodoppler.timing.dt_seconds)
    segments = tuple(replace(value, filtered_signal=lowpass_velocity_signals(
        value.profile.segment_signal, dt_seconds=dt)) if value is not None else None
        for value in (artery, vein))
    ctx.dag.handoff("segments", segments)


def _per_beat(has_segments, ctx):
    retinal, source, harmonics, attrs = ctx.dag.require("inputs")
    artery, vein = ctx.dag.require("segments") if has_segments else (None, None)
    with _logged_stage("shared per-beat velocity analysis"):
        result = _run_velocity_per_beat_analysis(retinal, source, artery, vein, harmonic_count=harmonics)
    analysis = VelocityAnalysis(retinal, InvocationResource(source), artery, vein, result, attrs)
    accepted = ctx.dag.handoff("analysis", analysis)
    ctx.state.set(VELOCITY_ANALYSIS_STATE, accepted)


def _continuous(ctx):
    ctx.dag.handoff("continuous", pack_continuous_velocity_outputs(ctx.dag.require("analysis").velocity))


def _segment_signals(ctx):
    analysis = ctx.dag.require("analysis")
    ctx.dag.handoff("segment_signals", pack_segment_velocity_outputs(analysis.artery_segments,
        analysis.vein_segments, source_data=analysis.source_data, velocity=analysis.velocity))


def _per_beat_outputs(ctx):
    analysis = ctx.dag.require("analysis")
    values = pack_velocity_per_beat_inputs(analysis.per_beat_result, velocity=analysis.velocity)
    ctx.dag.handoff("waveform_inputs", values)


def _publish_analysis(ctx):
    analysis = ctx.dag.require("analysis")
    # All downstream waveform families reuse this accepted representation.
    ctx.state.set(VELOCITY_ANALYSIS_STATE, replace(analysis, waveform_inputs=ctx.dag.require("waveform_inputs")))


def _pack_per_beat(values, ctx):
    duration = EyeFlowOutputPaths.active().cardiac_cycle.systolic_cycle_duration_seconds
    return {path: value for path, value in values.items() if path != duration}


def _maps(ctx):
    from input_output.writers.h5_storage import velocity_map_h5_options
    analysis = ctx.dag.require("analysis")
    with ctx.output.h5.build_arrays() as allocate:
        artery, vein = prepare_segment_velocity_maps_per_beat(analysis.artery_segments,
            analysis.vein_segments, analysis.cycle_boundary_indexes, index_base=0,
            array_allocator=lambda shape: allocate(shape, h5_options=velocity_map_h5_options(shape)))
    ctx.dag.handoff("maps", pack_segment_map_outputs(analysis.artery_segments,
        analysis.vein_segments, artery, vein, velocity=analysis.velocity))


def _profiles(ctx):
    analysis = ctx.dag.require("analysis")
    ctx.dag.handoff("profiles", pack_cross_section_profile_outputs(analysis.artery_segments,
        analysis.vein_segments, analysis.cycle_boundary_indexes, index_base=0, velocity=analysis.velocity))


def _profile_fits(ctx):
    ctx.dag.handoff("profile_fits", pack_velocity_profile_analysis_outputs(ctx.dag.require("profiles")))


def _profile_fft(ctx):
    analysis = ctx.dag.require("analysis")
    ctx.dag.handoff("profile_fft", pack_velocity_profile_fft_outputs(analysis.artery_segments, analysis.vein_segments))


def _quadrants(ctx):
    analysis = ctx.dag.require("analysis")
    ctx.dag.handoff("quadrants", pack_quadrant_velocity_outputs(ctx.dag.require("waveform_inputs"),
        analysis.source_data, analysis.artery_segments, analysis.vein_segments, velocity=analysis.velocity))


def _segment_figures(ctx):
    if not ctx.output.available:
        return
    _, source, _, _ = ctx.dag.require("inputs")
    for name, segments in zip(("artery", "vein"), ctx.dag.require("segments")):
        if segments is None:
            continue
        export_rotated_mean_pngs(ctx.output, segments, "arteries" if name == "artery" else "veins")
        _export_branch_identity_debug(ctx, segments, source.source.segmentation.optic_disc.center,
            segments.profile.topology.native.ring_settings, name)


def _pulse_figures(ctx):
    analysis = ctx.dag.require("analysis")
    with _optional_artifact_export("pulse-analysis PNG"):
        _export_pulse_pngs(ctx, analysis, analysis.per_beat_result)


def _map_videos(ctx):
    analysis = ctx.dag.require("analysis")
    with _optional_artifact_export("segment velocity-map AVI"):
        export_segment_velocity_map_avis(ctx.output, analysis.artery_segments,
            analysis.vein_segments, ctx.dag.require("maps"))


def _signal_figures(ctx):
    values = ctx.dag.require("waveform_inputs")
    paths = EyeFlowOutputPaths.active()
    artery = paths.artery_per_beat_filtered.segment_velocity_signal_band_limited
    vein = paths.vein_per_beat_filtered.segment_velocity_signal_band_limited
    if artery in values and vein in values:
        with _optional_artifact_export("velocity-signal figure"):
            export_velocity_signals(ctx.output, values[artery], values[vein])


def velocity_analysis(ctx) -> VelocityAnalysis:
    """Return the typed result produced by the velocity-analysis pipeline."""

    value = ctx.state.get(VELOCITY_ANALYSIS_STATE)
    if not isinstance(value, VelocityAnalysis):
        raise RuntimeError(
            "Velocity analysis state is unavailable; check the pipeline DAG dependency."
        )
    return value


@contextmanager
def _optional_artifact_export(label):
    try:
        yield
    except Exception as exc:
        Logger.log_warning(f"Failed {label} export; calculated HDF5 results retained: {exc}")


def _segments_required(ctx) -> bool:
    """Return whether the canonical segment analysis must run."""
    if ctx.pipeline_scheduled("blood_volume_rate") and ctx.options_for("blood_volume_rate"):
        return True
    if ctx.pipeline_scheduled("lowrank_waveform_decomposition"):
        return True
    if ctx.pipeline_scheduled("velocity_analysis") and {
        "segments",
        "segment_velocity_maps",
        "velocity_profiles",
        "velocity_profile_analysis",
        "velocity_profile_fft",
        "quadrants",
    } & ctx.options_for("velocity_analysis"):
        return True
    if ctx.pipeline_scheduled("waveform_shape_metrics") and {
        "segments",
        "quadrants",
    } & ctx.options_for("waveform_shape_metrics"):
        return True
    return bool(
        ctx.pipeline_scheduled("absolute_waveform_metrics")
        and {"segments", "quadrants"} & ctx.options_for("absolute_waveform_metrics")
    )


def _pulse_pngs_required(ctx) -> bool:
    return ctx.pipeline_targeted("velocity_analysis") or ctx.pipeline_scheduled("pdf_report")


def _number_of_radii_in_fov(ctx) -> int:
    value = read_int_setting(
        ctx,
        default=NUMBER_OF_RADII_IN_FOV,
        keys=(
            "number_of_radii_in_FOV",
            "number_of_radii_in_fov",
            "NumberOfRadiiInFOV",
            # Accept the earlier spelling for compatibility.
            "number_of_radii_over_FOV",
            "number_of_radii_over_fov",
            "NumberOfRadiiOverFOV",
        ),
    )
    if value < 1:
        raise ValueError("number_of_radii_in_FOV must be positive.")
    return value


def _build_segments(
    ctx,
    retinal: RetinalVelocity,
    source_data: VelocityAnalysisSourceData,
    *,
    segments_required: bool,
) -> tuple[
    VelocitySegmentResult | None,
    VelocitySegmentResult | None,
]:
    if not segments_required:
        Logger.log("Skipping segment velocity extraction; no selected output requires it.")
        return None, None

    velocity_map = retinal.maps.velocity
    if velocity_map is None:
        raise ValueError("velocity_map is required for segment extraction.")
    shared_topologies = prepared_topologies(ctx)
    ring_settings = shared_topologies["artery"].native.ring_settings
    if not isinstance(ring_settings, AnnulusGeometry):
        raise RuntimeError("Prepared topology has no annulus geometry.")
    return _segment_velocity_inputs(
        velocity_map,
        source_data,
        ring_settings,
        ctx,
        cycle_boundary_indexes=retinal.cycle_boundary_indexes,
        prepared_topologies=shared_topologies,
    )


def _run_velocity_per_beat_analysis(
    retinal: RetinalVelocity,
    source_data: VelocityAnalysisSourceData,
    artery_segments: VelocitySegmentResult | None,
    vein_segments: VelocitySegmentResult | None,
    *,
    harmonic_count: int,
) -> PerBeatAnalysisResult:
    timing = source_data.source.holodoppler.timing
    arterial_velocity_signal, venous_velocity_signal = _raw_velocity_signals_for_per_beat(retinal)
    artery_segment_signal = _waveform_segment_input(artery_segments)
    vein_segment_signal = _waveform_segment_input(vein_segments)
    return run_per_beat_analysis(
        PerBeatAnalysisInput(
            arterial_velocity_signal=arterial_velocity_signal,
            venous_velocity_signal=venous_velocity_signal,
            cycle_boundary_indexes=np.asarray(
                retinal.cycle_boundary_indexes,
                dtype=np.int32,
            ),
            band_limited_harmonics=harmonic_count,
            cardiac_cycle=retinal.cardiac_cycle.spectral,
            dt_seconds=timing.dt_seconds,
            arterial_velocity_segments=artery_segment_signal,
            venous_velocity_segments=vein_segment_signal,
            arterial_filtered_velocity_signal=retinal.continuous("artery"),
            arterial_filtered_velocity_segments=_filtered_segment_input(
                artery_segments, timing.dt_seconds
            ),
            venous_filtered_velocity_signal=retinal.continuous("vein"),
            venous_filtered_velocity_segments=_filtered_segment_input(
                vein_segments, timing.dt_seconds
            ),
            index_base=source_data.provenance["beat_index_base"],
        )
    )


def _filtered_segment_input(segments, dt_seconds):
    if segments is None:
        return None
    cached = getattr(segments, "filtered_signal", None)
    if cached is None:
        cached = lowpass_velocity_signals(
            segments.profile.segment_signal, dt_seconds=float(dt_seconds)
        )
    return cached


def _raw_velocity_signals_for_per_beat(
    retinal: RetinalVelocity,
) -> tuple[np.ndarray, np.ndarray]:
    return (
        retinal.continuous("artery", raw=True),
        retinal.continuous("vein", raw=True),
    )


def _segment_velocity_inputs(
    velocity_map,
    source_data: VelocityAnalysisSourceData,
    ring_settings: AnnulusGeometry,
    ctx,
    *,
    cycle_boundary_indexes,
    prepared_topologies=None,
) -> tuple[VelocitySegmentResult, VelocitySegmentResult]:
    velocity_analysis_scheduled = ctx.pipeline_scheduled("velocity_analysis")
    retain_velocity_maps = bool(
        velocity_analysis_scheduled
        and ctx.option_enabled(
            "segment_velocity_maps",
            pipeline="velocity_analysis",
        )
    )
    velocity_profile_fft = bool(
        velocity_analysis_scheduled
        and ctx.option_enabled(
            "velocity_profile_fft",
            pipeline="velocity_analysis",
        )
    )
    # All three options export both cross-section profile directions. Gradient-
    # edge blood-volume rate additionally consumes transverse velocity profiles;
    # other segment consumers only need signals, means, or maps.
    compute_profile_outputs = bool(
        velocity_analysis_scheduled
        and any(
            ctx.option_enabled(option, pipeline="velocity_analysis")
            for option in (
                "velocity_profiles",
                "velocity_profile_analysis",
                "velocity_profile_fft",
            )
        )
    )
    compute_transverse_profiles = bool(
        compute_profile_outputs
        or (
            ctx.pipeline_scheduled("blood_volume_rate")
            and ctx.option_enabled("gradient_edges", pipeline="blood_volume_rate")
        )
    )
    compute_longitudinal_profiles = compute_profile_outputs
    with _logged_stage("segment velocity extraction"), ctx.output.h5.build_arrays() as allocate:
        results = analyze_velocity_segment_profiles(
            velocity_map,
            {
                "artery": source_data.source.segmentation.vessels.artery,
                "vein": source_data.source.segmentation.vessels.vein,
            },
            source_data.source.segmentation.optic_disc,
            ring_settings,
            source_data.profile_settings,
            prepared_topologies=prepared_topologies,
            retain_velocity_maps=retain_velocity_maps,
            compute_transverse_profiles=compute_transverse_profiles,
            compute_longitudinal_profiles=compute_longitudinal_profiles,
            cycle_boundary_indexes=cycle_boundary_indexes,
            velocity_profile_fft=velocity_profile_fft,
            index_base=int(source_data.provenance["beat_index_base"]),
            map_allocator=allocate if retain_velocity_maps else None,
        )
    return results["artery"], results["vein"]


def _waveform_segment_input(
    result: VelocitySegmentResult | None,
) -> np.ndarray | None:
    if result is None:
        return None
    return result.profile.segment_signal


def _export_branch_identity_debug(
    ctx,
    result: VelocitySegmentResult,
    optic_disc_center,
    ring_settings: AnnulusGeometry,
    prefix: str,
) -> None:
    if not ctx.output.available:
        return
    export_branch_identity_stage_pngs(
        ctx.output,
        result.profile.topology.native.branches.require_stages(),
        prefix,
        optic_disc_center,
        ring_settings,
        segment_centers_xy=result.profile.topology.native.segment_centers_xy,
        profile_window_bounds_xyxy=result.profile.topology.native.window_bounds_xyxy,
    )


def _export_pulse_pngs(ctx, analysis: VelocityAnalysis, per_beat_result) -> None:
    if not ctx.output.available:
        return
    with _logged_stage("pulse-analysis PNG export"):
        from input_output.writers.png import REPORT_IMAGES_STATE

        report_images = {}
        export_pulse_pngs(
            ctx.output, analysis, per_beat_result, report_images=report_images
        )
        ctx.state.set(REPORT_IMAGES_STATE, report_images)


@contextmanager
def _logged_stage(label: str):
    started = perf_counter()
    Logger.log_perf(f"START pipeline=velocity_analysis stage={label} at={started:.6f}")
    yield
    Logger.log_perf(
        f"END pipeline=velocity_analysis stage={label} at={perf_counter():.6f} "
        f"duration={perf_counter() - started:.3f}s"
    )


def _context_attrs(
    source_data: VelocityAnalysisSourceData,
    velocity: RetinalVelocity,
    harmonic_count: int,
    analysis_source: str,
    number_of_radii_in_fov: int,
    cardiac_cycle_detection_source: str,
) -> dict[str, object]:
    output_paths = EyeFlowOutputPaths.active()
    analysis_paths = output_paths.analysis
    velocity_input_dependency = (
        "holodoppler.h5.frequency_bands"
        if source_data.source.velocity_estimation_method == "frequency_bands"
        else "holodoppler.h5.moment0_moment2"
    )
    dependency_chain = (
        ["external_velocity_analysis"]
        if analysis_source == "external_velocity_analysis"
        else [
            velocity_input_dependency,
            "retinal_segmentation_input",
            "eyeflow.velocity.recomputed",
        ]
    )
    attrs = {
        "dependency_chain": dependency_chain
        + [
            "blood_flow_velocity.signal_analysis.cardiac_cycle.spectral",
            "blood_flow_velocity.signal_analysis.per_beat.signal",
            "blood_flow_velocity.signal_analysis.per_beat.runner",
        ],
        "analysis_source": analysis_source,
        "output_schema": output_paths.name,
        "velocity_section_geometry": "optic_disc_centered_frame_fraction",
        "velocity_section_inner_radius_fraction": float(SEGMENT_INNER_RADIUS_FRAC),
        "velocity_section_outer_radius_fraction": float(SEGMENT_OUTER_RADIUS_FRAC),
        "velocity_section_outer_to_disc_radius": float(
            SEGMENT_OUTER_RADIUS_FRAC / SEGMENT_INNER_RADIUS_FRAC
        ),
        "number_of_radii_in_FOV": int(number_of_radii_in_fov),
        "arterial_velocity_signal_path": (analysis_paths.retinal_artery_velocity_signal),
        "venous_velocity_signal_path": (analysis_paths.retinal_vein_velocity_signal),
        "systolic_peak_indexes_path": analysis_paths.beat_indices,
        "cycle_duration_seconds_path": (output_paths.cardiac_cycle.systolic_cycle_duration_seconds),
        "band_limited_harmonics": int(harmonic_count),
        "filter_velocity_signals": True,
        "velocity_signal_lowpass_hz": float(DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ),
        "cardiac_cycle_detection_source": cardiac_cycle_detection_source,
    }
    for key, value in source_data.provenance.items():
        if value is not None:
            attrs[key] = value
    for key, value in velocity.provenance.items():
        if value is not None:
            attrs[key] = value
    return attrs
