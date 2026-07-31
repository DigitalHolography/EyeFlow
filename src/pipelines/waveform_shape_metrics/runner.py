"""Orchestrate the waveform-shape metrics sandbox pipeline."""

from collections.abc import Mapping
from dataclasses import dataclass

from calculations.blood_flow_velocity import (
    CrossSectionSignalResult,
    HeartbeatAnalysisResult,
    PerBeatAnalysisInput,
    SegmentRingSettings,
    segment_velocity_results,
    spectral_heartbeat_analysis,
)
from input_output import EyeFlowOutputPaths
from pipeline_engine.imports import (
    HolodopplerTiming,
    np,
    read_int_setting,
)

from .dopplerview.constants import (
    LEGACY_FILTER_VELOCITY_SIGNALS,
    LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
)
from .velocity.constants import (
    LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
    SEGMENT_RING_COUNT,
    SEGMENT_INNER_RADIUS_FRAC,
    SEGMENT_LENGTH_FRAC,
    SEGMENT_OUTER_RADIUS_FRAC,
)
from .dopplerview.outputs import pack_dopplerview_analysis_outputs
from .dopplerview.runner import run_dopplerview_analysis
from .outputs import pack_waveform_shape_outputs
from .scratch import waveform_scratch_h5
from .sources import WaveformShapeSourceData, WaveformShapeSources
from .velocity.branch_identity_debug import export_branch_identity_stage_pngs
from .velocity.figures import export_pulse_pngs
from .velocity.profiles import pack_cross_section_profile_outputs
from .velocity.runner import run_velocity_per_beat_metrics


@dataclass(frozen=True)
class WaveformShapeMetricsContext:
    source_data: WaveformShapeSourceData
    per_beat_analysis: PerBeatAnalysisInput
    artery_segment_result: CrossSectionSignalResult
    vein_segment_result: CrossSectionSignalResult
    dopplerview_analysis: dict[str, object]
    attrs: dict[str, object]


def run_waveform_shape_metrics(ctx) -> tuple[dict[str, object], dict[str, object]]:
    ctx.require_inputs("hd", "dv")

    with waveform_scratch_h5(ctx) as scratch_h5:
        _log(ctx, "Starting waveform-shape metrics context build...")
        context = _build_waveform_shape_metrics_context(ctx, scratch_h5)

        _log(ctx, "Starting per-beat analysis...")
        per_beat_result, velocity_metrics = run_velocity_per_beat_metrics(context)
        metrics = pack_dopplerview_analysis_outputs(context.dopplerview_analysis)
        metrics.update(velocity_metrics)
        metrics.update(_pack_meta_outputs(context))
        metrics.update(
            pack_cross_section_profile_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                per_beat_result.cycle_boundary_indexes,
                index_base=0,
            )
        )
        ctx.state.set("velocity_per_beat_result", per_beat_result)
        _export_pulse_pngs(ctx, context, per_beat_result)

        _log(ctx, "Starting waveform-shape metric calculation...")
        waveform_metric_outputs = pack_waveform_shape_outputs(
            metrics,
            context.source_data,
            context.artery_segment_result,
            context.vein_segment_result,
        )
        metrics.update(waveform_metric_outputs)
        ctx.state.set("waveform_shape_metric_outputs", waveform_metric_outputs)
        attrs = context.attrs

    return metrics, attrs


def _build_waveform_shape_metrics_context(ctx, scratch_h5) -> WaveformShapeMetricsContext:
    _log(ctx, "Starting waveform source loading...")
    source_data = WaveformShapeSources.from_context(ctx).load()
    timing = source_data.timing
    dopplerview_analysis, analysis_source = _resolve_dopplerview_analysis(
        source_data,
        ctx,
        scratch_h5,
    )
    harmonic_count = _band_limited_harmonic_count(ctx)
    per_beat_analysis, artery_segments, vein_segments = _per_beat_input_from_analysis(
        dopplerview_analysis,
        source_data,
        timing,
        harmonic_count,
        ctx,
    )

    return WaveformShapeMetricsContext(
        source_data=source_data,
        per_beat_analysis=per_beat_analysis,
        artery_segment_result=artery_segments,
        vein_segment_result=vein_segments,
        dopplerview_analysis=dopplerview_analysis,
        attrs=_context_attrs(
            source_data,
            timing,
            harmonic_count,
            analysis_source,
            per_beat_analysis.heartbeat,
        ),
    )


def _resolve_dopplerview_analysis(
    source_data: WaveformShapeSourceData,
    ctx,
    scratch_h5,
) -> tuple[dict[str, object], str]:
    _log(
        ctx,
        "Building EyeFlow velocity analysis from HD moments and DV segmentation.",
    )
    return (
        run_dopplerview_analysis(
            source_data,
            scratch_h5,
            log=getattr(ctx, "log", None),
        ),
        "eyeflow_recomputed_dopplerview_analysis",
    )


def _band_limited_harmonic_count(ctx) -> int:
    return read_int_setting(
        ctx,
        default=LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
        keys=("BandLimitedSignalHarmonicCount", "band_limited_signal_harmonic_count"),
    )


def _per_beat_input_from_analysis(
    analysis: Mapping[str, object],
    source_data: WaveformShapeSourceData,
    timing: HolodopplerTiming,
    harmonic_count: int,
    ctx,
) -> tuple[
    PerBeatAnalysisInput,
    CrossSectionSignalResult,
    CrossSectionSignalResult,
]:
    ring_settings = _segment_ring_settings()
    artery_segments, vein_segments = _segment_velocity_inputs(
        analysis,
        source_data,
        ring_settings,
        ctx,
    )
    arterial_velocity_signal, venous_velocity_signal = (
        _filtered_velocity_signals_for_per_beat(analysis)
    )
    beat_indexes = np.asarray(
        analysis["beat_indices"],
        dtype=np.int32,
    )
    cached_heartbeat = analysis.get("_heartbeat_analysis_result")
    heartbeat = (
        cached_heartbeat.spectral
        if isinstance(cached_heartbeat, HeartbeatAnalysisResult)
        else spectral_heartbeat_analysis(
            arterial_velocity_signal,
            timing.dt_seconds,
            beat_indexes.size,
        )
    )
    inputs = PerBeatAnalysisInput(
        arterial_velocity_signal=arterial_velocity_signal,
        venous_velocity_signal=venous_velocity_signal,
        cycle_boundary_indexes=beat_indexes,
        band_limited_signal_harmonic_count=harmonic_count,
        heartbeat=heartbeat,
        dt_seconds=timing.dt_seconds,
        arterial_velocity_segments=_waveform_segment_input(artery_segments),
        venous_velocity_segments=_waveform_segment_input(vein_segments),
        index_base=source_data.provenance["beat_index_base"],
    )
    return inputs, artery_segments, vein_segments


def _filtered_velocity_signals_for_per_beat(
    analysis: Mapping[str, object],
) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.asarray(
            analysis["retinal_artery_velocity_signal_filtered"],
            dtype=np.float32,
        ),
        np.asarray(
            analysis["retinal_vein_velocity_signal_filtered"],
            dtype=np.float32,
        ),
    )


def _segment_velocity_inputs(
    analysis: Mapping[str, object],
    source_data: WaveformShapeSourceData,
    ring_settings: SegmentRingSettings,
    ctx,
) -> tuple[CrossSectionSignalResult, CrossSectionSignalResult]:
    _log(ctx, "Starting segment velocity extraction...")
    artery, vein = segment_velocity_results(
        analysis["retinal_vessel_velocity"],
        source_data.retinal_artery_mask,
        source_data.retinal_vein_mask,
        source_data.optic_disc_center,
        ring_settings,
        source_data.cross_section_settings,
        log=getattr(ctx, "log", None),
    )
    _export_branch_identity_debug(
        ctx,
        artery.branch_identity.stages,
        source_data.optic_disc_center,
        ring_settings,
        "artery",
    )
    _export_branch_identity_debug(
        ctx,
        vein.branch_identity.stages,
        source_data.optic_disc_center,
        ring_settings,
        "vein",
    )
    return artery, vein


def _waveform_segment_input(result: CrossSectionSignalResult) -> np.ndarray | None:
    if result.branch_ids.size == 0:
        return None
    return result.velocity


def _log(ctx, message: str) -> None:
    log = getattr(ctx, "log", None)
    if callable(log):
        log(message)


def _export_branch_identity_debug(
    ctx,
    stages,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    prefix: str,
) -> None:
    if not ctx.output.available:
        return
    export_branch_identity_stage_pngs(
        ctx.output,
        stages,
        prefix,
        optic_disc_center,
        ring_settings,
    )


def _export_pulse_pngs(ctx, context: WaveformShapeMetricsContext, per_beat_result) -> None:
    if not ctx.output.available:
        return
    _log(ctx, "Exporting pulse-analysis PNG artifacts...")
    export_pulse_pngs(ctx.output, context, per_beat_result, log=getattr(ctx, "log", None))


def _segment_ring_settings() -> SegmentRingSettings:
    inner = SEGMENT_INNER_RADIUS_FRAC
    outer = SEGMENT_OUTER_RADIUS_FRAC
    count = SEGMENT_RING_COUNT
    width = _ring_width(inner, outer, count)
    return SegmentRingSettings(
        inner,
        outer,
        width,
        count,
        SEGMENT_LENGTH_FRAC,
    )


def _ring_width(inner: float, outer: float, count: int) -> float:
    return max((outer - inner) / float(count), np.finfo(np.float32).eps)


def _context_attrs(
    source_data: WaveformShapeSourceData,
    timing: HolodopplerTiming,
    harmonic_count: int,
    analysis_source: str,
    heartbeat,
) -> dict[str, object]:
    output_paths = EyeFlowOutputPaths.active()
    analysis_paths = output_paths.analysis
    dependency_chain = (
        ["dopplerview.h5.analysis"]
        if analysis_source == "dopplerview_h5_analysis"
        else [
            "holodoppler.h5.moment0_moment2",
            "dopplerview.h5.segmentation",
            "eyeflow.dopplerview_analysis.recomputed",
        ]
    )
    return {
        "dependency_chain": dependency_chain + [
            "blood_flow_velocity.signal_analysis.heartbeat.spectral",
            "blood_flow_velocity.signal_analysis.per_beat.signal",
            "blood_flow_velocity.signal_analysis.per_beat.runner",
        ],
        "analysis_source": analysis_source,
        "output_schema": output_paths.name,
        "moment0_flat_field_source": "dopplerview_recomputed_from_raw",
        "moment2_flat_field_source": "dopplerview_recomputed_from_raw",
        "flat_field_gaussian_ratio": float(source_data.flat_field_gaussian_ratio),
        "flat_field_border": float(source_data.flat_field_border),
        "arterial_velocity_signal_path": (
            analysis_paths.retinal_artery_velocity_signal_filtered
        ),
        "venous_velocity_signal_path": (
            analysis_paths.retinal_vein_velocity_signal_filtered
        ),
        "systolic_peak_indexes_path": analysis_paths.beat_indices,
        "beat_period_seconds_path": output_paths.beat_period_seconds,
        "heart_rate_hz": float(heartbeat.heart_rate_hz),
        "heart_rate_bpm": float(heartbeat.heart_rate_bpm),
        "heart_rate_ste_hz": float(heartbeat.heart_rate_ste_hz),
        "heart_rate_ste_bpm": float(heartbeat.heart_rate_ste_bpm),
        "sampling_freq": float(timing.sampling_freq),
        "batch_stride": float(timing.batch_stride),
        "dt_seconds": float(timing.dt_seconds),
        "band_limited_signal_harmonic_count": int(harmonic_count),
        "filter_velocity_signals": bool(LEGACY_FILTER_VELOCITY_SIGNALS),
        "velocity_signal_lowpass_hz": float(LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ),
    }


def _pack_meta_outputs(context: WaveformShapeMetricsContext) -> dict[str, object]:
    schema = EyeFlowOutputPaths.active()
    timing = context.source_data.timing
    return {
        f"{schema.meta_root}/SamplingFrequencyHz/value": (
            np.float32(timing.sampling_freq),
            {"unit": "Hz"},
        ),
        f"{schema.meta_root}/BatchStride/value": np.float32(timing.batch_stride),
        f"{schema.meta_root}/FrameIntervalSeconds/value": (
            np.float32(timing.dt_seconds),
            {"unit": "s"},
        ),
    }
