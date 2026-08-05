"""Build shared DopplerView, spatial, and segment-analysis state."""

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
from utils.logger import Logger

from .dopplerview.constants import (
    LEGACY_FILTER_VELOCITY_SIGNALS,
    LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
)
from .constants import (
    LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
    SEGMENT_RING_COUNT,
    SEGMENT_INNER_RADIUS_FRAC,
    SEGMENT_LENGTH_FRAC,
    SEGMENT_OUTER_RADIUS_FRAC,
)
from .dopplerview.outputs import (
    pack_dopplerview_shared_outputs,
)
from .dopplerview.runner import run_dopplerview_analysis
from .scratch import waveform_scratch_h5
from .sources import WaveformVelocitySourceData, WaveformVelocitySources
from .branch_identity_debug import export_branch_identity_stage_pngs
from .figures import export_pulse_pngs
from .per_beat import run_velocity_per_beat_metrics
from .segmentation import pack_segmentation_outputs


WAVEFORM_CONTEXT_STATE = "waveform_velocity_context"
VELOCITY_PER_BEAT_RESULT_STATE = "velocity_per_beat_result"
VELOCITY_PER_BEAT_OUTPUTS_STATE = "velocity_per_beat_outputs"


@dataclass(frozen=True)
class WaveformVelocityCoreContext:
    source_data: WaveformVelocitySourceData
    per_beat_analysis: PerBeatAnalysisInput
    artery_segment_result: CrossSectionSignalResult
    vein_segment_result: CrossSectionSignalResult
    dopplerview_analysis: dict[str, object]
    attrs: dict[str, object]


def run_waveform_velocity_core(
    ctx,
) -> tuple[dict[str, object], dict[str, object]]:
    """Run the shared DopplerView and spatial velocity foundation once."""
    ctx.require_inputs("hd", "dv")

    with waveform_scratch_h5(ctx) as scratch_h5:
        Logger.log("Starting waveform velocity core context build...")
        segments_required = _segments_required(ctx)
        context = _build_waveform_velocity_core_context(
            ctx,
            scratch_h5,
            segments_required=segments_required,
        )
        metrics = pack_dopplerview_shared_outputs(context.dopplerview_analysis)
        metrics.update(_pack_meta_outputs(context))
        metrics.update(
            pack_segmentation_outputs(
                context.source_data,
                context.artery_segment_result,
                context.vein_segment_result,
            )
        )
        ctx.state.set(WAVEFORM_CONTEXT_STATE, context)

        if _per_beat_required(ctx):
            Logger.log("Starting shared per-beat velocity analysis...")
            per_beat_result, velocity_outputs = run_velocity_per_beat_metrics(context)
            ctx.state.set(VELOCITY_PER_BEAT_RESULT_STATE, per_beat_result)
            ctx.state.set(VELOCITY_PER_BEAT_OUTPUTS_STATE, velocity_outputs)
            if _pulse_pngs_required(ctx):
                _export_pulse_pngs(ctx, context, per_beat_result)

    return metrics, context.attrs


def _per_beat_required(ctx) -> bool:
    velocity_options = ctx.options_for("waveform_velocity")
    metric_options = ctx.options_for("waveform_shape_metrics")
    if ctx.pipeline_scheduled("pdf_report"):
        return True
    if ctx.pipeline_scheduled("waveform_velocity"):
        return bool(
            {"per_beat", "hemifield"} & velocity_options
            or "hemifield" in metric_options
        )
    return bool(
        metric_options and ctx.pipeline_scheduled("waveform_shape_metrics")
    )


def _segments_required(ctx) -> bool:
    """Return whether any selected product needs spatial vessel segments."""
    velocity_options = ctx.options_for("waveform_velocity")
    metric_options = ctx.options_for("waveform_shape_metrics")
    if ctx.pipeline_scheduled("waveform_velocity"):
        return bool(
            {"segments", "velocity_profiles", "hemifield"} & velocity_options
            or "hemifield" in metric_options
        )
    return bool(
        {"segments", "hemifield"} & metric_options
    )


def _pulse_pngs_required(ctx) -> bool:
    return bool(
        ctx.pipeline_scheduled("pdf_report")
        or (
            ctx.pipeline_scheduled("waveform_velocity")
            and ctx.option_enabled("per_beat", pipeline="waveform_velocity")
        )
    )


def _build_waveform_velocity_core_context(
    ctx,
    scratch_h5,
    *,
    segments_required: bool,
) -> WaveformVelocityCoreContext:
    Logger.log("Starting waveform source loading...")
    source_data = WaveformVelocitySources.from_context(ctx).load()
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
        segments_required=segments_required,
    )

    return WaveformVelocityCoreContext(
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
    source_data: WaveformVelocitySourceData,
    ctx,
    scratch_h5,
) -> tuple[dict[str, object], str]:
    Logger.log(
        "Building EyeFlow velocity analysis from HD moments and DV segmentation.",
    )
    return (
        run_dopplerview_analysis(
            source_data,
            scratch_h5,
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
    source_data: WaveformVelocitySourceData,
    timing: HolodopplerTiming,
    harmonic_count: int,
    ctx,
    *,
    segments_required: bool,
) -> tuple[
    PerBeatAnalysisInput,
    CrossSectionSignalResult | None,
    CrossSectionSignalResult | None,
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
        arterial_velocity_segments=_waveform_segment_input(
            artery_segments,
            include_segments=segments_required,
        ),
        venous_velocity_segments=_waveform_segment_input(
            vein_segments,
            include_segments=segments_required,
        ),
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
    source_data: WaveformVelocitySourceData,
    ring_settings: SegmentRingSettings,
    ctx,
) -> tuple[CrossSectionSignalResult, CrossSectionSignalResult]:
    Logger.log("Starting segment velocity extraction...")
    artery, vein = segment_velocity_results(
        analysis["retinal_vessel_velocity"],
        source_data.retinal_artery_mask,
        source_data.retinal_vein_mask,
        source_data.optic_disc_center,
        ring_settings,
        source_data.cross_section_settings,
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


def _waveform_segment_input(
    result: CrossSectionSignalResult | None,
    *,
    include_segments: bool,
) -> np.ndarray | None:
    if not include_segments or result is None or result.branch_ids.size == 0:
        return None
    return result.velocity


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


def _export_pulse_pngs(ctx, context: WaveformVelocityCoreContext, per_beat_result) -> None:
    if not ctx.output.available:
        return
    Logger.log("Exporting pulse-analysis PNG artifacts...")
    export_pulse_pngs(ctx.output, context, per_beat_result)


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
    source_data: WaveformVelocitySourceData,
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
            analysis_paths.retinal_artery_velocity_signal_band_limited
        ),
        "venous_velocity_signal_path": (
            analysis_paths.retinal_vein_velocity_signal_band_limited
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


def _pack_meta_outputs(context: WaveformVelocityCoreContext) -> dict[str, object]:
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
