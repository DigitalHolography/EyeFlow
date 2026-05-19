"""Orchestrate the waveform-shape metrics sandbox pipeline."""

from collections.abc import Mapping

from calculations.blood_flow_velocity import (
    PerBeatAnalysisInput,
    SegmentRingSettings,
    run_per_beat_analysis,
    segment_velocity_results,
)
from input_output import EyeFlowOutputPaths
from pipeline_engine.imports import (
    HolodopplerTiming,
    np,
    read_int_setting,
)

from .constants import (
    LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
    LEGACY_FILTER_VELOCITY_SIGNALS,
    LEGACY_SEGMENT_INNER_RADIUS_FRAC,
    LEGACY_SEGMENT_OUTER_RADIUS_FRAC,
    LEGACY_SEGMENT_RING_COUNT,
    LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
)
from .dopplerview import run_dopplerview_analysis
from .branch_identity_debug import export_branch_identity_stage_pngs
from .models import WaveformShapeMetricsContext
from .outputs import pack_dopplerview_analysis_outputs, pack_velocity_per_beat_outputs
from .sources import WaveformShapeSourceData, WaveformShapeSources


def run_waveform_shape_metrics(ctx) -> tuple[dict[str, object], dict[str, object]]:
    ctx.require_inputs("hd", "dv")

    context = _build_waveform_shape_metrics_context(ctx)
    ctx.set_var("waveform_shape_metrics_context", context)
    ctx.set_var("dopplerview_analysis", context.dopplerview_analysis)

    per_beat_result = run_per_beat_analysis(context.per_beat_analysis)
    metrics = pack_dopplerview_analysis_outputs(context.dopplerview_analysis)
    metrics.update(pack_velocity_per_beat_outputs(per_beat_result))
    ctx.set_var("velocity_per_beat_result", per_beat_result)

    return metrics, context.attrs


def _build_waveform_shape_metrics_context(ctx) -> WaveformShapeMetricsContext:
    source_data = WaveformShapeSources.from_context(ctx).load()
    timing = source_data.timing
    dopplerview_analysis = run_dopplerview_analysis(source_data)
    harmonic_count = _band_limited_harmonic_count(ctx)

    return WaveformShapeMetricsContext(
        per_beat_analysis=_per_beat_input_from_analysis(
            dopplerview_analysis,
            source_data,
            timing,
            harmonic_count,
            ctx,
        ),
        dopplerview_analysis=dopplerview_analysis,
        attrs=_context_attrs(source_data, timing, harmonic_count),
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
) -> PerBeatAnalysisInput:
    ring_settings = _segment_ring_settings(source_data)
    artery_segments, vein_segments = _segment_velocity_inputs(
        analysis,
        source_data,
        ring_settings,
        ctx,
    )
    return PerBeatAnalysisInput(
        arterial_velocity_signal=np.asarray(
            analysis["retinal_artery_velocity_signal"],
            dtype=np.float32,
        ),
        venous_velocity_signal=np.asarray(
            analysis["retinal_vein_velocity_signal"],
            dtype=np.float32,
        ),
        systolic_acceleration_peak_indexes=np.asarray(
            analysis["beat_indices"],
            dtype=np.int32,
        ),
        band_limited_signal_harmonic_count=harmonic_count,
        dt_seconds=timing.dt_seconds,
        arterial_velocity_segments=artery_segments,
        venous_velocity_segments=vein_segments,
        beat_period_seconds=np.asarray(analysis["time_per_beat"], dtype=np.float32),
        index_base=source_data.provenance["beat_index_base"],
    )


def _segment_velocity_inputs(
    analysis: Mapping[str, object],
    source_data: WaveformShapeSourceData,
    ring_settings: SegmentRingSettings,
    ctx,
) -> tuple[np.ndarray, np.ndarray]:
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
    return artery.velocity, vein.velocity


def _export_branch_identity_debug(
    ctx,
    stages,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    prefix: str,
) -> None:
    if ctx.output is None:
        return
    export_branch_identity_stage_pngs(
        ctx.output,
        stages,
        prefix,
        optic_disc_center,
        ring_settings,
    )


def _segment_ring_settings(source_data: WaveformShapeSourceData) -> SegmentRingSettings:
    inner = _positive_float(source_data.peripapillary_inner_radius, LEGACY_SEGMENT_INNER_RADIUS_FRAC)
    outer = _positive_float(source_data.peripapillary_outer_radius, LEGACY_SEGMENT_OUTER_RADIUS_FRAC)
    count = _positive_int(source_data.peripapillary_ring_count, LEGACY_SEGMENT_RING_COUNT)
    width = _ring_width(source_data.peripapillary_ring_width, inner, outer, count)
    length = _positive_float(source_data.segment_length, width)
    return SegmentRingSettings(inner, outer, width, count, length)


def _ring_width(value, inner: float, outer: float, count: int) -> float:
    width = _positive_float(value, 0.0)
    if width > 0:
        return width
    return max((outer - inner) / float(count), np.finfo(np.float32).eps)


def _positive_float(value, default: float) -> float:
    if value is None:
        return float(default)
    scalar = float(np.asarray(value).reshape(-1)[0])
    return scalar if scalar > 0 else float(default)


def _positive_int(value, default: int) -> int:
    if value is None:
        return int(default)
    scalar = int(np.asarray(value).reshape(-1)[0])
    return scalar if scalar > 0 else int(default)


def _context_attrs(
    source_data: WaveformShapeSourceData,
    timing: HolodopplerTiming,
    harmonic_count: int,
) -> dict[str, object]:
    analysis_paths = EyeFlowOutputPaths.active().analysis
    return {
        "dependency_chain": [
            "dopplerview.vessel_velocity_estimator",
            "dopplerview.arterial_waveform_analysis",
            "blood_flow_velocity.signal_analysis.signal.per_beat_signal",
            "blood_flow_velocity.signal_analysis.signal.per_beat",
        ],
        "analysis_source": "computed_dopplerview_steps",
        "arterial_velocity_signal_path": analysis_paths.retinal_artery_velocity_signal,
        "venous_velocity_signal_path": analysis_paths.retinal_vein_velocity_signal,
        "systolic_peak_indexes_path": analysis_paths.beat_indices,
        "beat_period_seconds_path": analysis_paths.time_per_beat,
        "segment_velocity_source": (
            "CrossSection/labelVesselBranches.m-style branch labeling; "
            "CrossSection/generateCrossSectionSignals.m-style branch/circle "
            "velocity extraction from retinal_velocity_array"
        ),
        "sampling_freq": float(timing.sampling_freq),
        "batch_stride": float(timing.batch_stride),
        "dt_seconds": float(timing.dt_seconds),
        "band_limited_signal_harmonic_count": int(harmonic_count),
        "filter_velocity_signals": bool(LEGACY_FILTER_VELOCITY_SIGNALS),
        "velocity_signal_lowpass_hz": float(LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ),
    }
