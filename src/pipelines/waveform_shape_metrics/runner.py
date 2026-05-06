"""Orchestrate the waveform-shape metrics sandbox pipeline."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from calculations.blood_flow_velocity import (
    PerBeatAnalysisInput,
    run_per_beat_analysis,
)
from input_output import (
    DOPPLER_VIEW_ANALYSIS_SCHEMA,
    pack_dopplerview_analysis_outputs,
    pack_velocity_per_beat_outputs,
    systolic_index_base_for_path,
)
from input_output.input_access import (
    HolodopplerTiming,
    read_int_setting,
    resolve_holodoppler_timing,
)

from .constants import LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT
from .dopplerview import run_dopplerview_analysis
from .models import WaveformShapeMetricsContext


def run_waveform_shape_metrics(ctx) -> tuple[dict[str, object], dict[str, object]]:
    if ctx.hd.h5file is None or ctx.dv.h5file is None:
        raise ValueError("waveform_shape_metrics requires both HD and DV inputs.")

    context = _build_waveform_shape_metrics_context(ctx)
    per_beat_result = run_per_beat_analysis(context.per_beat_analysis)
    metrics = pack_dopplerview_analysis_outputs(context.dopplerview_analysis)
    metrics.update(pack_velocity_per_beat_outputs(per_beat_result))
    return metrics, context.attrs


def _build_waveform_shape_metrics_context(ctx) -> WaveformShapeMetricsContext:
    timing = resolve_holodoppler_timing(ctx)
    dopplerview_analysis = run_dopplerview_analysis(ctx, timing)
    harmonic_count = _band_limited_harmonic_count(ctx)

    return WaveformShapeMetricsContext(
        per_beat_analysis=_per_beat_input_from_analysis(
            dopplerview_analysis,
            timing,
            harmonic_count,
        ),
        dopplerview_analysis=dopplerview_analysis,
        attrs=_context_attrs(timing, harmonic_count),
    )


def _band_limited_harmonic_count(ctx) -> int:
    return read_int_setting(
        ctx,
        default=LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
        keys=("BandLimitedSignalHarmonicCount", "band_limited_signal_harmonic_count"),
    )


def _per_beat_input_from_analysis(
    analysis: Mapping[str, object],
    timing: HolodopplerTiming,
    harmonic_count: int,
) -> PerBeatAnalysisInput:
    artery_segments, vein_segments = _segment_velocity_inputs(analysis)
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
        index_base=systolic_index_base_for_path(
            DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path("beat_indices")
        ),
    )


def _segment_velocity_inputs(
    analysis: Mapping[str, object],
) -> tuple[np.ndarray, np.ndarray]:
    velocity = np.asarray(analysis["retinal_vessel_velocity"], dtype=np.float32)
    labels = analysis.get("retinal_labeled_vessels")
    artery_mask = np.asarray(analysis["retinal_artery_mask"], dtype=bool)
    vein_mask = np.asarray(analysis["retinal_vein_mask"], dtype=bool)
    if labels is None:
        return (
            _single_vessel_segment_signal(velocity, artery_mask),
            _single_vessel_segment_signal(velocity, vein_mask),
        )
    label_array = np.asarray(labels, dtype=np.int32)
    return (
        _labeled_segment_velocity_signals(velocity, label_array, artery_mask),
        _labeled_segment_velocity_signals(velocity, label_array, vein_mask),
    )


def _single_vessel_segment_signal(
    velocity: np.ndarray,
    vessel_mask: np.ndarray,
) -> np.ndarray:
    if not np.any(vessel_mask):
        return np.full((1, 1, velocity.shape[0]), np.nan, dtype=np.float32)
    signal = np.nanmean(velocity[:, vessel_mask], axis=1)
    return signal.astype(np.float32, copy=False).reshape(1, 1, -1)


def _labeled_segment_velocity_signals(
    velocity: np.ndarray,
    labels: np.ndarray,
    vessel_mask: np.ndarray,
) -> np.ndarray:
    branch_ids = _vessel_branch_ids(labels, vessel_mask)
    if branch_ids.size == 0:
        return _single_vessel_segment_signal(velocity, vessel_mask)
    segments = np.full((1, branch_ids.size, velocity.shape[0]), np.nan, dtype=np.float32)
    for branch_index, branch_id in enumerate(branch_ids):
        branch_mask = (labels == int(branch_id)) & vessel_mask
        segments[0, branch_index, :] = np.nanmean(
            velocity[:, branch_mask],
            axis=1,
        ).astype(np.float32, copy=False)
    return segments


def _vessel_branch_ids(labels: np.ndarray, vessel_mask: np.ndarray) -> np.ndarray:
    branch_ids = np.unique(labels[vessel_mask & (labels > 0)])
    return branch_ids.astype(np.int32, copy=False)


def _context_attrs(
    timing: HolodopplerTiming,
    harmonic_count: int,
) -> dict[str, object]:
    return {
        "dependency_chain": [
            "dopplerview.vessel_velocity_estimator",
            "dopplerview.arterial_waveform_analysis",
            "blood_flow_velocity.per_beat_signal",
            "blood_flow_velocity.per_beat",
        ],
        "analysis_source": "computed_dopplerview_steps",
        "arterial_velocity_signal_path": DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path(
            "retinal_artery_velocity_signal"
        ),
        "venous_velocity_signal_path": DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path(
            "retinal_vein_velocity_signal"
        ),
        "systolic_peak_indexes_path": DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path(
            "beat_indices"
        ),
        "beat_period_seconds_path": DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path(
            "time_per_beat"
        ),
        "segment_velocity_source": (
            "retinal_velocity_array averaged by DV labeled branches; "
            "falls back to one whole-vessel segment when labels are absent"
        ),
        "sampling_freq": float(timing.sampling_freq),
        "batch_stride": float(timing.batch_stride),
        "dt_seconds": float(timing.dt_seconds),
        "band_limited_signal_harmonic_count": int(harmonic_count),
    }
