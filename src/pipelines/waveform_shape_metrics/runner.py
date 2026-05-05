from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from domain.blood_flow_velocity import PerBeatAnalysisInput, run_per_beat_analysis
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
    return PerBeatAnalysisInput(
        arterial_velocity_signal=np.asarray(
            analysis["retinal_artery_velocity_signal"],
            dtype=np.float64,
        ),
        venous_velocity_signal=np.asarray(
            analysis["retinal_vein_velocity_signal"],
            dtype=np.float64,
        ),
        systolic_acceleration_peak_indexes=np.asarray(
            analysis["beat_indices"],
            dtype=np.int64,
        ),
        band_limited_signal_harmonic_count=harmonic_count,
        dt_seconds=timing.dt_seconds,
        beat_period_seconds=np.asarray(analysis["time_per_beat"], dtype=np.float64),
        index_base=systolic_index_base_for_path(
            DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path("beat_indices")
        ),
    )


def _context_attrs(
    timing: HolodopplerTiming,
    harmonic_count: int,
) -> dict[str, object]:
    return {
        "dependency_chain": [
            "DopplerView vessel_velocity_estimator",
            "DopplerView arterial_waveform_analysis",
            "perBeatSignalAnalysis",
            "perBeatAnalysis",
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
        "sampling_freq": float(timing.sampling_freq),
        "batch_stride": float(timing.batch_stride),
        "dt_seconds": float(timing.dt_seconds),
        "band_limited_signal_harmonic_count": int(harmonic_count),
    }
