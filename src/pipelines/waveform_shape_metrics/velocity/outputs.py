"""Output packing for per-beat velocity calculations."""

from collections.abc import Iterable

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisResult
from calculations.blood_flow_velocity.analysis_preparation.beat_detection import (
    RECOVERY_GAP_RATIO_MAX,
    RECOVERY_GAP_RATIO_MIN,
    RECOVERY_MEDIAN_PRIMARY_HEIGHT_FRACTION,
    RECOVERY_MEDIAN_PRIMARY_PROMINENCE_FRACTION,
    RECOVERY_ORIGINAL_THRESHOLD_FRACTION,
    RECOVERY_SEARCH_RADIUS_PERIOD_FRACTION,
)
from input_output.schema import EyeFlowOutputPaths, VelocityPerBeatOutputPaths


def pack_velocity_per_beat_outputs(
    result: PerBeatAnalysisResult,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    schema = _resolve_output_paths(output_paths)
    metrics = {
        schema.beat_period_idx: metric_value(
            _matlab_row_vector(result.beat_period_idx),
            unit="frame",
            dim_desc=("row", "beat"),
        ),
        schema.beat_period_seconds: metric_value(
            _matlab_row_vector(result.beat_period_seconds),
            unit="s",
            dim_desc=("row", "beat"),
        ),
    }
    metrics.update(_pack_vessel_outputs(schema.artery_per_beat, result.artery))
    metrics.update(_pack_vessel_outputs(schema.vein_per_beat, result.vein))
    metrics.update(_pack_quality_outputs(schema, result))
    return metrics


def _pack_quality_outputs(
    schema: EyeFlowOutputPaths,
    result: PerBeatAnalysisResult,
) -> dict[str, object]:
    paths = schema.beat_quality
    quality = result.quality
    if paths is None or quality is None:
        return {}

    accepted_indexes = np.flatnonzero(quality.accepted_mask).astype(np.int32)
    boundary_index_attrs = (
        {}
        if result.index_base is None
        else {"indexBase": np.int32(result.index_base)}
    )
    settings_attrs = {
        "scoreAcceptanceThreshold": np.float32(1.0),
        "rawBandlimitedResidualLimit": np.float32(
            quality.settings.raw_bandlimited_residual_limit
        ),
        "pairedTemplateDistanceLimit": np.float32(
            quality.settings.paired_template_distance_limit
        ),
        "minimumTemplateBeats": np.int32(quality.settings.minimum_template_beats),
        "minimumPeriodRatio": np.float32(quality.settings.minimum_period_ratio),
        "maximumPeriodRatio": np.float32(quality.settings.maximum_period_ratio),
    }
    flag_attrs = {
        **settings_attrs,
        "flagDefinitions": [
            "1=duration_too_short",
            "2=duration_too_long",
            "4=artery_raw_bandlimited_residual",
            "8=vein_raw_bandlimited_residual",
            "16=artery_template_mismatch",
            "32=vein_template_mismatch",
            "64=insufficient_or_nonfinite_waveform",
        ],
    }
    return {
        paths.candidate_start_index: metric_value(
            _matlab_row_vector(result.candidate_start_indexes),
            unit="frame",
            dim_desc=("row", "candidate_beat"),
            extra_attrs=boundary_index_attrs,
        ),
        paths.candidate_stop_index: metric_value(
            _matlab_row_vector(result.candidate_stop_indexes),
            unit="frame",
            dim_desc=("row", "candidate_beat"),
            extra_attrs=boundary_index_attrs,
        ),
        paths.accepted_mask: metric_value(
            _matlab_row_vector(quality.accepted_mask),
            dim_desc=("row", "candidate_beat"),
        ),
        paths.accepted_candidate_index: metric_value(
            _matlab_row_vector(accepted_indexes),
            dim_desc=("row", "accepted_beat"),
            extra_attrs={"indexBase": np.int32(0)},
        ),
        paths.rejection_flags: metric_value(
            _matlab_row_vector(quality.rejection_flags),
            dim_desc=("row", "candidate_beat"),
            extra_attrs=flag_attrs,
        ),
        paths.quality_score: metric_value(
            _matlab_row_vector(quality.quality_score),
            dim_desc=("row", "candidate_beat"),
            extra_attrs=settings_attrs,
        ),
        paths.period_ratio: metric_value(
            _matlab_row_vector(quality.period_ratio),
            dim_desc=("row", "candidate_beat"),
        ),
        paths.boundary_recovered_mask: metric_value(
            _matlab_row_vector(result.recovered_boundary_mask),
            dim_desc=("row", "boundary"),
            extra_attrs={
                "gapRatioRange": np.asarray(
                    [RECOVERY_GAP_RATIO_MIN, RECOVERY_GAP_RATIO_MAX],
                    dtype=np.float32,
                ),
                "midpointSearchRadiusPeriodFraction": np.float32(
                    RECOVERY_SEARCH_RADIUS_PERIOD_FRACTION
                ),
                "minimumOriginalThresholdFraction": np.float32(
                    RECOVERY_ORIGINAL_THRESHOLD_FRACTION
                ),
                "minimumMedianPrimaryHeightFraction": np.float32(
                    RECOVERY_MEDIAN_PRIMARY_HEIGHT_FRACTION
                ),
                "minimumMedianPrimaryProminenceFraction": np.float32(
                    RECOVERY_MEDIAN_PRIMARY_PROMINENCE_FRACTION
                ),
            },
        ),
        paths.nominal_period_samples: metric_value(
            result.nominal_period_samples,
            unit="frame",
        ),
        paths.artery_raw_bandlimited_residual: metric_value(
            _matlab_row_vector(quality.artery.raw_bandlimited_residual),
            dim_desc=("row", "candidate_beat"),
        ),
        paths.vein_raw_bandlimited_residual: metric_value(
            _matlab_row_vector(quality.vein.raw_bandlimited_residual),
            dim_desc=("row", "candidate_beat"),
        ),
        paths.artery_template_correlation: metric_value(
            _matlab_row_vector(quality.artery.template_correlation),
            dim_desc=("row", "candidate_beat"),
        ),
        paths.vein_template_correlation: metric_value(
            _matlab_row_vector(quality.vein.template_correlation),
            dim_desc=("row", "candidate_beat"),
        ),
        paths.candidate_count: metric_value(int(quality.accepted_mask.size)),
        paths.accepted_count: metric_value(int(np.sum(quality.accepted_mask))),
        paths.rejected_count: metric_value(int(np.sum(~quality.accepted_mask))),
        paths.recovered_boundary_count: metric_value(
            int(np.sum(result.recovered_boundary_mask))
        ),
    }


def _pack_vessel_outputs(
    paths: VelocityPerBeatOutputPaths,
    vessel,
) -> dict[str, object]:
    signal = vessel.signal
    metrics = {
        paths.velocity_signal: metric_value(
            signal.velocity_signal_per_beat,
            unit="mm/s",
            dim_desc=("beat", "sample"),
        ),
        paths.velocity_signal_fft_abs: metric_value(
            np.abs(signal.velocity_signal_per_beat_fft),
            unit="a.u.",
            dim_desc=("beat", "frequency_bin"),
        ),
        paths.velocity_signal_fft_arg: metric_value(
            np.angle(signal.velocity_signal_per_beat_fft),
            unit="rad",
            dim_desc=("beat", "frequency_bin"),
        ),
        paths.velocity_signal_band_limited: metric_value(
            signal.velocity_signal_per_beat_band_limited,
            unit="mm/s",
            dim_desc=("beat", "sample"),
        ),
    }
    if vessel.segments is not None:
        metrics.update(_pack_vessel_segment_outputs(paths, vessel.segments))
    return metrics


def _pack_vessel_segment_outputs(
    paths: VelocityPerBeatOutputPaths,
    segments,
) -> dict[str, object]:
    if paths.segment_velocity_signal is None:
        return {}
    if paths.segment_velocity_signal_band_limited is None:
        return {}
    return {
        paths.segment_velocity_signal: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment,
            unit="mm/s",
        ),
        paths.segment_velocity_signal_band_limited: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment_band_limited,
            unit="mm/s",
        ),
    }


def _segment_metric_value(data, *, unit: str):
    value = metric_data(data)
    if value.ndim != 4:
        raise ValueError(
            "segment per-beat outputs must have shape "
            "(sample, beat, branch, radius)."
        )
    return metric_value(
        value,
        unit=unit,
        dim_desc=("sample", "beat", "branch", "radius"),
    )


def _matlab_row_vector(data) -> np.ndarray:
    return np.asarray(data).reshape(1, -1)


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


def metric_value(
    data,
    *,
    unit: str | None = None,
    dim_desc: Iterable[str] | None = None,
    extra_attrs: dict[str, object] | None = None,
):
    attrs: dict[str, object] = dict(extra_attrs or {})
    if unit:
        attrs["unit"] = unit
    if dim_desc:
        attrs["dimDesc"] = list(dim_desc)
    data = metric_data(data)
    return (data, attrs) if attrs else data


def metric_data(data):
    if isinstance(data, bool):
        return data
    if isinstance(data, float):
        return np.float32(data)
    if isinstance(data, int):
        return np.int32(data)
    if isinstance(data, complex):
        return np.complex64(data)
    value = np.asarray(data)
    if value.dtype.kind == "f":
        return value.astype(np.float32, copy=False)
    if value.dtype.kind == "c":
        return value.astype(np.complex64, copy=False)
    if value.dtype.kind == "i":
        return value.astype(np.int32, copy=False)
    if value.dtype.kind == "u":
        return value.astype(np.uint32, copy=False)
    return value
