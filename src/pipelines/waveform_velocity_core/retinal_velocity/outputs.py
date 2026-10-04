"""Output packing for shared retinal velocity products."""

from collections.abc import Iterable, Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths


def pack_retinal_velocity_outputs(
    velocity_analysis: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack shared frequency-map, heartbeat, and provenance analysis outputs."""

    schema = _resolve_output_paths(output_paths)
    paths = schema.analysis
    frequency_attrs = _frequency_attrs(velocity_analysis)
    metrics = {
        paths.fRMS_avg: metric_value(
            velocity_analysis["fRMS_avg"],
            attrs=frequency_attrs,
        ),
        paths.fRMS_bkg_avg: metric_value(
            velocity_analysis["fRMS_bkg_avg"],
            attrs=frequency_attrs,
        ),
        paths.beat_indices: metric_value(velocity_analysis["beat_indices"]),
        paths.time_per_beat: metric_value(
            velocity_analysis["time_per_beat"],
            unit="s",
        ),
    }
    heartbeat = velocity_analysis.get("_heartbeat_analysis_result")
    spectral = getattr(heartbeat, "spectral", None)
    if spectral is not None:
        heartbeat_paths = schema.heartbeat
        metrics.update(
            {
                heartbeat_paths.spectral_fundamental_frequency_hz: metric_value(
                    spectral.fundamental_hz,
                    unit="Hz",
                ),
                heartbeat_paths.spectral_heart_rate_bpm: metric_value(
                    spectral.heart_rate_bpm,
                    unit="bpm",
                ),
                heartbeat_paths.spectral_heart_rate_standard_error_bpm: metric_value(
                    spectral.heart_rate_ste_bpm,
                    unit="bpm",
                ),
                heartbeat_paths.spectral_period_seconds: metric_value(
                    spectral.period_seconds,
                    unit="s",
                ),
            }
        )
    return metrics


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
    attrs: Mapping[str, object] | None = None,
):
    output_attrs: dict[str, object] = dict(attrs or {})
    if unit:
        output_attrs["unit"] = unit
    if dim_desc:
        output_attrs["dimDesc"] = list(dim_desc)
    data = metric_data(data)
    return (data, output_attrs) if output_attrs else data


def _frequency_attrs(
    velocity_analysis: Mapping[str, object],
) -> dict[str, object]:
    from velocity_calibration import calibration_attrs_from_metadata

    attrs = calibration_attrs_from_metadata(velocity_analysis)
    attrs.update(
        {
            "unit": "Hz",
            "quantity": "rms_frequency",
            "velocity_estimation_method": str(
                velocity_analysis.get(
                    "velocity_estimation_method",
                    "doppler_moments",
                )
            ),
        }
    )
    return attrs


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
