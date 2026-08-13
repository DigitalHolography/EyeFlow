"""Output packing for shared DopplerView analysis products."""

from collections.abc import Iterable, Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths


def pack_dopplerview_shared_outputs(
    analysis: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack shared frequency-map, heartbeat, and provenance analysis outputs."""

    schema = _resolve_output_paths(output_paths)
    paths = schema.analysis
    metrics = {
        paths.fRMS_avg: metric_value(analysis["fRMS_avg"]),
        paths.fRMS_bkg_avg: metric_value(analysis["fRMS_bkg_avg"]),
        paths.beat_indices: metric_value(analysis["beat_indices"]),
        paths.time_per_beat: metric_value(
            analysis["time_per_beat"],
            unit="s",
        ),
    }
    if paths.velocity_map_avg is not None:
        metrics[paths.velocity_map_avg] = metric_value(
            analysis["velocity_map_avg"],
            unit="mm/s",
        )
    if paths.retinal_velocity_array is not None:
        metrics[paths.retinal_velocity_array] = metric_value(
            analysis["retinal_vessel_velocity"],
            unit="mm/s",
        )
    heartbeat = analysis.get("_heartbeat_analysis_result")
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
):
    attrs: dict[str, object] = {}
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
