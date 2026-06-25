"""Output packing for the waveform-shape metrics pipeline."""

from collections.abc import Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths


def pack_dopplerview_analysis_outputs(
    analysis: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    paths = _resolve_output_paths(output_paths).analysis
    metrics = {
        paths.retinal_artery_velocity_signal: metric_value(
            analysis["retinal_artery_velocity_signal"],
            unit="mm/s",
        ),
        paths.retinal_vein_velocity_signal: metric_value(
            analysis["retinal_vein_velocity_signal"],
            unit="mm/s",
        ),
        paths.velocity_map_avg: metric_value(analysis["velocity_map_avg"]),
        paths.fRMS_avg: metric_value(analysis["fRMS_avg"]),
        paths.fRMS_bkg_avg: metric_value(analysis["fRMS_bkg_avg"]),
        paths.velocitysignal_per_beat: metric_value(
            analysis["retinal_artery_velocity_signal_filtered_perbeat"],
            unit="mm/s",
        ),
        paths.velocitysignal_filtered: metric_value(
            analysis["retinal_artery_velocity_signal_filtered"],
            unit="mm/s",
        ),
        paths.beat_indices: metric_value(analysis["beat_indices"]),
        paths.time_per_beat: metric_value(
            analysis["time_per_beat"],
            unit="s",
        ),
    }
    if paths.retinal_velocity_array is not None:
        metrics[paths.retinal_velocity_array] = metric_value(
            analysis["retinal_vessel_velocity"],
            unit="mm/s",
        )
    return metrics


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


def metric_value(data, *, unit: str | None = None):
    attrs: dict[str, object] = {}
    if unit:
        attrs["unit"] = unit
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
