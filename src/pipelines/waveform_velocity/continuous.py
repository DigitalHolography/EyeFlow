"""Whole-vessel continuous velocity output packing."""

from collections.abc import Mapping

from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_velocity_core.dopplerview.outputs import metric_value


def pack_continuous_velocity_outputs(
    analysis: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack raw and band-limited artery and vein velocity signals."""
    schema = _resolve_output_paths(output_paths)
    paths = schema.analysis
    metrics = {
        paths.retinal_artery_velocity_signal: metric_value(
            analysis["retinal_artery_velocity_signal"],
            unit="mm/s",
        ),
        paths.retinal_vein_velocity_signal: metric_value(
            analysis["retinal_vein_velocity_signal"],
            unit="mm/s",
        ),
        paths.retinal_artery_velocity_signal_band_limited: metric_value(
            analysis["retinal_artery_velocity_signal_filtered"],
            unit="mm/s",
        ),
        paths.retinal_vein_velocity_signal_band_limited: metric_value(
            analysis["retinal_vein_velocity_signal_filtered"],
            unit="mm/s",
        ),
    }
    if paths.velocitysignal_per_beat is not None:
        metrics[paths.velocitysignal_per_beat] = metric_value(
            analysis["retinal_artery_velocity_signal_filtered_perbeat"],
            unit="mm/s",
        )
    if paths.velocitysignal_filtered is not None:
        metrics[paths.velocitysignal_filtered] = metric_value(
            analysis["retinal_artery_velocity_signal_filtered"],
            unit="mm/s",
        )
    return metrics


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
