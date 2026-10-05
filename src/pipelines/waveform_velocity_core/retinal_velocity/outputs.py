"""Output packing for shared retinal velocity products."""

from collections.abc import Mapping

from input_output.schema import EyeFlowOutputPaths
from input_output.writers.h5 import metric_data, metric_value


def pack_retinal_velocity_outputs(
    velocity_analysis: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack shared frequency-map, heartbeat, and provenance analysis outputs."""

    schema = EyeFlowOutputPaths.active(output_paths)
    paths = schema.analysis
    metrics = {
        paths.fRMS_avg: metric_value(velocity_analysis["fRMS_avg"]),
        paths.fRMS_bkg_avg: metric_value(velocity_analysis["fRMS_bkg_avg"]),
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
