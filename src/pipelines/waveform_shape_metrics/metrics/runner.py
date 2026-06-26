"""In-memory adapter for waveform-shape metric calculations."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths, VelocityPerBeatOutputPaths
from pipeline_engine import DatasetValue

from .calculator import WaveformShapeMetricsCalculator
from .models import VesselWaveformInputs, WaveformShapeMetricInputs


def run_waveform_shape_metric_calculations(
    metrics: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Calculate waveform-shape metrics from already-packed EyeFlow outputs."""

    schema = _resolve_output_paths(output_paths)
    inputs = waveform_shape_metric_inputs_from_outputs(metrics, schema)
    computed = WaveformShapeMetricsCalculator().compute(inputs)
    return {
        f"{schema.waveform_shape_metrics_root}/{key}": value
        for key, value in computed.items()
    }


def waveform_shape_metric_inputs_from_outputs(
    metrics: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> WaveformShapeMetricInputs:
    """Build calculator inputs from a metrics dictionary, without HDF5 reads."""

    schema = _resolve_output_paths(output_paths)
    return WaveformShapeMetricInputs(
        beat_period_seconds=_required_array(metrics, schema.beat_period_seconds),
        artery=_read_vessel(metrics, schema.artery_per_beat),
        vein=_read_vessel(metrics, schema.vein_per_beat),
    )


def _read_vessel(
    metrics: Mapping[str, object],
    paths: VelocityPerBeatOutputPaths,
) -> VesselWaveformInputs:
    raw_global, bandlimited_global = _optional_pair(
        metrics,
        paths.velocity_signal,
        paths.velocity_signal_band_limited,
        "global waveform",
    )
    raw_segments, bandlimited_segments = _optional_pair(
        metrics,
        paths.segment_velocity_signal,
        paths.segment_velocity_signal_band_limited,
        "segment waveform",
    )
    return VesselWaveformInputs(
        raw_global=raw_global,
        bandlimited_global=bandlimited_global,
        raw_segments=raw_segments,
        bandlimited_segments=bandlimited_segments,
    )


def _required_array(metrics: Mapping[str, object], path: str) -> np.ndarray:
    try:
        value = metrics[path]
    except KeyError as exc:
        raise KeyError(f"Missing waveform-shape metric input '{path}'.") from exc
    return np.asarray(_payload_data(value), dtype=np.float32)


def _optional_array(
    metrics: Mapping[str, object],
    path: str | None,
) -> np.ndarray | None:
    if path is None or path not in metrics:
        return None
    return np.asarray(_payload_data(metrics[path]), dtype=np.float32)


def _optional_pair(
    metrics: Mapping[str, object],
    raw_path: str | None,
    bandlimited_path: str | None,
    label: str,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    raw = _optional_array(metrics, raw_path)
    bandlimited = _optional_array(metrics, bandlimited_path)
    if (raw is None) != (bandlimited is None):
        raise ValueError(f"Incomplete {label} inputs.")
    return raw, bandlimited


def _payload_data(value: object) -> object:
    if isinstance(value, DatasetValue):
        return value.data
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        return value[0]
    return value


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
