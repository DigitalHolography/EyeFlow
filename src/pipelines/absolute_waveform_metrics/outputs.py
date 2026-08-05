"""Output packing for absolute waveform metrics."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths, VelocityPerBeatOutputPaths
from pipeline_engine import DatasetValue

from .calculator import AbsoluteWaveformMetricsCalculator
from .hemifield import pack_hemifield_metrics
from .models import (
    AbsoluteWaveformMetricInputs,
    VesselAbsoluteWaveformInputs,
)


def pack_absolute_waveform_metric_calculations(
    inputs: AbsoluteWaveformMetricInputs,
    *,
    include_segments: bool = True,
) -> dict[str, object]:
    """Calculate and pack unprefixed absolute metric groups."""
    periods = np.asarray(inputs.beat_period_seconds, dtype=float)
    calculator = AbsoluteWaveformMetricsCalculator()
    metrics: dict[str, object] = {}

    for vessel_name, vessel in (
        ("artery", inputs.artery),
        ("vein", inputs.vein),
    ):
        if (
            include_segments
            and vessel.raw_segments is not None
            and vessel.bandlimited_segments is not None
        ):
            calculator._pack_segment_outputs(
                metrics,
                vessel_name,
                vessel.raw_segments,
                vessel.bandlimited_segments,
                periods,
            )
        if vessel.raw_global is not None and vessel.bandlimited_global is not None:
            calculator._pack_global_outputs(
                metrics,
                vessel_name,
                vessel.raw_global,
                vessel.bandlimited_global,
                periods,
            )

    return metrics


def pack_absolute_waveform_outputs(
    metrics: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    source_data=None,
    artery_segments=None,
    vein_segments=None,
    include_per_beat: bool = True,
    include_segments: bool = True,
    include_hemifield: bool = False,
) -> dict[str, object]:
    """Build the selected absolute waveform metric output groups."""
    if not include_per_beat and not include_hemifield:
        return {}

    schema = _resolve_output_paths(output_paths)
    inputs = absolute_waveform_metric_inputs_from_outputs(metrics, schema)
    computed = pack_absolute_waveform_metric_calculations(
        inputs,
        include_segments=include_segments or include_hemifield,
    )
    computed_outputs = {
        f"{schema.absolute_waveform_metrics_root}/{key}": value
        for key, value in computed.items()
    }
    metrics_with_absolute_outputs = dict(metrics)
    metrics_with_absolute_outputs.update(computed_outputs)
    hemifield_outputs = (
        pack_hemifield_metrics(
            metrics_with_absolute_outputs,
            source_data,
            artery_segments,
            vein_segments,
            schema,
        )
        if include_hemifield
        else {}
    )
    return {
        **{
            key: value
            for key, value in computed_outputs.items()
            if include_per_beat
            and (include_segments or "/by_segment/" not in key)
        },
        **hemifield_outputs,
    }


def absolute_waveform_metric_inputs_from_outputs(
    metrics: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> AbsoluteWaveformMetricInputs:
    """Build calculator inputs from packed EyeFlow per-beat outputs."""
    schema = _resolve_output_paths(output_paths)
    return AbsoluteWaveformMetricInputs(
        beat_period_seconds=_required_array(metrics, schema.beat_period_seconds),
        artery=_read_vessel(metrics, schema.artery_per_beat),
        vein=_read_vessel(metrics, schema.vein_per_beat),
    )


def _read_vessel(
    metrics: Mapping[str, object],
    paths: VelocityPerBeatOutputPaths,
) -> VesselAbsoluteWaveformInputs:
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
    return VesselAbsoluteWaveformInputs(
        raw_global=raw_global,
        bandlimited_global=bandlimited_global,
        raw_segments=raw_segments,
        bandlimited_segments=bandlimited_segments,
    )


def _required_array(metrics: Mapping[str, object], path: str) -> np.ndarray:
    try:
        value = metrics[path]
    except KeyError as exc:
        raise KeyError(
            f"Missing absolute waveform metric input '{path}'."
        ) from exc
    return np.asarray(_payload_data(value), dtype=float)


def _optional_array(
    metrics: Mapping[str, object],
    path: str | None,
) -> np.ndarray | None:
    if path is None or path not in metrics:
        return None
    return np.asarray(_payload_data(metrics[path]), dtype=float)


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
