"""EyeFlow adapter for the AngioEye waveform-shape calculator."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths, VelocityPerBeatOutputPaths

from .calculator import WaveformShapeMetricsCalculator
from .models import VesselWaveformInputs, WaveformShapeMetricInputs


def run_waveform_shape_metrics_angioeye(ctx) -> dict[str, object]:
    """Read EyeFlow outputs, calculate metrics, and prefix output paths."""

    schema = EyeFlowOutputPaths.active()
    _log(ctx, "Starting AngioEye waveform-shape input loading...")
    inputs = WaveformShapeMetricInputs(
        beat_period_seconds=_required_array(ctx, schema.beat_period_seconds),
        artery=_read_vessel(ctx, schema.artery_per_beat),
        vein=_read_vessel(ctx, schema.vein_per_beat),
    )
    _log(ctx, "Starting AngioEye waveform-shape metric calculation...")
    metrics = WaveformShapeMetricsCalculator().compute(inputs)
    return {
        f"{schema.waveform_shape_metrics_root}/{key}": value
        for key, value in metrics.items()
    }


def _read_vessel(ctx, paths: VelocityPerBeatOutputPaths) -> VesselWaveformInputs:
    raw_global, bandlimited_global = _optional_pair(
        ctx,
        paths.velocity_signal,
        paths.velocity_signal_band_limited,
        "global waveform",
    )
    raw_segments, bandlimited_segments = _optional_pair(
        ctx,
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


def _required_array(ctx, path: str) -> np.ndarray:
    return np.asarray(ctx.ef[path], dtype=np.float32)


def _optional_array(ctx, path: str | None) -> np.ndarray | None:
    if path is None or path not in ctx.ef:
        return None
    return np.asarray(ctx.ef[path], dtype=np.float32)


def _optional_pair(
    ctx,
    raw_path: str | None,
    bandlimited_path: str | None,
    label: str,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    raw = _optional_array(ctx, raw_path)
    bandlimited = _optional_array(ctx, bandlimited_path)
    if (raw is None) != (bandlimited is None):
        raise ValueError(f"Incomplete {label} inputs.")
    return raw, bandlimited


def _log(ctx, message: str) -> None:
    log = getattr(ctx, "log", None)
    if callable(log):
        log(message)
