"""Region-indexed velocity and per-beat velocity outputs."""

from __future__ import annotations

import numpy as np

from calculations.math import nanmedian
from calculations.topology import QUADRANT_NAMES, quadrant_membership
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue, with_attrs
from pipelines.velocity.models import RetinalVelocity
from pipelines.velocity.semantics import velocity_dataset_attrs

from ..analysis.filtering import lowpass_velocity_signals
from .paths import resolve_output_paths

QUADRANTS_GROUP_NAME = "Quadrants"


def pack_quadrant_velocity_outputs(
    metrics: dict[str, object],
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    velocity: RetinalVelocity | None = None,
) -> dict[str, object]:
    """Pack quadrant-level continuous and per-beat velocity signals."""
    schema = resolve_output_paths(output_paths)
    velocity_attrs = velocity_dataset_attrs(velocity)
    result: dict[str, object] = {}

    for vessel_name, segments, per_beat_paths in (
        ("artery", artery_segments, schema.artery_per_beat),
        ("vein", vein_segments, schema.vein_per_beat),
    ):
        if segments is None or np.asarray(segments.profile.topology.native.branch_ids).size == 0:
            continue

        membership = quadrant_membership(segments.profile.topology.native)
        result.update(
            _pack_region_velocity_outputs(
                schema,
                source_data,
                vessel_name,
                per_beat_paths,
                segments,
                membership,
                metrics,
                velocity_attrs,
            )
        )

    return result


def _quadrant_path(path: str, region: str) -> str:
    parent, _, signal = path.rpartition("/")
    return f"{parent}/Quadrants/{region}/{signal}"


def _pack_region_velocity_outputs(
    schema,
    source_data,
    vessel_name,
    per_beat_paths,
    segments,
    membership,
    metrics,
    velocity_attrs,
) -> dict[str, object]:
    values = np.asarray(segments.profile.segment_signal, dtype=np.float32)
    if values.ndim != 3 or values.shape[:2] != membership.shape[1:][::-1]:
        raise ValueError("Segment velocity must have matching (radius, branch, frame) dimensions.")
    filtered = getattr(segments, "filtered_signal", None)
    if filtered is None:
        filtered = lowpass_velocity_signals(
            values, dt_seconds=float(source_data.source.holodoppler.timing.dt_seconds)
        )
    raw_path = getattr(schema.analysis, f"retinal_{vessel_name}_velocity_signal")
    filtered_path = getattr(schema.analysis, f"retinal_{vessel_name}_velocity_signal_filtered")
    output = {}
    for region_index, region_name in enumerate(QUADRANT_NAMES):
        selected = membership[region_index]
        for path, signals, signal_type in (
            (raw_path, values, "raw"),
            (filtered_path, filtered, "filtered"),
        ):
            attrs = _velocity_region_attrs(region_name, signal_type, ("frame",), velocity_attrs)
            if signal_type == "filtered":
                attrs.update(
                    filter_method="butterworth_filtfilt", filter_order=4, filter_cutoff_hz=15.0
                )
            output[_quadrant_path(path, region_name)] = with_attrs(
                _reduce_segment_velocity(signals, selected), attrs
            )

    for paths in (per_beat_paths, getattr(schema, f"{vessel_name}_per_beat_filtered")):
        for segment_path, global_path in (
            (paths.segment_velocity_signal, paths.velocity_signal),
            (paths.segment_velocity_signal_band_limited, paths.velocity_signal_band_limited),
        ):
            payload = _optional_payload(metrics, segment_path)
            if payload is None:
                continue
            data = np.asarray(payload, dtype=np.float32)
            if data.ndim != 4 or data.shape[2:] != membership.shape[1:]:
                raise ValueError(
                    "Per-beat segment velocity must have matching (sample, beat, branch, radius) dimensions."
                )
            original = metrics[segment_path]
            attrs = dict(
                original.attrs
                if isinstance(original, DatasetValue)
                else original[1]
                if isinstance(original, tuple)
                else {}
            )
            signal_type = (
                "bandlimited"
                if "/BandLimited/" in global_path
                else "filtered"
                if "/Filtered/" in global_path
                else "raw"
            )
            for region_index, region_name in enumerate(QUADRANT_NAMES):
                output[_quadrant_path(global_path, region_name)] = with_attrs(
                    _reduce_per_beat_segment_velocity(data, membership[region_index]),
                    {
                        **attrs,
                        **_velocity_region_attrs(
                            region_name, signal_type, ("beat", "sample"), velocity_attrs
                        ),
                    },
                )
    return output


def _reduce_segment_velocity(
    values: np.ndarray,
    selected: np.ndarray,
) -> np.ndarray:
    frame_by_segment = np.asarray(values, dtype=np.float32).transpose(2, 1, 0)
    if not np.any(selected):
        return np.full(frame_by_segment.shape[0], np.nan, dtype=np.float32)
    return np.asarray(
        nanmedian(frame_by_segment[:, selected], axis=1),
        dtype=np.float32,
    )


def _reduce_per_beat_segment_velocity(
    values: np.ndarray,
    selected: np.ndarray,
) -> np.ndarray:
    if not np.any(selected):
        return np.full(
            (values.shape[1], values.shape[0]),
            np.nan,
            dtype=np.float32,
        )
    selected_indexes = np.flatnonzero(selected.reshape(-1))
    sample_by_beat = nanmedian(
        np.take(
            values.reshape(values.shape[0], values.shape[1], -1),
            selected_indexes,
            axis=2,
        ),
        axis=2,
    )
    return np.asarray(sample_by_beat.T, dtype=np.float32)


def _velocity_region_attrs(
    region_name: str,
    signal_type: str,
    dim_desc: tuple[str, ...],
    velocity_attrs: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        **(velocity_attrs or {"unit": "mm/s"}),
        "aggregation": "median over selected branch-radius segment velocities",
        "dimDesc": list(dim_desc),
        "region": region_name,
        "signal_type": signal_type,
    }


def _optional_payload(metrics: dict[str, object], path: str | None) -> object | None:
    if path is None or path not in metrics:
        return None
    return _payload_data(metrics[path])


def _payload_data(value: object) -> object:
    if isinstance(value, DatasetValue):
        return value.data
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        return value[0]
    return value
