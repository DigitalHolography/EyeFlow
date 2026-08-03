"""Region-indexed velocity and per-beat velocity outputs."""

from __future__ import annotations

import numpy as np

from calculations.math import butter_lowpass_filtfilt, nanmedian
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue, with_attrs
from pipelines.waveform_velocity_core.dopplerview.constants import (
    LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
)
from pipelines.waveform_velocity_core.regions import (
    REGION_NAMES,
    normalize_spatial_frame,
    optic_disc_center_xy,
    region_membership,
)


def pack_hemifield_velocity_outputs(
    metrics: dict[str, object],
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack eight-region continuous and per-beat velocity signals."""
    schema = _resolve_output_paths(output_paths)
    result: dict[str, object] = {}

    for vessel_name, segments, per_beat_paths in (
        ("artery", artery_segments, schema.artery_per_beat),
        ("vein", vein_segments, schema.vein_per_beat),
    ):
        if segments is None or np.asarray(segments.branch_ids).size == 0:
            continue

        branch_ids = np.asarray(segments.branch_ids, dtype=np.int32).reshape(-1)
        labels = np.asarray(segments.labels, dtype=np.int32)
        centers = np.asarray(segments.segment_center_xy, dtype=float)
        center_xy = optic_disc_center_xy(source_data, labels.shape)
        labels, center_xy = normalize_spatial_frame(labels, center_xy)
        membership = region_membership(branch_ids, labels, centers, center_xy)
        result.update(
            _pack_region_velocity_outputs(
                schema,
                source_data,
                vessel_name,
                per_beat_paths,
                segments,
                membership,
                metrics,
            )
        )

    return result


def _pack_region_velocity_outputs(
    schema: EyeFlowOutputPaths,
    source_data,
    vessel_name: str,
    per_beat_paths,
    segments,
    membership: np.ndarray,
    metrics: dict[str, object],
) -> dict[str, object]:
    segment_velocity = np.asarray(segments.velocity, dtype=np.float32)
    if segment_velocity.ndim != 3:
        raise ValueError(
            "Segment velocity must have shape (radius, branch, frame), got "
            f"{segment_velocity.shape}."
        )
    if segment_velocity.shape[:2] != membership.shape[1:][::-1]:
        raise ValueError(
            "Segment velocity and hemifield membership dimensions do not match: "
            f"{segment_velocity.shape[:2]} versus {membership.shape[1:][::-1]}."
        )

    velocity_signal_path = (
        schema.analysis.retinal_artery_velocity_signal
        if vessel_name == "artery"
        else schema.analysis.retinal_vein_velocity_signal
    )
    velocity_filtered_path = (
        schema.analysis.retinal_artery_velocity_signal_filtered
        if vessel_name == "artery"
        else schema.analysis.retinal_vein_velocity_signal_filtered
    )
    velocity_root = _dataset_group(velocity_signal_path)
    output: dict[str, object] = {}
    for region_index, region_name in enumerate(REGION_NAMES):
        selected = membership[region_index]
        raw = _reduce_segment_velocity(segment_velocity, selected)
        filtered = _lowpass_velocity(raw, source_data)
        root = f"{velocity_root}/hemifield/{region_name}"
        output[f"{root}/{_path_variant(velocity_signal_path)}/value"] = with_attrs(
            raw,
            _velocity_region_attrs(region_name, "raw", ("frame",)),
        )
        output[f"{root}/{_path_variant(velocity_filtered_path)}/value"] = with_attrs(
            filtered,
            _velocity_region_attrs(region_name, "filtered", ("frame",)),
        )

    output.update(
        _pack_region_per_beat_velocity_outputs(
            per_beat_paths,
            metrics,
            membership,
            vessel_name,
        )
    )
    return output


def _pack_region_per_beat_velocity_outputs(
    paths,
    metrics: dict[str, object],
    membership: np.ndarray,
    vessel_name: str,
) -> dict[str, object]:
    if paths.segment_velocity_signal is None:
        return {}
    if paths.segment_velocity_signal_band_limited is None:
        return {}

    raw = _optional_payload(metrics, paths.segment_velocity_signal)
    band_limited = _optional_payload(
        metrics,
        paths.segment_velocity_signal_band_limited,
    )
    if (raw is None) != (band_limited is None):
        raise ValueError(
            f"Incomplete {vessel_name} hemifield per-beat segment outputs."
        )
    if raw is None:
        return {}

    raw = np.asarray(raw, dtype=np.float32)
    band_limited = np.asarray(band_limited, dtype=np.float32)
    expected_shape = tuple(membership.shape[1:])
    for name, values in (("raw", raw), ("band-limited", band_limited)):
        if values.ndim != 4 or values.shape[2:] != expected_shape:
            raise ValueError(
                f"Per-beat {name} segment velocity must have shape "
                f"(sample, beat, branch, radius), got {values.shape}; "
                f"expected branch/radius {expected_shape}."
            )

    root = _dataset_group(paths.velocity_signal)
    result: dict[str, object] = {}
    for region_index, region_name in enumerate(REGION_NAMES):
        selected = membership[region_index]
        region_raw = _reduce_per_beat_segment_velocity(raw, selected)
        region_band_limited = _reduce_per_beat_segment_velocity(
            band_limited,
            selected,
        )
        region_root = f"{root}/hemifield/{region_name}"
        result[
            f"{region_root}/{_path_variant(paths.velocity_signal)}/value"
        ] = with_attrs(
            region_raw,
            _velocity_region_attrs(region_name, "raw", ("beat", "sample")),
        )
        result[
            f"{region_root}/{_path_variant(paths.velocity_signal_band_limited)}/value"
        ] = with_attrs(
            region_band_limited,
            _velocity_region_attrs(
                region_name,
                "bandlimited",
                ("beat", "sample"),
            ),
        )
    return result


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


def _lowpass_velocity(values: np.ndarray, source_data) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    if not np.any(np.isfinite(values)):
        return np.full(values.shape, np.nan, dtype=np.float32)
    timing = getattr(source_data, "timing", None)
    if timing is None:
        raise ValueError("Hemifield filtered velocity requires source timing.")
    return butter_lowpass_filtfilt(
        values,
        dt_seconds=np.float32(timing.dt_seconds),
        lowpass_freq_hz=np.float32(LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ),
        order=4,
    )


def _velocity_region_attrs(
    region_name: str,
    signal_type: str,
    dim_desc: tuple[str, ...],
) -> dict[str, object]:
    return {
        "aggregation": "median over selected branch-radius segment velocities",
        "dimDesc": list(dim_desc),
        "region": region_name,
        "signal_type": signal_type,
        "unit": "mm/s",
    }


def _dataset_group(path: str) -> str:
    parts = path.split("/")
    if len(parts) >= 2 and parts[-1] == "value":
        return "/".join(parts[:-2])
    return "/".join(parts[:-1])


def _path_variant(path: str) -> str:
    parts = path.split("/")
    if len(parts) >= 2 and parts[-1] == "value":
        return parts[-2]
    return parts[-1]


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


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
