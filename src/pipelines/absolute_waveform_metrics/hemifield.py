"""Region-indexed absolute waveform metric outputs."""

from __future__ import annotations

import numpy as np

from calculations.math import nanmedian
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue, with_attrs
from pipelines.waveform_velocity_core.regions import (
    REGION_NAMES,
    normalize_spatial_frame,
    optic_disc_center_xy,
    region_membership,
)

from .calculator import AbsoluteWaveformMetricsCalculator


def pack_hemifield_metrics(
    metrics: dict[str, object],
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack regional raw and band-limited absolute metric groups."""
    schema = _resolve_output_paths(output_paths)
    calculator = AbsoluteWaveformMetricsCalculator()
    metric_names = tuple(
        [
            item[0]
            for item in (
                *calculator._scalar_metric_keys(),
                *calculator._array_metric_keys(),
            )
        ]
    )
    metric_meta = calculator._metric_meta()
    result: dict[str, object] = {}

    for vessel_name, segments in (
        ("artery", artery_segments),
        ("vein", vein_segments),
    ):
        if segments is None:
            continue

        segment_metrics = _read_segment_metrics(
            metrics,
            schema,
            vessel_name,
            segments,
            metric_names,
        )
        if segment_metrics is None:
            continue

        branch_ids = np.asarray(segments.branch_ids, dtype=np.int32).reshape(-1)
        labels = np.asarray(segments.labels, dtype=np.int32)
        centers = np.asarray(segments.segment_center_xy, dtype=float)
        center_xy = optic_disc_center_xy(source_data, labels.shape)
        labels, center_xy = normalize_spatial_frame(labels, center_xy)
        membership = region_membership(branch_ids, labels, centers, center_xy)
        result.update(
            _pack_region_metrics(
                schema,
                vessel_name,
                branch_ids,
                membership,
                segment_metrics,
                metric_meta,
            )
        )

    return result


def _read_segment_metrics(
    metrics: dict[str, object],
    schema: EyeFlowOutputPaths,
    vessel_name: str,
    segments,
    metric_names: tuple[str, ...],
) -> dict[str, dict[str, np.ndarray]] | None:
    branch_ids = np.asarray(segments.branch_ids).reshape(-1)
    if branch_ids.size == 0:
        return None

    expected_tail = (
        int(branch_ids.size),
        int(np.asarray(segments.velocity).shape[0]),
    )
    outputs: dict[str, dict[str, np.ndarray]] = {}
    for signal_type, group_name in (
        ("raw", "raw_segment"),
        ("bandlimited", "bandlimited_segment"),
    ):
        values: dict[str, np.ndarray] = {}
        for metric_name in metric_names:
            path = (
                f"{schema.absolute_waveform_metrics_root}/{vessel_name}/"
                f"by_segment/{group_name}/{metric_name}"
            )
            if path not in metrics:
                raise KeyError(f"Missing absolute waveform segment metric '{path}'.")
            value = np.asarray(_payload_data(metrics[path]), dtype=float)
            if value.ndim not in (3, 4) or value.shape[1:3] != expected_tail:
                raise ValueError(
                    f"Absolute waveform segment metric '{path}' has shape "
                    f"{value.shape}, expected (beat, {expected_tail[0]}, "
                    f"{expected_tail[1]}[, harmonic])."
                )
            values[metric_name] = value
        outputs[signal_type] = values
    return outputs


def _pack_region_metrics(
    schema: EyeFlowOutputPaths,
    vessel_name: str,
    branch_ids: np.ndarray,
    membership: np.ndarray,
    segment_metrics: dict[str, dict[str, np.ndarray]],
    metric_meta: dict[str, dict[str, object]],
) -> dict[str, object]:
    output: dict[str, object] = {}
    first_values = next(iter(segment_metrics.values()))
    first_metric = next(iter(first_values.values()))
    n_beats = first_metric.shape[0]
    root = f"{schema.absolute_waveform_metrics_root}/{vessel_name}/hemifield"

    for region_index, region_name in enumerate(REGION_NAMES):
        selected = membership[region_index]
        region_prefix = f"{root}/{region_name}"
        branch_indexes = np.flatnonzero(np.any(selected, axis=1))

        for signal_type, values_by_metric in segment_metrics.items():
            global_prefix = f"{region_prefix}/global/{signal_type}"
            for metric_name, values in values_by_metric.items():
                region_values = _median_selected(
                    values,
                    selected,
                    n_beats,
                )
                output[f"{global_prefix}/{metric_name}"] = with_attrs(
                    np.asarray(region_values, dtype=np.float32),
                    _region_attrs(
                        metric_meta[metric_name],
                        region_name,
                        signal_type,
                    ),
                )

            for branch_index in branch_indexes:
                branch_id = int(branch_ids[branch_index])
                branch_selected = selected[branch_index]
                branch_prefix = (
                    f"{region_prefix}/by_branch/branch_{branch_id}/{signal_type}"
                )
                for metric_name, values in values_by_metric.items():
                    branch_values = _median_branch(
                        values,
                        branch_index,
                        branch_selected,
                    )
                    attrs = _region_attrs(
                        metric_meta[metric_name],
                        region_name,
                        signal_type,
                    )
                    attrs.update({"branch_id": branch_id})
                    attrs["aggregation"] = (
                        "median over this branch's selected radius-segment "
                        "absolute metrics"
                    )
                    output[f"{branch_prefix}/{metric_name}"] = with_attrs(
                        np.asarray(branch_values, dtype=np.float32),
                        attrs,
                    )

    return output


def _median_selected(
    values: np.ndarray,
    selected: np.ndarray,
    n_beats: int,
) -> np.ndarray:
    if not np.any(selected):
        shape = (n_beats, *values.shape[3:])
        return np.full(shape, np.nan, dtype=np.float32)
    return nanmedian(values[:, selected, ...], axis=1)


def _median_branch(
    values: np.ndarray,
    branch_index: int,
    selected: np.ndarray,
) -> np.ndarray:
    return nanmedian(values[:, branch_index, selected, ...], axis=1)


def _region_attrs(
    metric_meta: dict[str, object],
    region_name: str,
    signal_type: str,
) -> dict[str, object]:
    attrs = dict(metric_meta)
    attrs.update(
        {
            "aggregation": "median over selected branch-radius segment metrics",
            "dimDesc": ["beat", *attrs.get("array_axis", [])],
            "region": region_name,
            "signal_type": signal_type,
        }
    )
    return attrs


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
