"""Pack the spatial topology used to build segment waveforms."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.analysis_preparation.segments.segment_geometry import (
    optic_disc_center_yx,
)
from input_output.schema import EyeFlowOutputPaths, VesselTopologyOutputPaths

from .outputs import metric_data


def pack_segment_topology_outputs(
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Return H5 values that preserve waveform branch/radius spatial identity."""

    schema = _resolve_output_paths(output_paths)
    image_shape = tuple(int(size) for size in source_data.retinal_artery_mask.shape)
    optic_disc_mask, mask_source = _optic_disc_mask(source_data, image_shape)
    center_xy = _used_optic_disc_center_xy(source_data.optic_disc_center, image_shape)

    metrics = {
        schema.topology.optic_disc.mask: _topology_value(
            optic_disc_mask,
            dim_desc=("y", "x"),
            coordinate_system="image_pixel",
            source=mask_source,
        ),
        schema.topology.optic_disc.center_xy: _topology_value(
            center_xy,
            dim_desc=("coordinate",),
            coordinate_order=("x", "y"),
            coordinate_system="image_pixel",
            unit="pixel",
        ),
    }
    metrics.update(_pack_vessel_topology(schema.topology.artery, artery_segments))
    metrics.update(_pack_vessel_topology(schema.topology.vein, vein_segments))
    return metrics


def _pack_vessel_topology(
    paths: VesselTopologyOutputPaths,
    segments,
) -> dict[str, object]:
    labels = np.asarray(segments.labels, dtype=np.int32)
    branch_ids = np.asarray(segments.branch_ids, dtype=np.int32)
    centers = np.asarray(segments.segment_center_xy, dtype=np.float32)
    velocity = np.asarray(segments.velocity)
    if velocity.ndim != 3:
        raise ValueError(
            "segment velocity must have shape (radius, branch, frame), "
            f"got {velocity.shape}."
        )
    if branch_ids.size and velocity.shape[1] != branch_ids.size:
        raise ValueError(
            f"BranchIds has {branch_ids.size} entries but the segment waveform "
            f"has {velocity.shape[1]} branches."
        )
    _validate_branch_labels(labels, branch_ids)
    expected_shape = (branch_ids.size, segments.velocity.shape[0], 2)
    if centers.shape != expected_shape:
        raise ValueError(
            f"segment_center_xy must have shape {expected_shape}, got {centers.shape}."
        )
    valid_centers = np.all(np.isfinite(centers), axis=-1)
    missing_centers = np.all(np.isnan(centers), axis=-1)
    if not np.all(valid_centers | missing_centers):
        raise ValueError(
            "Each segment center must contain finite [x, y] coordinates or two NaNs."
        )
    return {
        paths.branch_label_map: _topology_value(
            labels,
            dim_desc=("y", "x"),
            coordinate_system="image_pixel",
            background_label=np.int32(0),
        ),
        paths.branch_ids: _topology_value(
            branch_ids,
            dim_desc=("branch",),
            waveform_dim_desc=("sample", "beat", "branch", "radius"),
            waveform_branch_axis=np.int32(2),
        ),
        paths.segment_center_xy: _topology_value(
            centers,
            dim_desc=("branch", "radius", "coordinate"),
            coordinate_order=("x", "y"),
            coordinate_system="image_pixel",
            unit="pixel",
            waveform_dim_desc=("sample", "beat", "branch", "radius"),
            waveform_branch_axis=np.int32(2),
            waveform_radius_axis=np.int32(3),
        ),
    }


def _validate_branch_labels(labels: np.ndarray, branch_ids: np.ndarray) -> None:
    if labels.ndim != 2:
        raise ValueError(f"BranchLabelMap must be 2-D, got shape {labels.shape}.")
    if np.any(labels < 0):
        raise ValueError("BranchLabelMap must use 0 for background and positive labels.")
    if branch_ids.ndim != 1:
        raise ValueError(f"BranchIds must be 1-D, got shape {branch_ids.shape}.")
    if np.any(branch_ids <= 0) or np.unique(branch_ids).size != branch_ids.size:
        raise ValueError("BranchIds must contain unique positive label IDs.")
    label_ids = np.unique(labels[labels > 0])
    if not np.array_equal(np.sort(branch_ids), label_ids):
        raise ValueError(
            "BranchIds must contain exactly the positive IDs in BranchLabelMap."
        )


def _optic_disc_mask(source_data, image_shape: tuple[int, int]):
    if source_data.optic_disc_mask is not None:
        mask = np.asarray(source_data.optic_disc_mask, dtype=bool)
        if mask.shape != image_shape:
            raise ValueError(
                f"optic_disc_mask must have shape {image_shape}, got {mask.shape}."
            )
        return mask, "dopplerview_segmentation"

    mask = _ellipse_mask(
        image_shape,
        source_data.optic_disc_center,
        source_data.optic_disc_width,
        source_data.optic_disc_height,
    )
    if np.any(mask):
        return mask, "reconstructed_from_dopplerview_center_width_height"
    return mask, "unavailable"


def _ellipse_mask(
    image_shape: tuple[int, int],
    optic_disc_center,
    optic_disc_width,
    optic_disc_height,
) -> np.ndarray:
    width = _positive_scalar(optic_disc_width)
    height = _positive_scalar(optic_disc_height)
    if width is None or height is None:
        return np.zeros(image_shape, dtype=bool)

    center_x, center_y = _used_optic_disc_center_xy(optic_disc_center, image_shape)
    y, x = np.indices(image_shape, dtype=np.float32)
    x_radius = np.float32(width / 2.0)
    y_radius = np.float32(height / 2.0)
    return (
        ((x - center_x) / x_radius) ** 2
        + ((y - center_y) / y_radius) ** 2
        <= 1.0
    )


def _used_optic_disc_center_xy(
    optic_disc_center,
    image_shape: tuple[int, int],
) -> np.ndarray:
    center_y, center_x = optic_disc_center_yx(
        optic_disc_center,
        image_shape[0],
        image_shape[1],
    )
    return np.asarray([center_x, center_y], dtype=np.float32)


def _positive_scalar(value) -> float | None:
    if value is None:
        return None
    array = np.asarray(value, dtype=np.float32).reshape(-1)
    if array.size == 0 or not np.isfinite(array[0]) or array[0] <= 0:
        return None
    return float(array[0])


def _topology_value(
    data,
    *,
    dim_desc: tuple[str, ...],
    coordinate_order: tuple[str, ...] | None = None,
    **attrs,
):
    metadata = {"dimDesc": list(dim_desc), **attrs}
    if coordinate_order is not None:
        metadata["coordinate_order"] = list(coordinate_order)
    return metric_data(data), metadata


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
