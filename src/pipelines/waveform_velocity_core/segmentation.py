"""Pack segmentation products created by the waveform velocity core."""

from __future__ import annotations

import numpy as np

from calculations.topology.geometry import image_half_diagonal
from input_output.schema import EyeFlowOutputPaths

from .retinal_velocity.outputs import metric_data

OPTIC_DISC_LABEL = -1
REGION_AXIS_LABEL = -2


def pack_segmentation_outputs(
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack masks and enriched branch-label maps below ``Segmentation``.

    The source arrays use EyeFlow's image frame, whose Y direction is inverted
    relative to the lower-left image frame used by the published maps.  All
    maps are therefore flipped vertically and transposed to ``(x, y)`` before
    writing. Quadrant calculations continue to use the original in-memory
    arrays and do not consume these visualization overlays.
    """
    schema = _resolve_output_paths(output_paths)
    image_shape = tuple(int(size) for size in source_data.retinal_artery_mask.shape)
    optic_disc = source_data.optic_disc
    source_disc_mask = optic_disc.mask_for(image_shape)
    topology_disc_mask, topology_disc_radius = _topology_optic_disc_mask(
        optic_disc,
        image_shape,
        artery_segments,
        vein_segments,
    )
    if optic_disc.mask is not None:
        mask_source = "dopplerview_segmentation"
    else:
        mask_source = "reconstructed_from_dopplerview_center_width_height"
    center_xy = np.asarray(optic_disc.center, dtype=np.float32)

    segmentation = schema.segmentation
    metrics = {
        segmentation.optic_disc.mask: _segmentation_value(
            _serialize_spatial_image(source_disc_mask),
            _mask_attrs(mask_source),
        ),
    }
    metrics.update(
        _pack_vessel_segmentation(
            segmentation.artery,
            artery_segments,
            source_data.retinal_artery_mask,
            topology_disc_mask,
            center_xy,
            topology_disc_radius,
        )
    )
    metrics.update(
        _pack_vessel_segmentation(
            segmentation.vein,
            vein_segments,
            source_data.retinal_vein_mask,
            topology_disc_mask,
            center_xy,
            topology_disc_radius,
        )
    )
    return metrics


def _pack_vessel_segmentation(
    paths,
    segments,
    vessel_mask,
    optic_disc_mask: np.ndarray,
    center_xy: np.ndarray,
    topology_disc_radius: int,
) -> dict[str, object]:
    expected_shape = tuple(int(size) for size in vessel_mask.shape)
    labels = (
        np.zeros(expected_shape, dtype=np.int32)
        if segments is None
        else np.asarray(segments.labels, dtype=np.int32)
    )
    if labels.shape != expected_shape:
        raise ValueError(
            f"segment labels must have shape {expected_shape}, got {labels.shape}."
        )

    return {
        paths.mask: _segmentation_value(
            _serialize_spatial_image(np.asarray(vessel_mask, dtype=bool)),
            _mask_attrs("dopplerview_segmentation"),
        ),
        paths.branch_label_map: _segmentation_value(
            _serialize_branch_label_map(labels, optic_disc_mask, center_xy),
            _branch_label_attrs(
                _axis_thickness(labels.shape),
                topology_disc_radius,
            ),
        ),
    }


def _topology_optic_disc_mask(
    optic_disc,
    image_shape: tuple[int, int],
    artery_segments,
    vein_segments,
) -> tuple[np.ndarray, int]:
    fallback_radius = None
    for segments in (artery_segments, vein_segments):
        topology = getattr(segments, "topology", None)
        settings = getattr(topology, "ring_settings", None)
        if settings is None:
            continue
        fallback_radius = (
            float(settings.inner_radius_frac)
            * max(image_half_diagonal(*image_shape), 1.0)
        )
        break
    radius = optic_disc.centered_circle_radius_pixels(
        fallback_radius_pixels=fallback_radius,
    )
    return (
        optic_disc.centered_circle_mask_for(
            image_shape,
            fallback_radius_pixels=fallback_radius,
        ),
        radius,
    )


def _serialize_branch_label_map(
    labels: np.ndarray,
    optic_disc_mask: np.ndarray,
    center_xy: np.ndarray,
) -> np.ndarray:
    normalized_labels = np.flip(labels, axis=0).copy()
    normalized_optic_disc_mask = np.flip(optic_disc_mask, axis=0)
    normalized_center = center_xy.copy()
    normalized_center[1] = labels.shape[0] - 1 - normalized_center[1]

    image = normalized_labels
    image[normalized_optic_disc_mask] = OPTIC_DISC_LABEL
    axis_thickness = _axis_thickness(labels.shape)
    center_x = int(np.floor(normalized_center[0]))
    center_y = int(np.floor(normalized_center[1]))
    x_start = max(0, center_x - axis_thickness // 2)
    x_stop = min(labels.shape[1], x_start + axis_thickness)
    y_start = max(0, center_y - axis_thickness // 2)
    y_stop = min(labels.shape[0], y_start + axis_thickness)
    image[:, x_start:x_stop] = REGION_AXIS_LABEL
    image[y_start:y_stop, :] = REGION_AXIS_LABEL
    return image.T.copy()


def _serialize_spatial_image(image: np.ndarray) -> np.ndarray:
    return np.flip(np.asarray(image), axis=0).T.copy()


def _axis_thickness(image_shape: tuple[int, int]) -> int:
    return max(3, int(round(min(image_shape) / 128.0)))


def _mask_attrs(source: str) -> dict[str, object]:
    return {
        "dimDesc": ["x", "y"],
        "coordinate_system": "image_pixel",
        "image_origin": "lower_left",
        "source": source,
        "y_axis_direction": "increasing_toward_north",
    }


def _branch_label_attrs(
    axis_thickness: int,
    topology_disc_radius: int,
) -> dict[str, object]:
    return {
        "axis_label": REGION_AXIS_LABEL,
        "axis_thickness_pixels": axis_thickness,
        "background_label": 0,
        "branch_labels": "original in-memory branch labels",
        "coordinate_system": "image_pixel",
        "description": (
            "Two-dimensional branch label map with optic-disc and "
            "quadrant-axis overlays"
        ),
        "dimDesc": ["x", "y"],
        "image_origin": "lower_left",
        "optic_disc_label": OPTIC_DISC_LABEL,
        "optic_disc_overlay": "centered_bounding_circle_used_by_topology",
        "optic_disc_overlay_radius_pixels": np.int32(topology_disc_radius),
        "overlay_priority": "quadrant axes, optic disc, vessel branches",
        "y_axis_direction": "increasing_toward_north",
    }


def _segmentation_value(data, attrs: dict[str, object]):
    return metric_data(data), attrs


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
