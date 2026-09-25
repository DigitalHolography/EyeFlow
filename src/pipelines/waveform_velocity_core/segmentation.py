"""Pack generic vessel topology and segmentation products."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from calculations.topology import retinal_pixel_size_mm, segment_mask_areas_pixels
from calculations.topology.geometry import AnnulusGeometry, image_half_diagonal
from input_output.schema import EyeFlowOutputPaths

from .retinal_velocity.outputs import metric_data

BACKGROUND_LABEL = -1
ANNULUS_OUTLINE_LABEL = -2
INNER_R0_VESSEL_LABEL = -3
_FOUR_CONNECTED = ndi.generate_binary_structure(2, 1)


def pack_segmentation_outputs(
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack masks and branch/segment label maps below ``Segmentation``.

    The source arrays use EyeFlow's image frame, whose Y direction is inverted
    relative to the lower-left image frame used by the published maps.  All
    maps are flipped vertically and transposed to ``(x, y)`` before writing.
    Published branch IDs are contiguous and zero-based; negative values are
    reserved for the background, annulus outlines, and vessel pixels in R0.
    """
    return _pack_segmentation_outputs(
        source_data.retinal_artery_mask,
        source_data.retinal_vein_mask,
        source_data.optic_disc,
        artery_segments,
        vein_segments,
        output_paths,
    )


def _pack_segmentation_outputs(
    artery_mask,
    vein_mask,
    optic_disc,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None,
) -> dict[str, object]:
    schema = _resolve_output_paths(output_paths)
    artery_mask = np.asarray(artery_mask, dtype=bool)
    vein_mask = np.asarray(vein_mask, dtype=bool)
    image_shape = tuple(int(size) for size in artery_mask.shape)
    source_disc_mask = optic_disc.mask_for(image_shape)
    topology_disc_mask, topology_disc_radius, ring_settings = _topology_geometry(
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
    published_center_xy = np.asarray(
        [center_xy[0], np.float32(image_shape[0] - 1) - center_xy[1]],
        dtype=np.float32,
    )
    optic_disc_width = np.float32(
        np.nan if optic_disc.width is None else optic_disc.width
    )
    optic_disc_height = np.float32(
        np.nan if optic_disc.height is None else optic_disc.height
    )
    pixel_pitch_m = np.float32(retinal_pixel_size_mm(optic_disc) * 1e-3)

    segmentation = schema.segmentation
    metrics = {
        segmentation.optic_disc.mask: _segmentation_value(
            _serialize_spatial_image(source_disc_mask),
            _mask_attrs(mask_source),
        ),
        segmentation.optic_disc.height: _segmentation_value(
            optic_disc_height,
            _optic_disc_dimension_attrs("height"),
        ),
        segmentation.optic_disc.width: _segmentation_value(
            optic_disc_width,
            _optic_disc_dimension_attrs("width"),
        ),
        segmentation.optic_disc.center: _segmentation_value(
            published_center_xy,
            {
                "unit": "pixels",
                "dimDesc": ["coordinate"],
                "coordinates": ["x", "y"],
                "coordinate_system": "image_pixel",
                "image_origin": "lower_left",
                "y_axis_direction": "increasing_toward_north",
            },
        ),
        segmentation.pixel_pitch_m: _segmentation_value(
            pixel_pitch_m,
            {
                "unit": "m",
                "definition": "native retinal pixel pitch",
            },
        ),
    }
    metrics.update(
        _pack_vessel_segmentation(
            segmentation.artery,
            artery_segments,
            artery_mask,
            topology_disc_mask,
            center_xy,
            topology_disc_radius,
            ring_settings,
        )
    )
    metrics.update(
        _pack_vessel_segmentation(
            segmentation.vein,
            vein_segments,
            vein_mask,
            topology_disc_mask,
            center_xy,
            topology_disc_radius,
            ring_settings,
        )
    )
    return metrics


def pack_topology_outputs(
    artery_mask,
    vein_mask,
    optic_disc,
    prepared_topologies,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack segmentation products directly from shared prepared topology."""

    return _pack_segmentation_outputs(
        artery_mask,
        vein_mask,
        optic_disc,
        prepared_topologies.get("artery"),
        prepared_topologies.get("vein"),
        output_paths,
    )


def _pack_vessel_segmentation(
    paths,
    segments,
    vessel_mask,
    optic_disc_mask: np.ndarray,
    center_xy: np.ndarray,
    topology_disc_radius: int,
    ring_settings: AnnulusGeometry,
) -> dict[str, object]:
    expected_shape = tuple(int(size) for size in vessel_mask.shape)
    topology = _segment_topology(segments)
    labels = (
        np.zeros(expected_shape, dtype=np.int32)
        if topology is None
        else np.asarray(topology.labels, dtype=np.int32)
    )
    if labels.shape != expected_shape:
        raise ValueError(
            f"segment labels must have shape {expected_shape}, got {labels.shape}."
        )

    vessel = np.asarray(vessel_mask, dtype=bool)
    branch_map = _base_branch_label_map(labels, vessel, optic_disc_mask)
    r0_outline = _circle_outline(labels.shape, center_xy, topology_disc_radius)
    all_outlines = _annulus_outlines(
        labels.shape,
        center_xy,
        topology_disc_radius,
        ring_settings,
    )

    outputs = {
        paths.mask: _segmentation_value(
            _serialize_spatial_image(vessel),
            _mask_attrs("dopplerview_segmentation"),
        ),
        paths.branch_label_map: _segmentation_value(
            _serialize_label_map(_with_outlines(branch_map, r0_outline)),
            _label_map_attrs("innermost R0 outline", topology_disc_radius),
        ),
        paths.segment_map: _segmentation_value(
            _serialize_label_map(_with_outlines(branch_map, all_outlines)),
            _label_map_attrs("all calculated annulus outlines", topology_disc_radius),
        ),
    }
    if topology is not None and all(
        hasattr(topology, field)
        for field in ("branch_ids", "annulus_masks")
    ):
        segment_mask_area = segment_mask_areas_pixels(topology)
        outputs[paths.segment_mask_area] = _segmentation_value(
            segment_mask_area,
            {
                "unit": "pixels^2",
                "dimDesc": ["branch", "radius"],
                "definition": (
                    "count of native vessel-mask pixels belonging to each "
                    "branch inside each annular section"
                ),
                "branch_ids": np.asarray(topology.branch_ids, dtype=np.int32),
                "annulus_geometry": "native_pixel_center_section_mask",
                "annulus_pixel_coverage": "binary_pixel_center_membership",
            },
        )
        annulus_delta_radius = _topology_delta_radius(
            topology,
            segment_mask_area.shape[1],
        )
        delta_radius = np.full(
            segment_mask_area.shape[1],
            np.nan,
            dtype=np.float32,
        )
        lumen_diameter = np.full(
            segment_mask_area.shape,
            np.nan,
            dtype=np.float32,
        )
        valid_radius = (
            np.any(segment_mask_area > 0, axis=0)
            & np.isfinite(annulus_delta_radius)
            & (annulus_delta_radius != 0)
        )
        np.copyto(
            delta_radius,
            annulus_delta_radius,
            where=valid_radius,
        )
        valid_segment = (segment_mask_area > 0) & valid_radius[None, :]
        np.divide(
            segment_mask_area.astype(np.float32),
            delta_radius[None, :],
            out=lumen_diameter,
            where=valid_segment,
        )
        outputs[paths.lumen_diameter] = _segmentation_value(
            lumen_diameter,
            {
                "unit": "pixels",
                "dimDesc": ["branch", "radius"],
                "definition": "segment mask area divided by delta radius",
                "branch_ids": np.asarray(topology.branch_ids, dtype=np.int32),
            },
        )
        outputs[paths.delta_radius] = _segmentation_value(
            delta_radius,
            {
                "unit": "pixels",
                "dimDesc": ["radius"],
                "definition": (
                    "radial width of each annular section; NaN when no branch "
                    "segment or a finite nonzero width is available"
                ),
            },
        )
    return outputs


def _topology_delta_radius(topology, radius_count: int) -> np.ndarray:
    source = getattr(topology, "delta_radius", None)
    if source is None:
        return np.full(radius_count, np.nan, dtype=np.float32)
    delta_radius = np.asarray(source, dtype=np.float32).reshape(-1)
    if delta_radius.shape != (radius_count,):
        raise ValueError(
            "segment topology delta radius must have one value per annular "
            f"section, got {delta_radius.shape} for {radius_count} sections."
        )
    return delta_radius


def _topology_geometry(
    optic_disc,
    image_shape: tuple[int, int],
    artery_segments,
    vein_segments,
) -> tuple[np.ndarray, int, AnnulusGeometry]:
    settings = None
    fallback_radius = None
    for segments in (artery_segments, vein_segments):
        topology = _segment_topology(segments)
        candidate = getattr(topology, "ring_settings", None)
        if candidate is None:
            continue
        settings = candidate
        fallback_radius = (
            float(settings.inner_radius_frac)
            * max(image_half_diagonal(*image_shape), 1.0)
        )
        break
    if settings is None:
        settings = optic_disc.annulus_geometry(image_shape)
    radius = optic_disc.centered_circle_radius_pixels(
        fallback_radius_pixels=fallback_radius,
    )
    return (
        optic_disc.centered_circle_mask_for(
            image_shape,
            fallback_radius_pixels=fallback_radius,
        ),
        radius,
        settings,
    )


def _segment_topology(value):
    if value is None:
        return None
    topology = getattr(value, "topology", value)
    nested = getattr(topology, "topology", None)
    return nested if nested is not None else topology


def _base_branch_label_map(
    labels: np.ndarray,
    vessel_mask: np.ndarray,
    r0_mask: np.ndarray,
) -> np.ndarray:
    image = np.full(labels.shape, BACKGROUND_LABEL, dtype=np.int32)
    branch_ids = np.unique(labels[(labels > 0) & vessel_mask])
    for published_id, branch_id in enumerate(branch_ids):
        image[vessel_mask & (labels == int(branch_id))] = np.int32(published_id)
    image[vessel_mask & r0_mask] = INNER_R0_VESSEL_LABEL
    return image


def _with_outlines(label_map: np.ndarray, outlines: np.ndarray) -> np.ndarray:
    image = label_map.copy()
    image[np.asarray(outlines, dtype=bool)] = ANNULUS_OUTLINE_LABEL
    return image


def _circle_outline(
    image_shape: tuple[int, int],
    center_xy: np.ndarray,
    radius_pixels: float,
) -> np.ndarray:
    center_x, center_y = (float(value) for value in center_xy)
    y, x = np.ogrid[: image_shape[0], : image_shape[1]]
    circle = (
        (y - center_y) ** 2 + (x - center_x) ** 2
        <= float(radius_pixels) ** 2
    )
    return circle & ~ndi.binary_erosion(circle, structure=_FOUR_CONNECTED)


def _annulus_outlines(
    image_shape: tuple[int, int],
    center_xy: np.ndarray,
    r0_radius_pixels: int,
    settings: AnnulusGeometry,
) -> np.ndarray:
    scale = max(image_half_diagonal(*image_shape), 1.0)
    radii = [float(r0_radius_pixels)]
    radii.extend(
        (
            float(settings.inner_radius_frac)
            + (ring_index + 1) * float(settings.ring_width_frac)
        )
        * scale
        for ring_index in range(int(settings.ring_count))
    )
    outlines = np.zeros(image_shape, dtype=bool)
    for radius in dict.fromkeys(radii):
        outlines |= _circle_outline(image_shape, center_xy, radius)
    return outlines


def _serialize_label_map(image: np.ndarray) -> np.ndarray:
    return np.flip(np.asarray(image, dtype=np.int32), axis=0).T.copy()


def _serialize_spatial_image(image: np.ndarray) -> np.ndarray:
    return np.flip(np.asarray(image), axis=0).T.copy()


def _mask_attrs(source: str) -> dict[str, object]:
    return {
        "dimDesc": ["x", "y"],
        "coordinate_system": "image_pixel",
        "image_origin": "lower_left",
        "source": source,
        "y_axis_direction": "increasing_toward_north",
    }


def _optic_disc_dimension_attrs(dimension: str) -> dict[str, object]:
    return {
        "unit": "pixels",
        "definition": f"optic-disc {dimension} in the aligned image frame",
    }


def _label_map_attrs(
    outlines: str,
    topology_disc_radius: int,
) -> dict[str, object]:
    return {
        "annulus_outline_label": ANNULUS_OUTLINE_LABEL,
        "annulus_outlines": outlines,
        "background_label": BACKGROUND_LABEL,
        "branch_labels": "contiguous zero-based branch IDs",
        "coordinate_system": "image_pixel",
        "description": (
            "Two-dimensional vessel branch label map with thin annulus outlines"
        ),
        "dimDesc": ["x", "y"],
        "image_origin": "lower_left",
        "inner_r0_vessel_label": INNER_R0_VESSEL_LABEL,
        "r0_radius_pixels": np.int32(topology_disc_radius),
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
