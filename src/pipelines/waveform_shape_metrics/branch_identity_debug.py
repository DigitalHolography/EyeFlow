"""Export branch identity stages as PNG debug artifacts."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from calculations.blood_flow_velocity.context_builders.segments import (
    BranchIdentityStages,
    SegmentRingSettings,
)
from calculations.blood_flow_velocity.context_builders.segments.segment_geometry import (
    annulus_mask,
)
from input_output.output_manager import OutputManager, OutputType


def export_branch_identity_stage_pngs(
    output_manager: OutputManager,
    stages: BranchIdentityStages,
    prefix: str,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
) -> list[str]:
    paths = []
    for name, image in _stage_images(stages, optic_disc_center, ring_settings):
        path = output_manager.write(
            image,
            OutputType.PNG,
            f"branch_identity/{prefix}_{name}.png",
        )
        paths.append(str(path))
    return paths


def _stage_images(
    stages: BranchIdentityStages,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
):
    return (
        ("01_input", _mask_image(stages.vessel)),
        ("02_mask_section", _mask_image(stages.section)),
        ("03_skeleton", _mask_image(stages.skeleton)),
        ("04_branch_points", _overlay_points(stages.skeleton, stages.branch_points)),
        ("05_cleaned_skeleton", _mask_image(stages.cleaned_skeleton)),
        ("06_marker_labels", _label_image(stages.marker_labels)),
        ("07_distance_topography", _topography_image(stages.distance_topography)),
        ("08_imposed_minima_topography", _topography_image(stages.imposed_minima_topography)),
        ("09_watershed_labels", _label_image(stages.watershed_labels)),
        ("10_annulus_refined_labels", _label_image(stages.annulus_refined_labels)),
        ("11_per_circle_cleaned_labels", _label_image(stages.per_circle_cleaned_labels)),
        (
            "12_per_circle_cleaned_labels_with_rings",
            _labels_with_ring_overlay(
                stages.per_circle_cleaned_labels,
                optic_disc_center,
                ring_settings,
            ),
        ),
    )


def _mask_image(mask: np.ndarray) -> np.ndarray:
    image = np.zeros((*mask.shape, 3), dtype=np.uint8)
    image[np.asarray(mask, dtype=bool)] = (255, 255, 255)
    return image


def _overlay_points(base: np.ndarray, points: np.ndarray) -> np.ndarray:
    image = _mask_image(base)
    image[np.asarray(points, dtype=bool)] = (255, 0, 0)
    return image


def _label_image(labels: np.ndarray) -> np.ndarray:
    labels = np.asarray(labels, dtype=np.int32)
    image = np.zeros((*labels.shape, 3), dtype=np.uint8)
    for label_id in np.unique(labels[labels > 0]):
        image[labels == label_id] = _label_color(int(label_id))
    return image


def _labels_with_ring_overlay(
    labels: np.ndarray,
    optic_disc_center,
    settings: SegmentRingSettings,
) -> np.ndarray:
    image = _label_image(labels)
    image[_ring_boundaries(labels.shape, optic_disc_center, settings)] = (255, 255, 255)
    return image


def _ring_boundaries(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: SegmentRingSettings,
) -> np.ndarray:
    boundaries = np.zeros(image_shape, dtype=bool)
    for ring_index in range(settings.ring_count):
        ring_inner = settings.inner_radius_frac + ring_index * settings.ring_width_frac
        ring = annulus_mask(
            image_shape,
            optic_disc_center,
            ring_inner,
            ring_inner + settings.ring_width_frac,
        )
        boundaries |= ring & ~ndi.binary_erosion(ring)
    return boundaries


def _topography_image(topography: np.ndarray) -> np.ndarray:
    finite = np.isfinite(topography)
    image = np.zeros((*topography.shape, 3), dtype=np.uint8)
    image[np.isposinf(topography)] = (0, 0, 48)
    image[np.isneginf(topography)] = (255, 0, 0)
    if np.any(finite):
        values = topography[finite]
        span = np.max(values) - np.min(values)
        if span <= 0:
            gray = np.full(values.shape, 255, dtype=np.uint8)
        else:
            gray = np.rint((values - np.min(values)) / span * 255).astype(np.uint8)
        image[finite] = np.column_stack((gray, gray, gray))
    return image


def _label_color(label_id: int) -> tuple[int, int, int]:
    hue = (label_id * 0.61803398875) % 1.0
    return _hsv_to_rgb(hue, 0.75, 1.0)


def _hsv_to_rgb(hue: float, saturation: float, value: float) -> tuple[int, int, int]:
    chroma = value * saturation
    x = chroma * (1 - abs((hue * 6) % 2 - 1))
    match int(hue * 6):
        case 0:
            rgb = (chroma, x, 0)
        case 1:
            rgb = (x, chroma, 0)
        case 2:
            rgb = (0, chroma, x)
        case 3:
            rgb = (0, x, chroma)
        case 4:
            rgb = (x, 0, chroma)
        case _:
            rgb = (chroma, 0, x)
    m = value - chroma
    return tuple(int(round((channel + m) * 255)) for channel in rgb)
