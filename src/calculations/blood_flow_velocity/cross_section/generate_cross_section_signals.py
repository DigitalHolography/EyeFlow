"""Cross-section velocity signal port from CrossSection/generateCrossSectionSignals.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi

from calculations.math import nanmean_float32, rotate_array_threshold, rotate_image_with_nan

from .branch_identity import BranchIdentityResult, label_vessel_branches
from .segment_geometry import (
    SegmentRingSettings,
    annulus_mask,
    optic_disc_center_yx,
    section_masks,
)


@dataclass(frozen=True)
class CrossSectionSignalSettings:
    scale_factor_width: float
    hydrodynamic_diameters: bool
    velocity_profile_threshold: float
    rotate_from_mask: bool
    pixel_size_mm: float


@dataclass(frozen=True)
class CrossSectionSignalResult:
    velocity: np.ndarray
    safe_velocity: np.ndarray
    labels: np.ndarray
    branch_ids: np.ndarray
    segment_center_xy: np.ndarray
    branch_identity: BranchIdentityResult


@dataclass(frozen=True)
class _AdaptiveWindowGeometry:
    segment_width: int
    segment_height: int
    window_width: int
    center_x: int
    center_y: int


@dataclass(frozen=True)
class _CircleTiltGeometry:
    radius_inner: float
    radius_outer: float
    tilt_angle: float


def generate_cross_section_signals(
    velocity,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    log=None,
) -> CrossSectionSignalResult:
    vessel = np.asarray(vessel_mask, dtype=bool)
    branches = label_vessel_branches(vessel, optic_disc_center, ring_settings)
    if branches.branch_ids.size == 0:
        return _empty_result(velocity, vessel, ring_settings, branches)

    masks = section_masks(vessel.shape, optic_disc_center, ring_settings)
    shape = (ring_settings.ring_count, branches.branch_ids.size, velocity.shape[0])
    segment_v = np.full(shape, np.nan, dtype=np.float32)
    segment_safe = np.full(shape, np.nan, dtype=np.float32)
    segment_center_xy = np.full(
        (branches.branch_ids.size, ring_settings.ring_count, 2),
        np.nan,
        dtype=np.float32,
    )
    _fill_cross_section_signals(
        segment_v,
        segment_safe,
        segment_center_xy,
        velocity,
        masks,
        branches,
        optic_disc_center,
        cross_section_settings,
        log,
    )
    return CrossSectionSignalResult(
        segment_v,
        segment_safe,
        branches.labels,
        branches.branch_ids,
        segment_center_xy,
        branches,
    )


def _empty_result(
    velocity: np.ndarray,
    vessel: np.ndarray,
    settings: SegmentRingSettings,
    branches: BranchIdentityResult,
) -> CrossSectionSignalResult:
    shape = (settings.ring_count, 1, velocity.shape[0])
    return CrossSectionSignalResult(
        np.full(shape, np.nan, dtype=np.float32),
        np.full(shape, np.nan, dtype=np.float32),
        np.zeros(vessel.shape, dtype=np.int32),
        branches.branch_ids,
        np.full((0, settings.ring_count, 2), np.nan, dtype=np.float32),
        branches,
    )


def _fill_cross_section_signals(
    segment_v: np.ndarray,
    segment_safe: np.ndarray,
    segment_center_xy: np.ndarray,
    velocity: np.ndarray,
    masks: np.ndarray,
    branches: BranchIdentityResult,
    optic_disc_center,
    settings: CrossSectionSignalSettings,
    log=None,
) -> None:
    for circle_index, section in enumerate(masks):
        for branch_index, branch_id in enumerate(branches.branch_ids):
            mask = section & (branches.labels == int(branch_id))
            loc = _centroid_xy(mask)
            if loc is None:
                continue
            segment_center_xy[branch_index, circle_index] = loc
            raw, safe = _cross_section_velocity(
                velocity,
                mask,
                section,
                loc,
                optic_disc_center,
                settings,
                log,
            )
            segment_v[circle_index, branch_index] = raw
            segment_safe[circle_index, branch_index] = safe


def _centroid_xy(mask: np.ndarray) -> tuple[int, int] | None:
    labeled, count = ndi.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    if count == 0:
        return None
    y, x = ndi.center_of_mass(mask, labeled, 1)
    return int(np.floor(x + 0.5)), int(np.floor(y + 0.5))


def _tilt_angle(
    mask: np.ndarray,
    section: np.ndarray,
    optic_disc_center,
) -> float:
    geometry = _circle_tilt_geometry(mask, section, optic_disc_center)
    return np.nan if geometry is None else geometry.tilt_angle


def _circle_tilt_geometry(
    mask: np.ndarray,
    section: np.ndarray,
    optic_disc_center,
) -> _CircleTiltGeometry | None:
    radii = _normalized_radius_grid(mask.shape, optic_disc_center)
    section_radii = radii[np.asarray(section, dtype=bool)]
    if section_radii.size == 0:
        return None
    radius_inner = float(np.min(section_radii))
    radius_outer = float(np.max(section_radii))
    step = np.float32(1.0 / max(float(np.mean(mask.shape)), 1.0))
    inner = annulus_mask(
        mask.shape,
        optic_disc_center,
        max(radius_inner - float(step), 0.0),
        radius_inner + float(step),
    )
    outer = annulus_mask(
        mask.shape,
        optic_disc_center,
        max(radius_outer - float(step), 0.0),
        radius_outer,
    )
    p_in = _centroid_float(mask & inner)
    p_out = _centroid_float(mask & outer)
    if p_in is None or p_out is None:
        return None
    tilt_angle = float(
        np.degrees(np.arctan2(p_out[1] - p_in[1], p_out[0] - p_in[0]))
    )
    return _CircleTiltGeometry(radius_inner, radius_outer, tilt_angle)


def _normalized_radius_grid(
    shape: tuple[int, int],
    optic_disc_center,
) -> np.ndarray:
    ny, nx = shape
    cy, cx = optic_disc_center_yx(optic_disc_center, ny, nx)
    y = np.linspace(0.0, 1.0, ny, dtype=np.float32)[:, None]
    x = np.linspace(0.0, 1.0, nx, dtype=np.float32)[None, :]
    return np.sqrt(
        (y - np.float32(cy / max(ny, 1))) ** 2
        + (x - np.float32(cx / max(nx, 1))) ** 2,
    )


def _centroid_float(mask: np.ndarray) -> tuple[float, float] | None:
    if not np.any(mask):
        return None
    y, x = ndi.center_of_mass(mask)
    return float(x), float(y)


def _cross_section_velocity(
    velocity: np.ndarray,
    mask: np.ndarray,
    section: np.ndarray,
    loc_xy: tuple[int, int],
    optic_disc_center,
    settings: CrossSectionSignalSettings,
    log=None,
) -> tuple[np.ndarray, np.ndarray]:
    sub_stack, sub_mask = _adaptive_subimage_stack(
        velocity,
        mask,
        log=log,
    )
    if sub_stack.size == 0:
        if callable(log):
            log("[ERROR] Adaptive cross-section failed; using the legacy window.")
        sub_stack, sub_mask = _subimage_stack(velocity, mask, loc_xy, settings, log)
        if sub_stack.size == 0:
            return _nan_signal(velocity), _nan_signal(velocity)
        tilt_angle_mask = _tilt_angle(mask, section, optic_disc_center)
        mean_image = nanmean_float32(sub_stack, axis=0)
        rotated_mean, angle = _rotated_mean_image(
            mean_image,
            sub_mask,
            loc_xy,
            optic_disc_center,
            tilt_angle_mask,
            settings,
        )
        c1, c2 = _cross_section_limits(rotated_mean, settings)
        return _frame_velocities(sub_stack, angle, c1, c2)

    circle = _circle_tilt_geometry(mask, section, optic_disc_center)
    if circle is None:
        if callable(log):
            log(
                "[ERROR] Adaptive cross-section circles do not intersect both "
                "ends of the vessel segment; using legacy rotation."
            )
        tilt_angle_mask = _tilt_angle(mask, section, optic_disc_center)
        mean_image = nanmean_float32(sub_stack, axis=0)
        rotated_mean, angle = _rotated_mean_image(
            mean_image,
            sub_mask,
            loc_xy,
            optic_disc_center,
            tilt_angle_mask,
            settings,
        )
        c1, c2 = _cross_section_limits(rotated_mean, settings)
        return _frame_velocities(sub_stack, angle, c1, c2)

    redress_angle = circle.tilt_angle + 90.0
    mean_image = nanmean_float32(sub_stack, axis=0)
    rotated_mean = _rotate_masked_image(mean_image, sub_mask, redress_angle)
    c1, c2 = _cross_section_limits(rotated_mean, settings)
    if callable(log):
        _log_circle_tilt(log, circle, redress_angle)
    return _frame_velocities(sub_stack, redress_angle, c1, c2)


def _subimage_stack(
    velocity: np.ndarray,
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    settings: CrossSectionSignalSettings,
    log=None,
) -> tuple[np.ndarray, np.ndarray]:
    window_width = int(
        np.floor(0.01 * mask.shape[0] * settings.scale_factor_width + 0.5)
    )
    half_width = int(np.floor(window_width / 2.0 + 0.5))
    x, y = loc_xy
    x_range = slice(max(x - half_width, 0), min(x + half_width + 1, mask.shape[1]))
    y_range = slice(max(y - half_width, 0), min(y + half_width + 1, mask.shape[0]))

    # Debug log if the sub-mask is too small for the vessel segment
    sub_mask = mask[y_range, x_range]
    if sub_mask.size and (
        np.any(sub_mask[0])
        or np.any(sub_mask[-1])
        or np.any(sub_mask[:, 0])
        or np.any(sub_mask[:, -1])
    ):
        segment_y, segment_x = np.nonzero(mask)
        segment_width = int(segment_x.max() - segment_x.min() + 1)
        segment_height = int(segment_y.max() - segment_y.min() + 1)
        log(
            "[ERROR] Cross-section window is too small for the vessel segment: "
            f"{x_range.stop - x_range.start}x{y_range.stop - y_range.start} px "
            f"window for a {segment_width}x{segment_height} px segment at ({x}, {y})."
        )

    sub_stack = velocity[:, y_range, x_range].astype(np.float32, copy=True)
    sub_stack[:, ~sub_mask] = np.nan
    return sub_stack, sub_mask


def _adaptive_subimage_stack(
    velocity: np.ndarray,
    mask: np.ndarray,
    *,
    log=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Crop a square window sized to contain the complete vessel segment."""

    window = _adaptive_window_geometry(mask)
    if window is None:
        if callable(log):
            log("[ERROR] Adaptive cross-section received an empty vessel segment.")
        return _empty_subimage_stack(velocity, mask)

    sub_stack, sub_mask, pad_width = _padded_adaptive_crop(
        velocity,
        mask,
        window,
    )

    if callable(log):
        _log_adaptive_geometry(log, window)
        if _has_spatial_padding(pad_width):
            log(
                "[ERROR] Adaptive cross-section window extends beyond the image "
                "FOV; missing pixels were padded."
            )

    return sub_stack, sub_mask


def _adaptive_window_geometry(mask: np.ndarray) -> _AdaptiveWindowGeometry | None:
    segment_y, segment_x = np.nonzero(mask)
    if segment_x.size == 0:
        return None

    segment_width = int(segment_x.max() - segment_x.min() + 1)
    segment_height = int(segment_y.max() - segment_y.min() + 1)
    window_width = int(np.ceil(np.hypot(segment_width, segment_height))) + 2
    if window_width % 2 == 0:
        window_width += 1
    center_x = int(np.floor((segment_x.min() + segment_x.max()) / 2.0 + 0.5))
    center_y = int(np.floor((segment_y.min() + segment_y.max()) / 2.0 + 0.5))
    return _AdaptiveWindowGeometry(
        segment_width,
        segment_height,
        window_width,
        center_x,
        center_y,
    )


def _padded_adaptive_crop(
    velocity: np.ndarray,
    mask: np.ndarray,
    geometry: _AdaptiveWindowGeometry,
) -> tuple[np.ndarray, np.ndarray, tuple[tuple[int, int], ...]]:
    half_width = geometry.window_width // 2
    x_start = geometry.center_x - half_width
    x_stop = geometry.center_x + half_width + 1
    y_start = geometry.center_y - half_width
    y_stop = geometry.center_y + half_width + 1
    clipped_x_start, clipped_x_stop = max(x_start, 0), min(x_stop, mask.shape[1])
    clipped_y_start, clipped_y_stop = max(y_start, 0), min(y_stop, mask.shape[0])
    pad_width = (
        (0, 0),
        (clipped_y_start - y_start, y_stop - clipped_y_stop),
        (clipped_x_start - x_start, x_stop - clipped_x_stop),
    )
    sub_stack = velocity[
        :, clipped_y_start:clipped_y_stop, clipped_x_start:clipped_x_stop
    ].astype(np.float32, copy=True)
    sub_mask = mask[
        clipped_y_start:clipped_y_stop, clipped_x_start:clipped_x_stop
    ]
    sub_stack = np.pad(sub_stack, pad_width, constant_values=np.nan)
    sub_mask = np.pad(sub_mask, pad_width[1:], constant_values=False)
    sub_stack[:, ~sub_mask] = np.nan
    return sub_stack, sub_mask, pad_width


def _empty_subimage_stack(
    velocity: np.ndarray,
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    return velocity[:, :0, :0].astype(np.float32), mask[:0, :0]


def _log_adaptive_geometry(
    log,
    window: _AdaptiveWindowGeometry,
) -> None:
    log(
        "[DEBUG] Adaptive cross-section geometry: "
        f"segment={window.segment_width}x{window.segment_height} px, "
        f"window={window.window_width}x{window.window_width} px."
    )


def _log_circle_tilt(
    log,
    circle: _CircleTiltGeometry,
    redress_angle: float,
) -> None:
    log(
        "[DEBUG] Cross-section tilt geometry: "
        f"circles={circle.radius_inner:.4f}/{circle.radius_outer:.4f}, "
        f"tilt={circle.tilt_angle:.2f} deg, redress={redress_angle:.2f} deg."
    )


def _has_spatial_padding(pad_width: tuple[tuple[int, int], ...]) -> bool:
    return any(value > 0 for axis in pad_width[1:] for value in axis)


def _rotated_mean_image(
    mean_image: np.ndarray,
    sub_mask: np.ndarray,
    loc_xy: tuple[int, int],
    optic_disc_center,
    tilt_angle_mask: float,
    settings: CrossSectionSignalSettings,
) -> tuple[np.ndarray, float]:
    if settings.rotate_from_mask and np.isfinite(tilt_angle_mask):
        angle = tilt_angle_mask + 90.0
    else:
        angle = _estimate_orientation(mean_image, loc_xy, optic_disc_center)
    return _rotate_masked_image(mean_image, sub_mask, angle), float(angle)


def _rotate_masked_image(
    image: np.ndarray,
    mask: np.ndarray,
    angle: float,
) -> np.ndarray:
    rotated = rotate_image_with_nan(image, angle)
    rotated_mask = rotate_array_threshold(mask, angle)
    rotated[~rotated_mask] = np.nan
    return rotated.astype(np.float32, copy=False)


def _estimate_orientation(mean_image: np.ndarray, loc_xy, optic_disc_center) -> float:
    cropped = _crop_circle(mean_image)
    cropped = np.nan_to_num(cropped, nan=0.0)
    cropped[cropped < 0] = 0
    beta = _radial_beta_deg(loc_xy, optic_disc_center, mean_image.shape)
    angles = np.mod(np.arange(beta - 90, beta + 91), 180)
    scores = np.asarray([_projection_score(cropped, angle) for angle in angles])
    return float(angles[int(np.nanargmax(scores))])


def _radial_beta_deg(loc_xy, optic_disc_center, shape: tuple[int, int]) -> int:
    cy, cx = optic_disc_center_yx(optic_disc_center, shape[0], shape[1])
    alpha = np.degrees(np.arctan2(loc_xy[1] - cy, loc_xy[0] - cx))
    alpha_360 = float(np.mod(alpha, 360.0))
    matlab_rounded = int(np.floor(alpha_360 + 0.5))
    return int(np.mod(90 + matlab_rounded, 180))


def _projection_score(image: np.ndarray, angle: float) -> float:
    rotated = ndi.rotate(image, angle, reshape=False, order=0, mode="constant", cval=0.0)
    n_rows = rotated.shape[0]
    start = max(int(np.floor(n_rows / 3)) - 1, 0)
    stop = int(np.ceil(2 * n_rows / 3))
    central = rotated[start:stop, :]
    proj = np.sum(central, axis=0, dtype=np.float32)
    total = np.sum(proj, dtype=np.float32)
    return 0.0 if total <= 0 else float(np.max(proj) / total)


def _crop_circle(image: np.ndarray) -> np.ndarray:
    mask = annulus_mask(image.shape, None, 0.0, 0.5)
    cropped = image.copy()
    cropped[~mask] = np.nan
    return cropped


def _cross_section_limits(
    image: np.ndarray,
    settings: CrossSectionSignalSettings,
) -> tuple[int, int]:
    if not settings.hydrodynamic_diameters:
        return 0, max(image.shape[1] - 1, 0)
    profile = nanmean_float32(image, axis=0)
    limits = _hydrodynamic_limits(profile, settings)
    return limits if limits is not None else (0, max(profile.size - 1, 0))


def _hydrodynamic_limits(
    profile: np.ndarray,
    settings: CrossSectionSignalSettings,
) -> tuple[int, int] | None:
    if profile.size == 0 or np.all(~np.isfinite(profile)):
        return None
    threshold = settings.velocity_profile_threshold * float(np.nanmax(profile))
    central = np.flatnonzero(profile > threshold)
    if central.size < 3:
        return None
    center = float(np.mean(central))
    roots = _half_height_roots(profile, central, center, threshold, settings)
    if roots is None:
        return None
    c1 = max(int(np.ceil(center + roots[0] / settings.pixel_size_mm)), 0)
    c2 = min(int(np.floor(center + roots[1] / settings.pixel_size_mm)), profile.size - 1)
    return (c1, c2) if c1 <= c2 else None


def _half_height_roots(
    profile: np.ndarray,
    central: np.ndarray,
    center: float,
    threshold: float,
    settings: CrossSectionSignalSettings,
) -> tuple[float, float] | None:
    x = (central.astype(np.float32) - np.float32(center)) * np.float32(settings.pixel_size_mm)
    coeff = np.polyfit(x.astype(np.float64), profile[central].astype(np.float64), 2)
    p1, p2, p3 = coeff[0], coeff[1], coeff[2] - threshold
    disc = p2**2 - 4.0 * p1 * p3
    if disc < 0 or p1 == 0:
        return None
    roots = np.sort(((-p2 + np.sqrt(disc)) / (2.0 * p1), (-p2 - np.sqrt(disc)) / (2.0 * p1)))
    return float(roots[0]), float(roots[1])


def _frame_velocities(
    sub_stack: np.ndarray,
    angle: float,
    c1: int,
    c2: int,
) -> tuple[np.ndarray, np.ndarray]:
    raw = np.full((sub_stack.shape[0],), np.nan, dtype=np.float32)
    safe = np.full_like(raw, np.nan)
    for frame_index, frame in enumerate(sub_stack):
        rotated = rotate_image_with_nan(frame, angle)
        profile = nanmean_float32(rotated, axis=0)
        value = nanmean_float32(profile[c1 : c2 + 1])
        raw[frame_index] = np.float32(0.0) if np.isnan(value) else value
        safe[frame_index] = nanmean_float32(profile)
    return raw, safe


def _nan_signal(velocity: np.ndarray) -> np.ndarray:
    return np.full((velocity.shape[0],), np.nan, dtype=np.float32)
