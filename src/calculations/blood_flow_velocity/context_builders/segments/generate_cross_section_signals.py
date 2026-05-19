"""Cross-section velocity signal port from CrossSection/generateCrossSectionSignals.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi

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
    branch_identity: BranchIdentityResult


def generate_cross_section_signals(
    velocity,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
) -> CrossSectionSignalResult:
    velocity_array = np.asarray(velocity, dtype=np.float32)
    vessel = np.asarray(vessel_mask, dtype=bool)
    branches = label_vessel_branches(vessel, optic_disc_center, ring_settings)
    if branches.branch_ids.size == 0:
        return _empty_result(velocity_array, vessel, ring_settings, branches)

    masks = section_masks(vessel.shape, optic_disc_center, ring_settings)
    shape = (ring_settings.ring_count, branches.branch_ids.size, velocity_array.shape[0])
    segment_v = np.full(shape, np.nan, dtype=np.float32)
    segment_safe = np.full(shape, np.nan, dtype=np.float32)
    _fill_cross_section_signals(
        segment_v,
        segment_safe,
        velocity_array,
        masks,
        branches,
        optic_disc_center,
        cross_section_settings,
    )
    return CrossSectionSignalResult(
        segment_v,
        segment_safe,
        branches.labels,
        branches.branch_ids,
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
        branches,
    )


def _fill_cross_section_signals(
    segment_v: np.ndarray,
    segment_safe: np.ndarray,
    velocity: np.ndarray,
    masks: np.ndarray,
    branches: BranchIdentityResult,
    optic_disc_center,
    settings: CrossSectionSignalSettings,
) -> None:
    for circle_index, section in enumerate(masks):
        for branch_index, branch_id in enumerate(branches.branch_ids):
            mask = section & (branches.labels == int(branch_id))
            loc = _centroid_xy(mask)
            if loc is None:
                continue
            tilt = _tilt_angle(mask, optic_disc_center)
            raw, safe = _cross_section_velocity(
                velocity,
                mask,
                loc,
                tilt,
                optic_disc_center,
                settings,
            )
            segment_v[circle_index, branch_index] = raw
            segment_safe[circle_index, branch_index] = safe


def _centroid_xy(mask: np.ndarray) -> tuple[int, int] | None:
    labeled, count = ndi.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    if count == 0:
        return None
    sizes = np.bincount(labeled.reshape(-1))
    label_id = int(np.argmax(sizes[1:]) + 1)
    y, x = ndi.center_of_mass(mask, labeled, label_id)
    return int(np.rint(x)), int(np.rint(y))


def _tilt_angle(mask: np.ndarray, optic_disc_center) -> float:
    ny, nx = mask.shape
    center = optic_disc_center_yx(optic_disc_center, ny, nx)
    radius = _mean_radius(mask, center)
    step = np.float32(1.0 / max(float(np.mean(mask.shape)), 1.0))
    inner = annulus_mask(mask.shape, optic_disc_center, max(radius - step, 0.0), radius)
    outer = annulus_mask(mask.shape, optic_disc_center, radius, radius + step)
    p_in = _centroid_float(mask & inner)
    p_out = _centroid_float(mask & outer)
    if p_in is None or p_out is None:
        return np.nan
    return float(np.degrees(np.arctan2(p_out[1] - p_in[1], p_out[0] - p_in[0])))


def _mean_radius(mask: np.ndarray, center: tuple[float, float]) -> float:
    y, x = np.nonzero(mask)
    if y.size == 0:
        return 0.0
    cy, cx = center
    y_norm = (y.astype(np.float32) - np.float32(cy)) / np.float32(mask.shape[0])
    x_norm = (x.astype(np.float32) - np.float32(cx)) / np.float32(mask.shape[1])
    return float(np.nanmean(np.sqrt(x_norm**2 + y_norm**2), dtype=np.float32))


def _centroid_float(mask: np.ndarray) -> tuple[float, float] | None:
    if not np.any(mask):
        return None
    y, x = ndi.center_of_mass(mask)
    return float(x), float(y)


def _cross_section_velocity(
    velocity: np.ndarray,
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    tilt_angle_mask: float,
    optic_disc_center,
    settings: CrossSectionSignalSettings,
) -> tuple[np.ndarray, np.ndarray]:
    sub_stack, sub_mask = _subimage_stack(velocity, mask, loc_xy, settings)
    if sub_stack.size == 0:
        return _nan_signal(velocity), _nan_signal(velocity)
    mean_image = _nanmean(sub_stack, axis=0)
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


def _subimage_stack(
    velocity: np.ndarray,
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    settings: CrossSectionSignalSettings,
) -> tuple[np.ndarray, np.ndarray]:
    half_width = int(round(0.01 * mask.shape[0] * settings.scale_factor_width / 2.0))
    x, y = loc_xy
    x_range = slice(max(x - half_width, 0), min(x + half_width + 1, mask.shape[1]))
    y_range = slice(max(y - half_width, 0), min(y + half_width + 1, mask.shape[0]))
    sub_stack = velocity[:, y_range, x_range].astype(np.float32, copy=True)
    sub_mask = mask[y_range, x_range]
    sub_stack[:, ~sub_mask] = np.nan
    return sub_stack, sub_mask


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
    rotated = _rotate_with_nan(mean_image, angle)
    rotated_mask = _rotate_mask(sub_mask, angle)
    rotated[~rotated_mask] = np.nan
    return rotated.astype(np.float32, copy=False), float(angle)


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
    return int(np.mod(90 + round(float(np.mod(alpha, 360.0))), 180))


def _projection_score(image: np.ndarray, angle: float) -> float:
    rotated = ndi.rotate(image, angle, reshape=False, order=0, mode="constant", cval=0.0)
    n_rows = rotated.shape[0]
    central = rotated[n_rows // 3 : int(np.ceil(2 * n_rows / 3)), :]
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
    profile = _nanmean(image, axis=0)
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
        rotated = _rotate_with_nan(frame, angle)
        profile = _nanmean(rotated, axis=0)
        raw[frame_index] = _nanmean(profile[c1 : c2 + 1])
        safe[frame_index] = _nanmean(profile)
    return raw, safe


def _nanmean(values: np.ndarray, axis=None):
    array = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(array)
    count = np.sum(finite, axis=axis)
    total = np.sum(np.where(finite, array, 0.0), axis=axis, dtype=np.float32)
    result = np.divide(
        total,
        count,
        out=np.full_like(total, np.nan, dtype=np.float32),
        where=count > 0,
    )
    if isinstance(result, np.ndarray):
        return result.astype(np.float32, copy=False)
    return np.float32(result)


def _rotate_with_nan(image: np.ndarray, angle: float) -> np.ndarray:
    valid = np.isfinite(image)
    fill = np.float32(np.nanmin(image[valid]) - 1.0) if np.any(valid) else np.float32(0.0)
    filled = np.where(valid, image, fill)
    rotated = ndi.rotate(filled, angle, reshape=False, order=1, mode="constant", cval=float(fill))
    mask = ndi.rotate(valid.astype(np.float32), angle, reshape=False, order=1, mode="constant", cval=0.0)
    rotated[mask < 0.5] = np.nan
    return rotated.astype(np.float32, copy=False)


def _rotate_mask(mask: np.ndarray, angle: float) -> np.ndarray:
    return ndi.rotate(
        mask.astype(np.float32),
        angle,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
    ) >= 0.5


def _nan_signal(velocity: np.ndarray) -> np.ndarray:
    return np.full((velocity.shape[0],), np.nan, dtype=np.float32)
