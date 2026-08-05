"""Cross-section velocity signal port from CrossSection/generateCrossSectionSignals.m."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from scipy import ndimage as ndi

from calculations.math import nanmean_float32, rotate_array_threshold, rotate_image_with_nan
from utils.logger import Logger

from .branch_identity import BranchIdentityResult, label_vessel_branches
from .segment_geometry import (
    SegmentRingSettings,
    annulus_mask,
    optic_disc_center_yx,
    section_masks,
)
from .profile_processing import ProfileData, process_velocity_profiles


@dataclass(frozen=True)
class CrossSectionSignalSettings:
    scale_factor_width: float
    hydrodynamic_diameters: bool
    velocity_profile_threshold: float
    rotate_from_mask: bool
    pixel_size_mm: float


@dataclass(frozen=True, kw_only=True)
class CrossSectionProfileOutputs:
    """Raw, interpolated, centered, and fitted transverse profile outputs."""

    velocity_profiles: np.ndarray
    profile_x_micrometers: np.ndarray
    profile_sample_count: np.ndarray
    profile_rotation_degrees: np.ndarray
    centered_velocity_profiles: np.ndarray
    centered_profile_x_micrometers: np.ndarray
    profile_center_micrometers: np.ndarray
    profile_lumen_edges_micrometers: np.ndarray
    profile_centering_fit_r_squared: np.ndarray
    poiseuille_coefficients: np.ndarray
    poiseuille_origin_micrometers: np.ndarray
    poiseuille_roots_micrometers: np.ndarray
    poiseuille_r_squared: np.ndarray
    poiseuille_profile_spatial_std: np.ndarray
    raw_profile: ProfileData
    interpolated_profile: ProfileData


@dataclass(frozen=True)
class CrossSectionSignalResult(CrossSectionProfileOutputs):
    """Segment waveforms and their transverse profile measurements."""

    velocity: np.ndarray
    safe_velocity: np.ndarray
    labels: np.ndarray
    branch_ids: np.ndarray
    segment_center_xy: np.ndarray
    branch_identity: BranchIdentityResult


@dataclass(frozen=True)
class _CircleTiltGeometry:
    radius_inner: float
    radius_outer: float
    tilt_angle: float


@dataclass
class _CrossSectionBuffers:
    velocity: np.ndarray
    safe_velocity: np.ndarray
    segment_center_xy: np.ndarray
    profile_cells: np.ndarray
    profile_std_cells: np.ndarray
    profile_rotation_degrees: np.ndarray

    @classmethod
    def allocate(
        cls,
        *,
        frame_count: int,
        ring_count: int,
        branch_count: int,
    ) -> _CrossSectionBuffers:
        signal_shape = (ring_count, branch_count, frame_count)
        profile_shape = (ring_count, branch_count)
        return cls(
            velocity=np.full(signal_shape, np.nan, dtype=np.float32),
            safe_velocity=np.full(signal_shape, np.nan, dtype=np.float32),
            segment_center_xy=np.full(
                (branch_count, ring_count, 2),
                np.nan,
                dtype=np.float32,
            ),
            profile_cells=_empty_profile_cells(profile_shape),
            profile_std_cells=_empty_profile_cells(profile_shape),
            profile_rotation_degrees=np.full(
                profile_shape,
                np.nan,
                dtype=np.float32,
            ),
        )


def generate_cross_section_signals(
    velocity,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
) -> CrossSectionSignalResult:
    vessel = np.asarray(vessel_mask, dtype=bool)
    branches = label_vessel_branches(vessel, optic_disc_center, ring_settings)
    if branches.branch_ids.size == 0:
        return _empty_result(velocity, vessel, ring_settings, branches)

    masks = section_masks(vessel.shape, optic_disc_center, ring_settings)
    buffers = _CrossSectionBuffers.allocate(
        frame_count=velocity.shape[0],
        ring_count=ring_settings.ring_count,
        branch_count=branches.branch_ids.size,
    )
    _fill_cross_section_buffers(
        buffers,
        velocity,
        masks,
        branches,
        optic_disc_center,
        cross_section_settings,
    )
    velocity_profiles, _, profile_sample_count = (
        _pack_transverse_profiles(
            buffers.profile_cells,
            frame_count=velocity.shape[0],
        )
    )
    packed_std, _, _ = _pack_transverse_profiles(
        buffers.profile_std_cells,
        frame_count=1,
    )
    poiseuille_profile_spatial_std = packed_std[:, :, 0, :]
    processed_profiles = process_velocity_profiles(
        velocity_profiles,
        pixel_size_mm=cross_section_settings.pixel_size_mm,
        velocity_profile_threshold=(
            cross_section_settings.velocity_profile_threshold
        ),
    )
    return CrossSectionSignalResult(
        velocity=buffers.velocity,
        safe_velocity=buffers.safe_velocity,
        labels=branches.labels,
        branch_ids=branches.branch_ids,
        segment_center_xy=buffers.segment_center_xy,
        branch_identity=branches,
        velocity_profiles=velocity_profiles,
        profile_x_micrometers=processed_profiles.raw_x_micrometers,
        profile_sample_count=profile_sample_count,
        profile_rotation_degrees=buffers.profile_rotation_degrees,
        centered_velocity_profiles=processed_profiles.centered_velocity,
        centered_profile_x_micrometers=(
            processed_profiles.centered_x_micrometers
        ),
        profile_center_micrometers=processed_profiles.center_micrometers,
        profile_lumen_edges_micrometers=(
            processed_profiles.lumen_edges_micrometers
        ),
        profile_centering_fit_r_squared=(
            processed_profiles.centering_fit_r_squared
        ),
        poiseuille_coefficients=processed_profiles.poiseuille_coefficients,
        poiseuille_origin_micrometers=(
            processed_profiles.poiseuille_origin_micrometers
        ),
        poiseuille_roots_micrometers=(
            processed_profiles.poiseuille_roots_micrometers
        ),
        poiseuille_r_squared=processed_profiles.poiseuille_r_squared,
        poiseuille_profile_spatial_std=poiseuille_profile_spatial_std,
        raw_profile=processed_profiles.raw_profile,
        interpolated_profile=processed_profiles.interpolated_profile,
    )


def _empty_result(
    velocity: np.ndarray,
    vessel: np.ndarray,
    settings: SegmentRingSettings,
    branches: BranchIdentityResult,
) -> CrossSectionSignalResult:
    shape = (settings.ring_count, 1, velocity.shape[0])
    empty_profiles = np.full(
        (settings.ring_count, 0, velocity.shape[0], 0),
        np.nan,
        dtype=np.float32,
    )
    processed = process_velocity_profiles(
        empty_profiles,
        pixel_size_mm=0.0,
        velocity_profile_threshold=0.5,
    )
    return CrossSectionSignalResult(
        velocity=np.full(shape, np.nan, dtype=np.float32),
        safe_velocity=np.full(shape, np.nan, dtype=np.float32),
        labels=np.zeros(vessel.shape, dtype=np.int32),
        branch_ids=branches.branch_ids,
        segment_center_xy=np.full(
            (0, settings.ring_count, 2),
            np.nan,
            dtype=np.float32,
        ),
        branch_identity=branches,
        velocity_profiles=empty_profiles,
        profile_x_micrometers=processed.raw_x_micrometers,
        profile_sample_count=np.zeros((settings.ring_count, 0), dtype=np.int32),
        profile_rotation_degrees=np.full(
            (settings.ring_count, 0),
            np.nan,
            dtype=np.float32,
        ),
        centered_velocity_profiles=processed.centered_velocity,
        centered_profile_x_micrometers=processed.centered_x_micrometers,
        profile_center_micrometers=processed.center_micrometers,
        profile_lumen_edges_micrometers=processed.lumen_edges_micrometers,
        profile_centering_fit_r_squared=processed.centering_fit_r_squared,
        poiseuille_coefficients=processed.poiseuille_coefficients,
        poiseuille_origin_micrometers=processed.poiseuille_origin_micrometers,
        poiseuille_roots_micrometers=processed.poiseuille_roots_micrometers,
        poiseuille_r_squared=processed.poiseuille_r_squared,
        poiseuille_profile_spatial_std=np.full(
            (settings.ring_count, 0, 0),
            np.nan,
            dtype=np.float32,
        ),
        raw_profile=processed.raw_profile,
        interpolated_profile=processed.interpolated_profile,
    )


def _fill_cross_section_buffers(
    buffers: _CrossSectionBuffers,
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
            buffers.segment_center_xy[branch_index, circle_index] = loc
            measurement = _cross_section_velocity(
                velocity,
                mask,
                section,
                loc,
                optic_disc_center,
                settings,
            )
            raw, safe = measurement[:2]
            buffers.velocity[circle_index, branch_index] = raw
            buffers.safe_velocity[circle_index, branch_index] = safe
            if len(measurement) < 5:
                continue
            buffers.profile_cells[circle_index, branch_index] = measurement[2]
            buffers.profile_rotation_degrees[circle_index, branch_index] = measurement[3]
            buffers.profile_std_cells[circle_index, branch_index] = measurement[4][None, :]


def _empty_profile_cells(shape: tuple[int, int]) -> np.ndarray:
    cells = np.empty(shape, dtype=object)
    cells.fill(None)
    return cells


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    sub_stack, sub_mask = _subimage_stack2(velocity, mask, loc_xy, settings)
    if sub_stack.size == 0:
        return (
            _nan_signal(velocity),
            _nan_signal(velocity),
            np.full((velocity.shape[0], 0), np.nan, dtype=np.float32),
            np.nan,
            np.empty((0,), dtype=np.float32),
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
    return _profile_measurement(sub_stack, angle, c1, c2, rotated_mean)


def _profile_measurement(
    sub_stack: np.ndarray,
    angle: float,
    c1: int,
    c2: int,
    rotated_mean: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    raw, safe, profiles = _frame_velocities(sub_stack, angle, c1, c2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spatial_std = np.nanstd(rotated_mean, axis=0, ddof=1)
    return (
        raw,
        safe,
        profiles,
        float(angle),
        spatial_std.astype(np.float32, copy=False),
    )


def _subimage_stack(
    velocity: np.ndarray,
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    settings: CrossSectionSignalSettings,
) -> tuple[np.ndarray, np.ndarray]:
    x_start, x_stop, y_start, y_stop = _substack_window_bounds(
        mask.shape,
        loc_xy,
        settings,
    )
    x, y = loc_xy
    x_range = slice(x_start, x_stop)
    y_range = slice(y_start, y_stop)

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
        Logger.log_warning(
            "Cross-section window is too small for the vessel segment: "
            f"{x_range.stop - x_range.start}x{y_range.stop - y_range.start} px "
            f"window for a {segment_width}x{segment_height} px segment at ({x}, {y})."
        )

    sub_stack = velocity[:, y_range, x_range].astype(np.float32, copy=True)
    sub_stack[:, ~sub_mask] = np.nan
    return sub_stack, sub_mask


def _subimage_stack2(
    velocity: np.ndarray,
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    settings: CrossSectionSignalSettings,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract an adaptive square that contains the complete segment mask.

    The segment mask is already restricted to one branch section between the
    relevant annulus boundaries. Its axis-aligned bounding-box diagonal is
    used as a conservative square size: that square can contain the segment
    after any in-plane rotation, with a one-pixel interpolation margin. The
    existing signal-based orientation estimation is intentionally unchanged.

    The square is shifted inward when it reaches an image boundary. Pixels
    outside the annular branch mask remain NaN, just as in ``_subimage_stack``.
    """
    x_start, x_stop, y_start, y_stop = _substack_window_bounds2(
        mask.shape,
        mask,
        loc_xy,
        settings,
    )
    x_range = slice(x_start, x_stop)
    y_range = slice(y_start, y_stop)
    sub_mask = mask[y_range, x_range]

    segment_y, segment_x = np.nonzero(mask)
    if segment_x.size and (
        segment_x.min() < x_start
        or segment_x.max() >= x_stop
        or segment_y.min() < y_start
        or segment_y.max() >= y_stop
    ):
        Logger.log_warning(
            "Adaptive cross-section window cannot contain the vessel segment: "
            f"{x_stop - x_start}x{y_stop - y_start} px window for a "
            f"{int(segment_x.max() - segment_x.min() + 1)}x"
            f"{int(segment_y.max() - segment_y.min() + 1)} px segment."
        )

    sub_stack = velocity[:, y_range, x_range].astype(np.float32, copy=True)
    sub_stack[:, ~sub_mask] = np.nan
    return sub_stack, sub_mask


def _substack_window_bounds(
    image_shape: tuple[int, int],
    loc_xy: tuple[int, int],
    settings: CrossSectionSignalSettings,
) -> tuple[int, int, int, int]:
    """Return the clipped pixel bounds used for one rotating substack.

    The bounds are ``x_start, x_stop, y_start, y_stop`` with exclusive stop
    coordinates, matching the slices used by ``_subimage_stack``.
    Keeping this calculation shared with the debug exporter ensures that the
    displayed boxes describe the actual rotating image area.
    """
    side = _substack_window_side(image_shape[0], settings)
    half_width = side // 2
    x, y = loc_xy
    return (
        max(x - half_width, 0),
        min(x + half_width + 1, image_shape[1]),
        max(y - half_width, 0),
        min(y + half_width + 1, image_shape[0]),
    )


def _substack_window_bounds2(
    image_shape: tuple[int, int],
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    settings: CrossSectionSignalSettings,
) -> tuple[int, int, int, int]:
    segment_y, segment_x = np.nonzero(mask)
    if segment_x.size == 0:
        return _substack_window_bounds(image_shape, loc_xy, settings)

    segment_width = int(segment_x.max() - segment_x.min() + 1)
    segment_height = int(segment_y.max() - segment_y.min() + 1)
    diagonal = int(
        np.ceil(np.hypot(float(segment_width), float(segment_height)))
    )
    side = max(
        _substack_window_side(image_shape[0], settings),
        diagonal + 2,
    )
    if side % 2 == 0:
        side += 1

    center_x = int(np.floor((segment_x.min() + segment_x.max()) / 2.0 + 0.5))
    center_y = int(np.floor((segment_y.min() + segment_y.max()) / 2.0 + 0.5))
    x_start, x_stop = _fit_window_to_axis(center_x, side, image_shape[1])
    y_start, y_stop = _fit_window_to_axis(center_y, side, image_shape[0])
    return x_start, x_stop, y_start, y_stop


def _substack_window_side(
    image_height: int,
    settings: CrossSectionSignalSettings,
) -> int:
    window_width = int(
        np.floor(0.01 * image_height * settings.scale_factor_width + 0.5)
    )
    half_width = int(np.floor(window_width / 2.0 + 0.5))
    return 2 * half_width + 1


def _fit_window_to_axis(
    center: int,
    side: int,
    axis_length: int,
) -> tuple[int, int]:
    if axis_length <= 0:
        return 0, 0
    side = min(side, axis_length)
    start = center - side // 2
    start = min(max(start, 0), axis_length - side)
    return start, start + side


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = np.full((sub_stack.shape[0],), np.nan, dtype=np.float32)
    safe = np.full_like(raw, np.nan)
    profiles = np.full(
        (sub_stack.shape[0], sub_stack.shape[2]),
        np.nan,
        dtype=np.float32,
    )
    for frame_index, frame in enumerate(sub_stack):
        rotated = rotate_image_with_nan(frame, angle)
        profile = nanmean_float32(rotated, axis=0)
        profiles[frame_index] = profile
        value = nanmean_float32(profile[c1 : c2 + 1])
        raw[frame_index] = np.float32(0.0) if np.isnan(value) else value
        safe[frame_index] = nanmean_float32(profile)
    return raw, safe, profiles


def _pack_transverse_profiles(
    profile_cells: np.ndarray,
    *,
    frame_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Center-pad variable-width measured profiles onto one pixel grid.

    The values remain on their measured transverse pixel support. No transverse
    interpolation or extrapolation is performed here; those paper-specific
    operations are intentionally left to the downstream asymmetry analysis.
    """

    widths = np.zeros(profile_cells.shape, dtype=np.int32)
    for index in np.ndindex(profile_cells.shape):
        profile = profile_cells[index]
        if isinstance(profile, np.ndarray) and profile.ndim == 2:
            widths[index] = np.int32(profile.shape[1])
    max_width = int(np.max(widths, initial=0))
    packed = np.full(
        (*profile_cells.shape, frame_count, max_width),
        np.nan,
        dtype=np.float32,
    )
    for index in np.ndindex(profile_cells.shape):
        profile = profile_cells[index]
        if not isinstance(profile, np.ndarray) or profile.ndim != 2:
            continue
        width = int(profile.shape[1])
        start = (max_width - width) // 2
        copy_frames = min(frame_count, int(profile.shape[0]))
        packed[index + (slice(0, copy_frames), slice(start, start + width))] = (
            profile[:copy_frames]
        )
    x_pixels = (
        np.arange(max_width, dtype=np.float32)
        - np.float32((max_width - 1) / 2.0)
    )
    return packed, x_pixels, widths




def _nan_signal(velocity: np.ndarray) -> np.ndarray:
    return np.full((velocity.shape[0],), np.nan, dtype=np.float32)
