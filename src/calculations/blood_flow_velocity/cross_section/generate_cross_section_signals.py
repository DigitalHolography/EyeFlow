"""Cross-section velocity signal port from CrossSection/generateCrossSectionSignals.m."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from scipy import special

from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32, rotate_array_threshold, rotate_image_with_nan
from runtime_limits import cap_parallel_jobs
from .branch_identity import BranchIdentityResult, label_vessel_branches
from .profile_processing import ProfileData, process_velocity_profiles
from .segment_geometry import (
    SegmentRingSettings,
    annulus_mask,
    optic_disc_center_yx,
    section_masks,
)


@dataclass(frozen=True)
class CrossSectionSignalSettings:
    hydrodynamic_diameters: bool
    velocity_profile_threshold: float
    rotate_from_mask: bool
    pixel_size_mm: float
    submask_size_percentile_kept: float = 0.95


@dataclass(frozen=True, kw_only=True)
class CrossSectionProfileOutputs:
    """Transverse and longitudinal cross-section profile outputs."""

    velocity_profiles: np.ndarray
    transverse_velocity_profiles_masked: np.ndarray
    longitudinal_velocity_profiles_unmasked: np.ndarray
    longitudinal_velocity_profiles_masked: np.ndarray
    profile_x_micrometers: np.ndarray
    profile_sample_count: np.ndarray
    profile_rotation_degrees: np.ndarray
    rotated_mean_images: np.ndarray
    rotated_mean_images_masked: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_window_side_pixels: int
    profile_pixel_size_mm: float
    profile_integration_limits_pixels: np.ndarray
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
    velocity_maps_per_segment: np.ndarray
    segment_masks: np.ndarray
    labels: np.ndarray
    branch_ids: np.ndarray
    segment_center_xy: np.ndarray
    branch_identity: BranchIdentityResult

    @property
    def projected_signal(self) -> np.ndarray:
        """Generic alias for ``velocity`` when projecting another cube type."""
        return self.velocity

    @property
    def full_profile_signal(self) -> np.ndarray:
        """Generic alias for the full-width projection signal."""
        return self.safe_velocity

    @property
    def transverse_profiles(self) -> np.ndarray:
        """Generic alias for the uncropped transverse profiles."""
        return self.velocity_profiles


@dataclass(frozen=True)
class _CircleTiltGeometry:
    radius_inner: float
    radius_outer: float
    tilt_angle: float


@dataclass(frozen=True)
class _PreparedCrossSectionGeometry:
    vessel: np.ndarray
    branches: BranchIdentityResult
    masks: np.ndarray


@dataclass(frozen=True)
class _CrossSectionWork:
    circle_index: int
    branch_index: int
    loc_xy: tuple[int, int]
    sub_stack: np.ndarray
    sub_mask: np.ndarray
    bounds_xyxy: tuple[int, int, int, int]
    tilt_angle_mask: float


@dataclass(frozen=True)
class _CrossSectionVelocityMeasurement:
    raw: np.ndarray
    safe_velocity: np.ndarray
    transverse_profiles: np.ndarray
    longitudinal_profiles: np.ndarray
    rotated_stack: np.ndarray
    angle: float
    spatial_std: np.ndarray


@dataclass(frozen=True)
class _CrossSectionMeasurement:
    unmasked: _CrossSectionVelocityMeasurement
    masked: _CrossSectionVelocityMeasurement
    rotated_mean: np.ndarray
    rotated_mean_masked: np.ndarray
    rotated_mask: np.ndarray
    limits: tuple[int, int]
    sample_count: int


_INTERPOLATED_SUBSTACK_SIDE = 128
_ROTATED_SUBSTACK_SIDE = int(_INTERPOLATED_SUBSTACK_SIDE * np.sqrt(2.0))
_MAX_PARALLEL_CROSS_SECTIONS = 8


@dataclass
class _CrossSectionBuffers:
    velocity: np.ndarray
    safe_velocity: np.ndarray
    velocity_maps_per_segment: np.ndarray
    segment_masks: np.ndarray
    segment_center_xy: np.ndarray
    velocity_profiles: np.ndarray
    transverse_velocity_profiles_masked: np.ndarray
    longitudinal_velocity_profiles_unmasked: np.ndarray
    longitudinal_velocity_profiles_masked: np.ndarray
    profile_sample_count: np.ndarray
    profile_spatial_std: np.ndarray
    profile_rotation_degrees: np.ndarray
    rotated_mean_images: np.ndarray
    rotated_mean_images_masked: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_integration_limits_pixels: np.ndarray

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
            velocity_maps_per_segment=np.full(
                (
                    *signal_shape,
                    _ROTATED_SUBSTACK_SIDE,
                    _ROTATED_SUBSTACK_SIDE,
                ),
                np.nan,
                dtype=np.float32,
            ),
            segment_masks=np.zeros(
                (
                    *profile_shape,
                    _ROTATED_SUBSTACK_SIDE,
                    _ROTATED_SUBSTACK_SIDE,
                ),
                dtype=bool,
            ),
            segment_center_xy=np.full(
                (branch_count, ring_count, 2),
                np.nan,
                dtype=np.float32,
            ),
            velocity_profiles=np.full(
                (*signal_shape, _ROTATED_SUBSTACK_SIDE),
                np.nan,
                dtype=np.float32,
            ),
            transverse_velocity_profiles_masked=np.full(
                (*signal_shape, _ROTATED_SUBSTACK_SIDE),
                np.nan,
                dtype=np.float32,
            ),
            longitudinal_velocity_profiles_unmasked=np.full(
                (*signal_shape, _ROTATED_SUBSTACK_SIDE),
                np.nan,
                dtype=np.float32,
            ),
            longitudinal_velocity_profiles_masked=np.full(
                (*signal_shape, _ROTATED_SUBSTACK_SIDE),
                np.nan,
                dtype=np.float32,
            ),
            profile_sample_count=np.zeros(profile_shape, dtype=np.int32),
            profile_spatial_std=np.full(
                (*profile_shape, _ROTATED_SUBSTACK_SIDE),
                np.nan,
                dtype=np.float32,
            ),
            profile_rotation_degrees=np.full(
                profile_shape,
                np.nan,
                dtype=np.float32,
            ),
            rotated_mean_images=np.full(
                (
                    *profile_shape,
                    _ROTATED_SUBSTACK_SIDE,
                    _ROTATED_SUBSTACK_SIDE,
                ),
                np.nan,
                dtype=np.float32,
            ),
            rotated_mean_images_masked=np.full(
                (
                    *profile_shape,
                    _ROTATED_SUBSTACK_SIDE,
                    _ROTATED_SUBSTACK_SIDE,
                ),
                np.nan,
                dtype=np.float32,
            ),
            profile_window_bounds_xyxy=np.full(
                (*profile_shape, 4),
                -1,
                dtype=np.int32,
            ),
            profile_integration_limits_pixels=np.full(
                (*profile_shape, 2),
                -1,
                dtype=np.int32,
            ),
        )


def generate_cross_section_signals(
    velocity,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
) -> CrossSectionSignalResult:
    geometry = _prepare_cross_section_geometry(
        vessel_mask,
        optic_disc_center,
        ring_settings,
    )
    substack_side_pixels = _fixed_substack_side_pixels(
        (geometry,),
        cross_section_settings.submask_size_percentile_kept,
    )
    return _generate_cross_section_signals_from_geometry(
        velocity,
        geometry,
        optic_disc_center,
        ring_settings,
        cross_section_settings,
        substack_side_pixels,
    )


def _prepare_cross_section_geometry(
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
) -> _PreparedCrossSectionGeometry:
    vessel = np.asarray(vessel_mask, dtype=bool)
    branches = label_vessel_branches(vessel, optic_disc_center, ring_settings)
    masks = section_masks(vessel.shape, optic_disc_center, ring_settings)
    return _PreparedCrossSectionGeometry(vessel, branches, masks)


def _generate_cross_section_signals_from_geometry(
    velocity,
    geometry: _PreparedCrossSectionGeometry,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    substack_side_pixels: int,
) -> CrossSectionSignalResult:
    branches = geometry.branches
    if branches.branch_ids.size == 0:
        return _empty_result(
            velocity,
            geometry.vessel,
            ring_settings,
            branches,
            substack_side_pixels=substack_side_pixels,
            profile_pixel_size_mm=_interpolated_pixel_size_mm(
                cross_section_settings.pixel_size_mm,
                substack_side_pixels,
            ),
        )

    buffers = _CrossSectionBuffers.allocate(
        frame_count=velocity.shape[0],
        ring_count=ring_settings.ring_count,
        branch_count=branches.branch_ids.size,
    )
    _fill_cross_section_buffers(
        buffers,
        velocity,
        geometry.masks,
        branches,
        optic_disc_center,
        cross_section_settings,
        substack_side_pixels,
    )
    return _result_from_buffers(
        buffers,
        branches,
        cross_section_settings,
        substack_side_pixels,
    )


def _result_from_buffers(
    buffers: _CrossSectionBuffers,
    branches: BranchIdentityResult,
    settings: CrossSectionSignalSettings,
    substack_side_pixels: int,
) -> CrossSectionSignalResult:
    profile_pixel_size_mm = _interpolated_pixel_size_mm(
        settings.pixel_size_mm,
        substack_side_pixels,
    )
    processed_profiles = process_velocity_profiles(
        buffers.transverse_velocity_profiles_masked,
        pixel_size_mm=profile_pixel_size_mm,
        velocity_profile_threshold=settings.velocity_profile_threshold,
    )
    return CrossSectionSignalResult(
        velocity=buffers.velocity,
        safe_velocity=buffers.safe_velocity,
        velocity_maps_per_segment=buffers.velocity_maps_per_segment,
        segment_masks=buffers.segment_masks,
        labels=branches.labels,
        branch_ids=branches.branch_ids,
        segment_center_xy=buffers.segment_center_xy,
        branch_identity=branches,
        velocity_profiles=buffers.velocity_profiles,
        transverse_velocity_profiles_masked=(
            buffers.transverse_velocity_profiles_masked
        ),
        longitudinal_velocity_profiles_unmasked=(
            buffers.longitudinal_velocity_profiles_unmasked
        ),
        longitudinal_velocity_profiles_masked=(
            buffers.longitudinal_velocity_profiles_masked
        ),
        profile_x_micrometers=processed_profiles.raw_x_micrometers,
        profile_sample_count=buffers.profile_sample_count,
        profile_rotation_degrees=buffers.profile_rotation_degrees,
        rotated_mean_images=buffers.rotated_mean_images,
        rotated_mean_images_masked=buffers.rotated_mean_images_masked,
        profile_window_bounds_xyxy=buffers.profile_window_bounds_xyxy,
        profile_window_side_pixels=int(substack_side_pixels),
        profile_pixel_size_mm=profile_pixel_size_mm,
        profile_integration_limits_pixels=(buffers.profile_integration_limits_pixels),
        centered_velocity_profiles=processed_profiles.centered_velocity,
        centered_profile_x_micrometers=(processed_profiles.centered_x_micrometers),
        profile_center_micrometers=processed_profiles.center_micrometers,
        profile_lumen_edges_micrometers=(processed_profiles.lumen_edges_micrometers),
        profile_centering_fit_r_squared=(processed_profiles.centering_fit_r_squared),
        poiseuille_coefficients=processed_profiles.poiseuille_coefficients,
        poiseuille_origin_micrometers=(processed_profiles.poiseuille_origin_micrometers),
        poiseuille_roots_micrometers=(processed_profiles.poiseuille_roots_micrometers),
        poiseuille_r_squared=processed_profiles.poiseuille_r_squared,
        poiseuille_profile_spatial_std=buffers.profile_spatial_std,
        raw_profile=processed_profiles.raw_profile,
        interpolated_profile=processed_profiles.interpolated_profile,
    )


def _empty_result(
    velocity: np.ndarray,
    vessel: np.ndarray,
    settings: SegmentRingSettings,
    branches: BranchIdentityResult,
    *,
    substack_side_pixels: int,
    profile_pixel_size_mm: float,
) -> CrossSectionSignalResult:
    shape = (settings.ring_count, 1, velocity.shape[0])
    empty_profiles = np.full(
        (
            settings.ring_count,
            0,
            velocity.shape[0],
            _ROTATED_SUBSTACK_SIDE,
        ),
        np.nan,
        dtype=np.float32,
    )
    empty_velocity_maps = np.full(
        (
            settings.ring_count,
            0,
            velocity.shape[0],
            _ROTATED_SUBSTACK_SIDE,
            _ROTATED_SUBSTACK_SIDE,
        ),
        np.nan,
        dtype=np.float32,
    )
    empty_segment_masks = np.zeros(
        (
            settings.ring_count,
            0,
            _ROTATED_SUBSTACK_SIDE,
            _ROTATED_SUBSTACK_SIDE,
        ),
        dtype=bool,
    )
    processed = process_velocity_profiles(
        empty_profiles,
        pixel_size_mm=0.0,
        velocity_profile_threshold=0.5,
    )
    return CrossSectionSignalResult(
        velocity=np.full(shape, np.nan, dtype=np.float32),
        safe_velocity=np.full(shape, np.nan, dtype=np.float32),
        velocity_maps_per_segment=empty_velocity_maps,
        segment_masks=empty_segment_masks,
        labels=np.zeros(vessel.shape, dtype=np.int32),
        branch_ids=branches.branch_ids,
        segment_center_xy=np.full(
            (0, settings.ring_count, 2),
            np.nan,
            dtype=np.float32,
        ),
        branch_identity=branches,
        velocity_profiles=empty_profiles,
        transverse_velocity_profiles_masked=empty_profiles.copy(),
        longitudinal_velocity_profiles_unmasked=empty_profiles.copy(),
        longitudinal_velocity_profiles_masked=empty_profiles.copy(),
        profile_x_micrometers=processed.raw_x_micrometers,
        profile_sample_count=np.zeros((settings.ring_count, 0), dtype=np.int32),
        profile_rotation_degrees=np.full(
            (settings.ring_count, 0),
            np.nan,
            dtype=np.float32,
        ),
        rotated_mean_images=np.full(
            (
                settings.ring_count,
                0,
                _ROTATED_SUBSTACK_SIDE,
                _ROTATED_SUBSTACK_SIDE,
            ),
            np.nan,
            dtype=np.float32,
        ),
        rotated_mean_images_masked=np.full(
            (
                settings.ring_count,
                0,
                _ROTATED_SUBSTACK_SIDE,
                _ROTATED_SUBSTACK_SIDE,
            ),
            np.nan,
            dtype=np.float32,
        ),
        profile_window_bounds_xyxy=np.full(
            (settings.ring_count, 0, 4),
            -1,
            dtype=np.int32,
        ),
        profile_window_side_pixels=int(substack_side_pixels),
        profile_pixel_size_mm=float(profile_pixel_size_mm),
        profile_integration_limits_pixels=np.full(
            (settings.ring_count, 0, 2),
            -1,
            dtype=np.int32,
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
            (settings.ring_count, 0, _ROTATED_SUBSTACK_SIDE),
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
    substack_side_pixels: int,
) -> None:
    work_items: list[_CrossSectionWork] = []
    for circle_index, section in enumerate(masks):
        for branch_index, branch_id in enumerate(branches.branch_ids):
            mask = section & (branches.labels == int(branch_id))
            loc = _centroid_xy(mask)
            if loc is None:
                continue
            buffers.segment_center_xy[branch_index, circle_index] = loc
            sub_stack, sub_mask, bounds_xyxy = _fixed_subimage_stack(
                velocity,
                mask,
                loc,
                substack_side_pixels,
            )
            work_items.append(
                _CrossSectionWork(
                    circle_index=circle_index,
                    branch_index=branch_index,
                    loc_xy=loc,
                    sub_stack=sub_stack,
                    sub_mask=sub_mask,
                    bounds_xyxy=bounds_xyxy,
                    tilt_angle_mask=_tilt_angle(mask, section, optic_disc_center),
                )
            )

    def measure(work: _CrossSectionWork):
        return _cross_section_velocity_from_substack(
            work.sub_stack,
            work.sub_mask,
            work.loc_xy,
            optic_disc_center,
            work.tilt_angle_mask,
            settings,
            substack_side_pixels,
        )

    worker_count = _cross_section_worker_count(len(work_items))
    if worker_count == 1:
        measurements = map(measure, work_items)
    else:
        executor = ThreadPoolExecutor(
            max_workers=worker_count,
            thread_name_prefix="cross-section",
        )
        measurements = executor.map(measure, work_items)
    try:
        for work, measurement in zip(work_items, measurements, strict=True):
            _store_cross_section_measurement(
                buffers,
                work.circle_index,
                work.branch_index,
                measurement,
                work.bounds_xyxy,
            )
    finally:
        if worker_count > 1:
            executor.shutdown(wait=True)


def _cross_section_worker_count(work_count: int) -> int:
    if work_count <= 1 or optional_cupy_backend() is not None:
        return 1
    return min(
        work_count,
        cap_parallel_jobs(_MAX_PARALLEL_CROSS_SECTIONS),
    )


def _store_cross_section_measurement(
    buffers: _CrossSectionBuffers,
    circle_index: int,
    branch_index: int,
    measurement: _CrossSectionMeasurement,
    bounds_xyxy: tuple[int, int, int, int],
) -> None:
    masked = measurement.masked
    buffers.velocity[circle_index, branch_index] = masked.raw
    buffers.safe_velocity[circle_index, branch_index] = masked.safe_velocity
    buffers.velocity_maps_per_segment[circle_index, branch_index] = (
        measurement.unmasked.rotated_stack
    )
    buffers.segment_masks[circle_index, branch_index] = measurement.rotated_mask
    buffers.velocity_profiles[circle_index, branch_index] = (
        measurement.unmasked.transverse_profiles
    )
    buffers.transverse_velocity_profiles_masked[circle_index, branch_index] = (
        masked.transverse_profiles
    )
    buffers.longitudinal_velocity_profiles_unmasked[circle_index, branch_index] = (
        measurement.unmasked.longitudinal_profiles
    )
    buffers.longitudinal_velocity_profiles_masked[circle_index, branch_index] = (
        masked.longitudinal_profiles
    )
    buffers.profile_sample_count[circle_index, branch_index] = measurement.sample_count
    buffers.profile_spatial_std[circle_index, branch_index] = masked.spatial_std
    buffers.profile_rotation_degrees[circle_index, branch_index] = np.float32(
        masked.angle
    )
    buffers.rotated_mean_images[circle_index, branch_index] = measurement.rotated_mean
    buffers.rotated_mean_images_masked[circle_index, branch_index] = (
        measurement.rotated_mean_masked
    )
    buffers.profile_window_bounds_xyxy[circle_index, branch_index] = bounds_xyxy
    buffers.profile_integration_limits_pixels[circle_index, branch_index] = (
        measurement.limits
    )


def _centroid_xy(mask: np.ndarray) -> tuple[int, int] | None:
    labeled, count = ndi.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    if count == 0:
        return None
    y, x = ndi.center_of_mass(mask, labeled, 1)
    return int(np.floor(x + 0.5)), int(np.floor(y + 0.5))


def _fixed_substack_side_pixels(
    geometries: tuple[_PreparedCrossSectionGeometry, ...],
    percentile_kept: float,
) -> int:
    """Return one odd square side from joint segment width/height percentiles."""
    percentile = float(percentile_kept)
    if not 0.0 < percentile <= 1.0:
        raise ValueError("submask_size_percentile_kept must be in (0, 1].")

    widths: list[int] = []
    heights: list[int] = []
    for geometry in geometries:
        for section in geometry.masks:
            for branch_id in geometry.branches.branch_ids:
                mask = section & (geometry.branches.labels == int(branch_id))
                segment_y, segment_x = np.nonzero(mask)
                if segment_x.size == 0:
                    continue
                widths.append(int(segment_x.max() - segment_x.min() + 1))
                heights.append(int(segment_y.max() - segment_y.min() + 1))

    if not widths:
        return 0
    width_kept = int(np.quantile(widths, percentile, method="higher"))
    height_kept = int(np.quantile(heights, percentile, method="higher"))
    side = max(width_kept, height_kept)
    return side if side % 2 == 1 else side + 1


def _interpolated_pixel_size_mm(
    native_pixel_size_mm: float,
    substack_side_pixels: int,
) -> float:
    if substack_side_pixels <= 0:
        return 0.0
    return float(
        native_pixel_size_mm * float(substack_side_pixels) / float(_INTERPOLATED_SUBSTACK_SIDE)
    )


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
    tilt_angle = float(np.degrees(np.arctan2(p_out[1] - p_in[1], p_out[0] - p_in[0])))
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
        (y - np.float32(cy / max(ny, 1))) ** 2 + (x - np.float32(cx / max(nx, 1))) ** 2,
    )


def _centroid_float(mask: np.ndarray) -> tuple[float, float] | None:
    if not np.any(mask):
        return None
    y, x = ndi.center_of_mass(mask)
    return float(x), float(y)


def _cross_section_velocity_from_substack(
    sub_stack: np.ndarray,
    sub_mask: np.ndarray,
    loc_xy: tuple[int, int],
    optic_disc_center,
    tilt_angle_mask: float,
    settings: CrossSectionSignalSettings,
    substack_side_pixels: int,
) -> _CrossSectionMeasurement:
    resized_stack = _resize_subimage_stack(sub_stack)
    resized_mask = _resize_submask(sub_mask)
    resized_stack_masked = resized_stack.copy()
    resized_stack_masked[:, ~resized_mask] = np.nan
    mean_image = nanmean_float32(resized_stack, axis=0)
    mean_image_masked = nanmean_float32(resized_stack_masked, axis=0)
    angle = _mean_image_rotation_angle(
        mean_image_masked,
        loc_xy,
        optic_disc_center,
        tilt_angle_mask,
        settings,
    )
    rotation_stack = _center_pad_for_rotation(resized_stack, np.nan)
    rotation_stack_masked = _center_pad_for_rotation(
        resized_stack_masked,
        np.nan,
    )
    rotation_mask = _center_pad_for_rotation(resized_mask, False)
    rotated_mean = _rotate_mean_image(
        _center_pad_for_rotation(mean_image, np.nan),
        angle,
    )
    rotated_mean_masked = _rotate_masked_image(
        _center_pad_for_rotation(mean_image_masked, np.nan),
        rotation_mask,
        angle,
    )
    rotated_mask = _rotate_mask(rotation_mask, angle)
    profile_pixel_size_mm = _interpolated_pixel_size_mm(
        settings.pixel_size_mm,
        substack_side_pixels,
    )
    c1, c2 = _cross_section_limits(
        rotated_mean_masked,
        settings,
        pixel_size_mm=profile_pixel_size_mm,
    )
    return _CrossSectionMeasurement(
        unmasked=_profile_measurement(
            rotation_stack,
            angle,
            c1,
            c2,
            rotated_mean,
        ),
        masked=_profile_measurement(
            rotation_stack_masked,
            angle,
            c1,
            c2,
            rotated_mean_masked,
        ),
        rotated_mean=rotated_mean,
        rotated_mean_masked=rotated_mean_masked,
        rotated_mask=rotated_mask,
        limits=(c1, c2),
        sample_count=_rotated_profile_sample_count(angle),
    )


def _profile_measurement(
    sub_stack: np.ndarray,
    angle: float,
    c1: int,
    c2: int,
    rotated_mean: np.ndarray,
) -> _CrossSectionVelocityMeasurement:
    (
        raw,
        safe_velocity,
        transverse_profiles,
        longitudinal_profiles,
        rotated_stack,
    ) = _frame_velocities(
        sub_stack,
        angle,
        c1,
        c2,
    )
    spatial_std = _sample_nanstd_axis0(rotated_mean)
    return _CrossSectionVelocityMeasurement(
        raw=raw,
        safe_velocity=safe_velocity,
        transverse_profiles=transverse_profiles,
        longitudinal_profiles=longitudinal_profiles,
        rotated_stack=rotated_stack,
        angle=float(angle),
        spatial_std=spatial_std,
    )


def _sample_nanstd_axis0(values: np.ndarray) -> np.ndarray:
    """Sample standard deviation without warning for columns with <= 1 value."""
    values_float64 = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(values_float64)
    counts = np.sum(finite, axis=0, dtype=np.int32)
    totals = np.sum(
        np.where(finite, values_float64, 0.0),
        axis=0,
        dtype=np.float64,
    )
    means = np.zeros_like(totals)
    np.divide(totals, counts, out=means, where=counts > 0)

    centered = np.where(finite, values_float64 - means, 0.0)
    squared_deviations = np.sum(centered * centered, axis=0, dtype=np.float64)
    variances = np.full_like(totals, np.nan)
    np.divide(
        squared_deviations,
        counts - 1,
        out=variances,
        where=counts > 1,
    )
    return np.sqrt(variances).astype(np.float32, copy=False)


def _fixed_subimage_stack(
    data_cube,
    mask: np.ndarray,
    loc_xy: tuple[int, int],
    side_pixels: int,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int, int]]:
    """Extract one fixed, centroid-centered square with padded periphery."""
    bounds_xyxy = _centered_substack_bounds(
        mask.shape,
        loc_xy,
        side_pixels,
    )
    sub_stack, sub_mask = _subimage_stack_from_bounds(
        data_cube,
        mask,
        bounds_xyxy,
        loc_xy=loc_xy,
        side_pixels=side_pixels,
    )
    return sub_stack, sub_mask, bounds_xyxy


def _subimage_stack_from_bounds(
    data_cube,
    mask: np.ndarray,
    bounds_xyxy: tuple[int, int, int, int],
    *,
    loc_xy: tuple[int, int],
    side_pixels: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a stored fixed window, padding beyond the image boundary."""
    if side_pixels <= 0 or side_pixels % 2 == 0:
        raise ValueError("side_pixels must be a positive odd integer.")
    x_start, x_stop, y_start, y_stop = bounds_xyxy
    x, y = loc_xy
    conceptual_x_start = x - side_pixels // 2
    conceptual_y_start = y - side_pixels // 2
    target_x_start = x_start - conceptual_x_start
    target_y_start = y_start - conceptual_y_start
    target_x_stop = target_x_start + (x_stop - x_start)
    target_y_stop = target_y_start + (y_stop - y_start)

    sub_mask = np.zeros((side_pixels, side_pixels), dtype=bool)
    sub_stack = np.full(
        (data_cube.shape[0], side_pixels, side_pixels),
        np.nan,
        dtype=np.float32,
    )
    if x_start < x_stop and y_start < y_stop:
        source = np.asarray(
            data_cube[:, y_start:y_stop, x_start:x_stop],
            dtype=np.float32,
        )
        source_mask = np.asarray(
            mask[y_start:y_stop, x_start:x_stop],
            dtype=bool,
        )
        target_y = slice(target_y_start, target_y_stop)
        target_x = slice(target_x_start, target_x_stop)
        sub_mask[target_y, target_x] = source_mask
        sub_stack[:, target_y, target_x] = source
    return sub_stack, sub_mask


def _centered_substack_bounds(
    image_shape: tuple[int, int],
    loc_xy: tuple[int, int],
    side_pixels: int,
) -> tuple[int, int, int, int]:
    """Return clipped source bounds for an odd square centered on ``loc_xy``."""
    if side_pixels <= 0 or side_pixels % 2 == 0:
        raise ValueError("side_pixels must be a positive odd integer.")
    half_width = side_pixels // 2
    x, y = loc_xy
    return (
        max(x - half_width, 0),
        min(x + half_width + 1, int(image_shape[1])),
        max(y - half_width, 0),
        min(y + half_width + 1, int(image_shape[0])),
    )


def _resize_subimage_stack(sub_stack: np.ndarray) -> np.ndarray:
    values = np.asarray(sub_stack, dtype=np.float32)
    if optional_cupy_backend() is not None:
        return _resize_values_with_nan(values)
    shared_validity = _shared_validity_mask(values)
    if shared_validity is not None:
        return _resize_stack_with_shared_validity_cpu(values, shared_validity)
    return np.stack(
        [_resize_values_with_nan_cpu(frame) for frame in values],
        axis=0,
    )


def _resize_stack_with_shared_validity_cpu(
    values: np.ndarray,
    shared_validity: np.ndarray,
) -> np.ndarray:
    """Resize a temporal stack while interpolating its shared validity once."""
    zoom_factors = (
        _INTERPOLATED_SUBSTACK_SIDE / values.shape[-2],
        _INTERPOLATED_SUBSTACK_SIDE / values.shape[-1],
    )
    resized_weights = ndi.zoom(
        shared_validity.astype(np.float32),
        zoom_factors,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ).astype(np.float32, copy=False)
    resized_values = np.stack(
        [
            ndi.zoom(
                np.where(shared_validity, frame, np.float32(0.0)),
                zoom_factors,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            ).astype(np.float32, copy=False)
            for frame in values
        ],
        axis=0,
    )
    resized = np.full(resized_values.shape, np.nan, dtype=np.float32)
    np.divide(
        resized_values,
        resized_weights[None, :, :],
        out=resized,
        where=resized_weights[None, :, :] > np.float32(1e-6),
    )
    return resized


def _resize_submask(mask: np.ndarray) -> np.ndarray:
    source = np.asarray(mask, dtype=np.float32)
    zoom_factors = (
        _INTERPOLATED_SUBSTACK_SIDE / source.shape[-2],
        _INTERPOLATED_SUBSTACK_SIDE / source.shape[-1],
    )
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            gpu_mask = backend.cupy.asarray(source)
            resized = backend.ndimage.zoom(
                gpu_mask,
                zoom_factors,
                order=0,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            )
            return backend.cupy.asnumpy(resized) >= np.float32(0.5)
        except Exception:
            pass
    return ndi.zoom(
        source,
        zoom_factors,
        order=0,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ) >= np.float32(0.5)


def _center_pad_for_rotation(
    values: np.ndarray,
    fill_value: float | bool,
) -> np.ndarray:
    """Center a 128-square array on the fixed diagonal-sized canvas."""
    source = np.asarray(values)
    expected_shape = (
        _INTERPOLATED_SUBSTACK_SIDE,
        _INTERPOLATED_SUBSTACK_SIDE,
    )
    if source.ndim < 2 or source.shape[-2:] != expected_shape:
        raise ValueError(
            "rotation input must end with the 128x128 interpolated shape."
        )
    total_padding = _ROTATED_SUBSTACK_SIDE - _INTERPOLATED_SUBSTACK_SIDE
    padding_before = total_padding // 2
    padding_after = total_padding - padding_before
    padding = [(0, 0)] * source.ndim
    padding[-2] = (padding_before, padding_after)
    padding[-1] = (padding_before, padding_after)
    return np.pad(
        source,
        padding,
        mode="constant",
        constant_values=fill_value,
    )


def _rotated_profile_sample_count(angle_degrees: float) -> int:
    """Return the unpadded width of a 128-square rotated by ``angle``."""
    if not np.isfinite(angle_degrees):
        return 0
    radians = np.deg2rad(float(angle_degrees))
    scale = abs(float(np.cos(radians))) + abs(float(np.sin(radians)))
    count = int(np.floor(_INTERPOLATED_SUBSTACK_SIDE * scale + 0.5))
    return min(
        max(count, _INTERPOLATED_SUBSTACK_SIDE),
        _ROTATED_SUBSTACK_SIDE,
    )


def _resize_values_with_nan(values: np.ndarray) -> np.ndarray:
    zoom_factors = (1.0,) * (values.ndim - 2) + (
        _INTERPOLATED_SUBSTACK_SIDE / values.shape[-2],
        _INTERPOLATED_SUBSTACK_SIDE / values.shape[-1],
    )
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            gpu_values = backend.cupy.asarray(values)
            valid = backend.cupy.isfinite(gpu_values)
            resized_values = backend.ndimage.zoom(
                backend.cupy.where(valid, gpu_values, backend.cupy.float32(0.0)),
                zoom_factors,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            )
            resized_weights = backend.ndimage.zoom(
                valid.astype(backend.cupy.float32),
                zoom_factors,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            )
            resized = backend.cupy.full(
                resized_values.shape,
                backend.cupy.nan,
                dtype=backend.cupy.float32,
            )
            backend.cupy.divide(
                resized_values,
                resized_weights,
                out=resized,
                where=resized_weights > backend.cupy.float32(1e-6),
            )
            return backend.cupy.asnumpy(resized)
        except Exception:
            pass

    return _resize_values_with_nan_cpu(values)


def _resize_values_with_nan_cpu(values: np.ndarray) -> np.ndarray:
    zoom_factors = (1.0,) * (values.ndim - 2) + (
        _INTERPOLATED_SUBSTACK_SIDE / values.shape[-2],
        _INTERPOLATED_SUBSTACK_SIDE / values.shape[-1],
    )
    valid = np.isfinite(values)
    resized_values = ndi.zoom(
        np.where(valid, values, np.float32(0.0)),
        zoom_factors,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ).astype(np.float32, copy=False)
    resized_weights = ndi.zoom(
        valid.astype(np.float32),
        zoom_factors,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ).astype(np.float32, copy=False)
    resized = np.full(resized_values.shape, np.nan, dtype=np.float32)
    np.divide(
        resized_values,
        resized_weights,
        out=resized,
        where=resized_weights > np.float32(1e-6),
    )
    return resized


def _mean_image_rotation_angle(
    mean_image_masked: np.ndarray,
    loc_xy: tuple[int, int],
    optic_disc_center,
    tilt_angle_mask: float,
    settings: CrossSectionSignalSettings,
) -> float:
    if settings.rotate_from_mask and np.isfinite(tilt_angle_mask):
        angle = tilt_angle_mask + 90.0
    else:
        angle = _estimate_orientation(mean_image_masked, loc_xy, optic_disc_center)
    return float(angle)


def _rotate_mean_image(image: np.ndarray, angle: float) -> np.ndarray:
    return rotate_image_with_nan(image, angle).astype(np.float32, copy=False)


def _rotate_masked_image(
    image: np.ndarray,
    mask: np.ndarray,
    angle: float,
) -> np.ndarray:
    rotated = _rotate_mean_image(image, angle)
    rotated_mask = _rotate_mask(mask, angle)
    rotated[~rotated_mask] = np.nan
    return rotated.astype(np.float32, copy=False)


def _rotate_mask(mask: np.ndarray, angle: float) -> np.ndarray:
    return rotate_array_threshold(mask, angle).astype(bool, copy=False)


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
    *,
    pixel_size_mm: float,
) -> tuple[int, int]:
    if not settings.hydrodynamic_diameters:
        return 0, max(image.shape[1] - 1, 0)
    profile = nanmean_float32(image, axis=0)
    limits = _hydrodynamic_limits(
        profile,
        settings,
        pixel_size_mm,
    )
    return limits if limits is not None else (0, max(profile.size - 1, 0))


def _hydrodynamic_limits(
    profile: np.ndarray,
    settings: CrossSectionSignalSettings,
    pixel_size_mm: float,
) -> tuple[int, int] | None:
    if profile.size == 0 or np.all(~np.isfinite(profile)):
        return None
    threshold = settings.velocity_profile_threshold * float(np.nanmax(profile))
    central = np.flatnonzero(profile > threshold)
    if central.size < 3:
        return None
    center = float(np.mean(central))
    roots = _half_height_roots(profile, central, center, threshold, pixel_size_mm)
    if roots is None:
        return None
    c1 = max(int(np.ceil(center + roots[0] / pixel_size_mm)), 0)
    c2 = min(int(np.floor(center + roots[1] / pixel_size_mm)), profile.size - 1)
    return (c1, c2) if c1 <= c2 else None


def _half_height_roots(
    profile: np.ndarray,
    central: np.ndarray,
    center: float,
    threshold: float,
    pixel_size_mm: float,
) -> tuple[float, float] | None:
    x = (central.astype(np.float32) - np.float32(center)) * np.float32(pixel_size_mm)
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rotated = _rotate_stack_with_nan(sub_stack, angle)
    transverse_profiles = nanmean_float32(rotated, axis=1)
    longitudinal_profiles = nanmean_float32(rotated, axis=2)
    raw = nanmean_float32(transverse_profiles[:, c1 : c2 + 1], axis=1)
    raw = np.where(np.isnan(raw), np.float32(0.0), raw).astype(
        np.float32,
        copy=False,
    )
    safe_velocity = nanmean_float32(transverse_profiles, axis=1)
    return (
        raw,
        safe_velocity,
        transverse_profiles,
        longitudinal_profiles,
        rotated,
    )


def _rotate_stack_with_nan(sub_stack: np.ndarray, angle: float) -> np.ndarray:
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            gpu_stack = backend.cupy.asarray(sub_stack, dtype=backend.cupy.float32)
            valid = backend.cupy.isfinite(gpu_stack)
            rotated_values = backend.ndimage.rotate(
                backend.cupy.where(valid, gpu_stack, backend.cupy.float32(0.0)),
                angle,
                axes=(-2, -1),
                reshape=False,
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=False,
            )
            rotated_weights = backend.ndimage.rotate(
                valid.astype(backend.cupy.float32),
                angle,
                axes=(-2, -1),
                reshape=False,
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=False,
            )
            rotated = backend.cupy.full(
                rotated_values.shape,
                backend.cupy.nan,
                dtype=backend.cupy.float32,
            )
            backend.cupy.divide(
                rotated_values,
                rotated_weights,
                out=rotated,
                where=rotated_weights >= backend.cupy.float32(0.5),
            )
            return backend.cupy.asnumpy(rotated)
        except Exception:
            pass

    shared_validity = _shared_validity_mask(sub_stack)
    if shared_validity is not None:
        return _rotate_stack_with_shared_validity_cpu(
            sub_stack,
            shared_validity,
            angle,
        )

    valid = np.isfinite(sub_stack)
    rotated_values = ndi.rotate(
        np.where(valid, sub_stack, np.float32(0.0)),
        angle,
        axes=(-2, -1),
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    )
    rotated_weights = ndi.rotate(
        valid.astype(np.float32),
        angle,
        axes=(-2, -1),
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    )
    rotated = np.full(rotated_values.shape, np.nan, dtype=np.float32)
    np.divide(
        rotated_values,
        rotated_weights,
        out=rotated,
        where=rotated_weights >= np.float32(0.5),
    )
    return rotated.astype(np.float32, copy=False)


def _rotate_stack_with_shared_validity_cpu(
    sub_stack: np.ndarray,
    shared_validity: np.ndarray,
    angle: float,
) -> np.ndarray:
    """Rotate a temporal stack while rotating its shared validity only once."""
    rotated_weights = ndi.rotate(
        shared_validity.astype(np.float32),
        angle,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    ).astype(np.float32, copy=False)
    finite_output = rotated_weights >= np.float32(0.5)
    finite_y, finite_x = np.nonzero(finite_output)
    rotated = np.full(sub_stack.shape, np.nan, dtype=np.float32)
    if finite_y.size == 0:
        return rotated

    y_start = int(finite_y.min())
    y_stop = int(finite_y.max()) + 1
    x_start = int(finite_x.min())
    x_stop = int(finite_x.max()) + 1
    output_shape = (y_stop - y_start, x_stop - x_start)
    rotation_matrix, rotation_offset = _rotation_affine(
        shared_validity.shape,
        angle,
    )
    cropped_offset = rotation_offset + rotation_matrix @ np.asarray(
        [y_start, x_start],
        dtype=np.float64,
    )
    rotated_values = np.empty(
        (sub_stack.shape[0], *output_shape),
        dtype=np.float32,
    )
    for frame_index, frame in enumerate(sub_stack):
        ndi.affine_transform(
            np.where(shared_validity, frame, np.float32(0.0)),
            rotation_matrix,
            cropped_offset,
            output_shape=output_shape,
            output=rotated_values[frame_index],
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        )
    cropped_weights = rotated_weights[y_start:y_stop, x_start:x_stop]
    np.divide(
        rotated_values,
        cropped_weights[None, :, :],
        out=rotated[:, y_start:y_stop, x_start:x_stop],
        where=cropped_weights[None, :, :] >= np.float32(0.5),
    )
    return rotated


def _rotation_affine(
    image_shape: tuple[int, int],
    angle: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return SciPy ``rotate(..., reshape=False)`` matrix and offset."""
    cosine = special.cosdg(angle)
    sine = special.sindg(angle)
    matrix = np.asarray(
        [[cosine, sine], [-sine, cosine]],
        dtype=np.float64,
    )
    center = (np.asarray(image_shape, dtype=np.float64) - 1.0) / 2.0
    return matrix, center - matrix @ center


def _shared_validity_mask(values: np.ndarray) -> np.ndarray | None:
    """Return a 2-D validity mask when every frame has the same finite pixels."""
    source = np.asarray(values)
    if source.ndim != 3 or source.shape[0] == 0:
        return None
    first = np.isfinite(source[0])
    for frame in source[1:]:
        if not np.array_equal(np.isfinite(frame), first):
            return None
    return first
