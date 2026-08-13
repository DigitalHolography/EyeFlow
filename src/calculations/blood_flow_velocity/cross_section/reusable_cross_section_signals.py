"""Reusable cross-section projections for spatially registered data cubes."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from calculations.math import nanmean_float32

from .branch_identity import BranchIdentityResult
from .generate_cross_section_signals import (
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    _cross_section_limits,
    _CrossSectionBuffers,
    _CrossSectionMeasurement,
    _center_pad_for_rotation,
    _empty_result,
    _profile_measurement,
    _resize_subimage_stack,
    _resize_submask,
    _result_from_buffers,
    _rotated_profile_sample_count,
    _rotate_mean_image,
    _rotate_mask,
    _rotate_masked_image,
    _store_cross_section_measurement,
    _subimage_stack_from_bounds,
    generate_cross_section_signals,
)
from .segment_geometry import SegmentRingSettings, section_masks

__all__ = [
    "CrossSectionProjectionPlan",
    "MultiCubeCrossSectionSignalResult",
    "fit_cross_section_plan",
    "generate_cross_section_signals_for_cubes",
    "project_cross_section_cube",
    "project_cross_section_cubes",
]


@dataclass(frozen=True)
class CrossSectionProjectionPlan:
    """Projection geometry fitted once on a reference data cube.

    Segment-indexed arrays use ``(ring, branch, ...)`` axes, except
    ``segment_center_xy``, which retains EyeFlow's existing
    ``(branch, ring, xy)`` layout. Window bounds are stored as
    clipped ``(x_start, x_stop, y_start, y_stop)`` source coordinates with
    exclusive stops. ``profile_window_side_pixels`` plus each segment center
    reconstruct the full centroid-centered square, including padded pixels.
    """

    spatial_shape: tuple[int, int]
    labels: np.ndarray
    branch_ids: np.ndarray
    segment_center_xy: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_window_side_pixels: int
    profile_pixel_size_mm: float
    profile_rotation_degrees: np.ndarray
    profile_integration_limits_pixels: np.ndarray
    valid_segments: np.ndarray
    branch_identity: BranchIdentityResult
    optic_disc_center: object
    ring_settings: SegmentRingSettings
    cross_section_settings: CrossSectionSignalSettings

    @property
    def rotation_angles_degrees(self) -> np.ndarray:
        """Return every fitted square-window rotation angle."""
        return self.profile_rotation_degrees


@dataclass(frozen=True)
class MultiCubeCrossSectionSignalResult:
    """One fitted projection plan and the result of every requested pass."""

    plan: CrossSectionProjectionPlan
    reference_name: str
    passes: dict[str, CrossSectionSignalResult]

    @property
    def reference(self) -> CrossSectionSignalResult:
        return self.passes[self.reference_name]


def fit_cross_section_plan(
    reference_cube,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
) -> tuple[CrossSectionProjectionPlan, CrossSectionSignalResult]:
    """Fit one fixed window size, angles, and limits on a reference cube."""
    vessel = np.asarray(vessel_mask, dtype=bool)
    if vessel.ndim != 2:
        raise ValueError("vessel_mask must be a 2-D array.")
    _validate_data_cube(reference_cube, vessel.shape, "reference_cube")

    result = generate_cross_section_signals(
        reference_cube,
        vessel,
        optic_disc_center,
        ring_settings,
        cross_section_settings,
    )
    if result.branch_ids.size == 0:
        return (
            _empty_projection_plan(
                vessel,
                result.branch_identity,
                optic_disc_center,
                ring_settings,
                cross_section_settings,
            ),
            result,
        )

    center_valid = np.all(np.isfinite(result.segment_center_xy), axis=-1).T
    bounds = result.profile_window_bounds_xyxy
    limits = result.profile_integration_limits_pixels
    valid_segments = (
        center_valid
        & np.isfinite(result.profile_rotation_degrees)
        & np.all(bounds >= 0, axis=-1)
        & (bounds[..., 0] < bounds[..., 1])
        & (bounds[..., 2] < bounds[..., 3])
        & (limits[..., 0] >= 0)
        & (limits[..., 0] <= limits[..., 1])
    )
    plan = CrossSectionProjectionPlan(
        spatial_shape=vessel.shape,
        labels=result.labels.copy(),
        branch_ids=result.branch_ids.copy(),
        segment_center_xy=result.segment_center_xy.copy(),
        profile_window_bounds_xyxy=bounds.copy(),
        profile_window_side_pixels=result.profile_window_side_pixels,
        profile_pixel_size_mm=result.profile_pixel_size_mm,
        profile_rotation_degrees=result.profile_rotation_degrees.copy(),
        profile_integration_limits_pixels=limits.copy(),
        valid_segments=valid_segments,
        branch_identity=result.branch_identity,
        optic_disc_center=_copy_optional_array(optic_disc_center),
        ring_settings=ring_settings,
        cross_section_settings=cross_section_settings,
    )
    return plan, result


def project_cross_section_cube(
    data_cube,
    plan: CrossSectionProjectionPlan,
    *,
    limits_mode: str = "reference",
) -> CrossSectionSignalResult:
    """Apply a fitted cross-section plan to one registered data cube.

    ``limits_mode='reference'`` reuses the reference cube's transverse
    integration limits. ``limits_mode='per_cube'`` keeps the reference window
    and angle but estimates new limits from the projected cube.
    """
    _validate_limits_mode(limits_mode)
    _validate_data_cube(data_cube, plan.spatial_shape, "data_cube")
    if plan.branch_ids.size == 0:
        vessel = np.zeros(plan.spatial_shape, dtype=bool)
        return _empty_result(
            data_cube,
            vessel,
            plan.ring_settings,
            plan.branch_identity,
            substack_side_pixels=plan.profile_window_side_pixels,
            profile_pixel_size_mm=plan.profile_pixel_size_mm,
        )

    masks = section_masks(
        plan.spatial_shape,
        plan.optic_disc_center,
        plan.ring_settings,
    )
    buffers = _CrossSectionBuffers.allocate(
        frame_count=data_cube.shape[0],
        ring_count=plan.ring_settings.ring_count,
        branch_count=plan.branch_ids.size,
    )
    buffers.segment_center_xy[...] = plan.segment_center_xy

    for circle_index, section in enumerate(masks):
        for branch_index, branch_id in enumerate(plan.branch_ids):
            if not plan.valid_segments[circle_index, branch_index]:
                continue
            mask = section & (plan.labels == int(branch_id))
            loc_values = plan.segment_center_xy[branch_index, circle_index]
            loc_xy = (int(loc_values[0]), int(loc_values[1]))
            bounds_xyxy = tuple(
                int(value)
                for value in plan.profile_window_bounds_xyxy[
                    circle_index,
                    branch_index,
                ]
            )
            sub_stack, sub_mask = _subimage_stack_from_bounds(
                data_cube,
                mask,
                bounds_xyxy,
                loc_xy=loc_xy,
                side_pixels=plan.profile_window_side_pixels,
            )
            resized_stack = _resize_subimage_stack(sub_stack)
            resized_mask = _resize_submask(sub_mask)
            resized_stack_masked = resized_stack.copy()
            resized_stack_masked[:, ~resized_mask] = np.nan
            angle = float(plan.profile_rotation_degrees[circle_index, branch_index])
            mean_image = nanmean_float32(resized_stack, axis=0)
            mean_image_masked = nanmean_float32(resized_stack_masked, axis=0)
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
            if limits_mode == "per_cube":
                limits = _cross_section_limits(
                    rotated_mean_masked,
                    plan.cross_section_settings,
                    pixel_size_mm=plan.profile_pixel_size_mm,
                )
            else:
                limits = tuple(
                    int(value)
                    for value in plan.profile_integration_limits_pixels[
                        circle_index,
                        branch_index,
                    ]
                )
            c1, c2 = limits
            measurement = _CrossSectionMeasurement(
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
                rotated_mask=_rotate_mask(rotation_mask, angle),
                limits=limits,
                sample_count=_rotated_profile_sample_count(angle),
            )
            _store_cross_section_measurement(
                buffers,
                circle_index,
                branch_index,
                measurement,
                bounds_xyxy,
            )

    return _result_from_buffers(
        buffers,
        plan.branch_identity,
        plan.cross_section_settings,
        plan.profile_window_side_pixels,
    )


def project_cross_section_cubes(
    cubes: Mapping[str, object],
    plan: CrossSectionProjectionPlan,
    *,
    limits_mode: str = "reference",
) -> dict[str, CrossSectionSignalResult]:
    """Apply one fitted plan sequentially to any number of named cubes."""
    _validate_limits_mode(limits_mode)
    results: dict[str, CrossSectionSignalResult] = {}
    for name, cube in cubes.items():
        _validate_pass_name(name)
        results[name] = project_cross_section_cube(
            cube,
            plan,
            limits_mode=limits_mode,
        )
    return results


def generate_cross_section_signals_for_cubes(
    reference_cube,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    additional_cubes: Mapping[str, object] | None = None,
    reference_name: str = "reference",
    limits_mode: str = "reference",
) -> MultiCubeCrossSectionSignalResult:
    """Fit one reference and project zero or more additional named cubes."""
    _validate_pass_name(reference_name, parameter="reference_name")
    _validate_limits_mode(limits_mode)
    extras = {} if additional_cubes is None else dict(additional_cubes)
    if reference_name in extras:
        raise ValueError(f"additional_cubes contains the reference name {reference_name!r}.")

    plan, reference_result = fit_cross_section_plan(
        reference_cube,
        vessel_mask,
        optic_disc_center,
        ring_settings,
        cross_section_settings,
    )
    passes = {reference_name: reference_result}
    passes.update(
        project_cross_section_cubes(
            extras,
            plan,
            limits_mode=limits_mode,
        )
    )
    return MultiCubeCrossSectionSignalResult(plan, reference_name, passes)


def _empty_projection_plan(
    vessel: np.ndarray,
    branches: BranchIdentityResult,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
) -> CrossSectionProjectionPlan:
    segment_shape = (ring_settings.ring_count, 0)
    return CrossSectionProjectionPlan(
        spatial_shape=vessel.shape,
        labels=branches.labels.copy(),
        branch_ids=branches.branch_ids.copy(),
        segment_center_xy=np.full(
            (0, ring_settings.ring_count, 2),
            np.nan,
            dtype=np.float32,
        ),
        profile_window_bounds_xyxy=np.full(
            (*segment_shape, 4),
            -1,
            dtype=np.int32,
        ),
        profile_window_side_pixels=0,
        profile_pixel_size_mm=0.0,
        profile_rotation_degrees=np.full(
            segment_shape,
            np.nan,
            dtype=np.float32,
        ),
        profile_integration_limits_pixels=np.full(
            (*segment_shape, 2),
            -1,
            dtype=np.int32,
        ),
        valid_segments=np.zeros(segment_shape, dtype=bool),
        branch_identity=branches,
        optic_disc_center=_copy_optional_array(optic_disc_center),
        ring_settings=ring_settings,
        cross_section_settings=cross_section_settings,
    )


def _validate_data_cube(
    data_cube,
    spatial_shape: tuple[int, int],
    name: str,
) -> None:
    shape = getattr(data_cube, "shape", None)
    if shape is None or len(shape) != 3:
        raise ValueError(f"{name} must have shape (frame, y, x), got {shape!r}.")
    if tuple(shape[-2:]) != tuple(spatial_shape):
        raise ValueError(
            f"{name} spatial shape must be {tuple(spatial_shape)}, got {tuple(shape[-2:])}."
        )


def _validate_limits_mode(limits_mode: str) -> None:
    if limits_mode not in {"reference", "per_cube"}:
        raise ValueError("limits_mode must be either 'reference' or 'per_cube'.")


def _validate_pass_name(name, *, parameter: str = "cube name") -> None:
    if not isinstance(name, str) or not name:
        raise ValueError(f"{parameter} must be a non-empty string.")


def _copy_optional_array(value):
    return None if value is None else np.asarray(value).copy()
