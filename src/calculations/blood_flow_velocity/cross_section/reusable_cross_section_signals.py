"""Reusable cross-section projections for spatially registered data cubes."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from calculations.topology import (
    AnnulusGeometry, BranchIdentityResult, OpticDisc, PreparedTopology,
    prepare_segment_chunks, prepare_topologies, resolve_segment_rotations,
)
from .generate_cross_section_signals import (
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    _cross_section_worker_count,
    _generate_cross_section_signals_from_prepared,
)


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

    prepared_topology: PreparedTopology
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
    optic_disc: OpticDisc
    ring_settings: AnnulusGeometry
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
    optic_disc: OpticDisc,
    ring_settings: AnnulusGeometry,
    cross_section_settings: CrossSectionSignalSettings,
) -> tuple[CrossSectionProjectionPlan, CrossSectionSignalResult]:
    """Prepare shared topology once and measure the reference cube."""
    vessel = np.asarray(vessel_mask, dtype=bool)
    if vessel.ndim != 2:
        raise ValueError("vessel_mask must be a 2-D array.")
    _validate_data_cube(reference_cube, vessel.shape, "reference_cube")
    prepared = prepare_topologies(
        {"vessel": vessel}, optic_disc, ring_settings, source_id="",
        window_size_percentile_kept=cross_section_settings.submask_size_percentile_kept,
    )["vessel"]
    prepared = resolve_segment_rotations(
        prepared,
        reference_cube,
        working_memory_mb=cross_section_settings.working_memory_mb,
    )
    result = _project_prepared(reference_cube, prepared, ring_settings, cross_section_settings)
    plan = CrossSectionProjectionPlan(
        prepared_topology=prepared,
        spatial_shape=vessel.shape,
        labels=result.labels.copy(),
        branch_ids=result.branch_ids.copy(),
        segment_center_xy=result.segment_center_xy.copy(),
        profile_window_bounds_xyxy=result.profile_window_bounds_xyxy.copy(),
        profile_window_side_pixels=result.profile_window_side_pixels,
        profile_pixel_size_mm=result.profile_pixel_size_mm,
        profile_rotation_degrees=result.profile_rotation_degrees.copy(),
        profile_integration_limits_pixels=result.profile_integration_limits_pixels.copy(),
        valid_segments=result.topology.valid_segments.copy(),
        branch_identity=result.branch_identity,
        optic_disc=optic_disc,
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
    """Project a registered cube through the same topology and transforms.

    Both limits modes use the full transverse width, matching the refactor's
    removal of hydrodynamic fitting from segment extraction.
    """
    _validate_limits_mode(limits_mode)
    _validate_data_cube(data_cube, plan.spatial_shape, "data_cube")
    return _project_prepared(
        data_cube, plan.prepared_topology, plan.ring_settings, plan.cross_section_settings,
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
    optic_disc: OpticDisc,
    ring_settings: AnnulusGeometry,
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
        optic_disc,
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


def _project_prepared(data_cube, prepared, ring_settings, settings):
    workers = _cross_section_worker_count(
        int(np.count_nonzero(prepared.topology.valid_segments)),
        frame_count=int(data_cube.shape[0]),
        working_memory_mb=settings.working_memory_mb,
    )
    return _generate_cross_section_signals_from_prepared(
        data_cube, prepared,
        prepare_segment_chunks(
            data_cube,
            prepared,
            worker_count=workers,
            working_memory_mb=settings.working_memory_mb,
            include_masked_before_rotation=True,
        ),
        ring_settings, settings,
    )


def _validate_data_cube(
    data_cube,
    spatial_shape: tuple[int, int],
    name: str,
) -> None:
    shape = getattr(data_cube, "shape", None)
    if shape is None or len(shape) != 3:
        raise ValueError(f"{name} must have shape (frame, y, x), got {shape!r}.")
    if any(int(size) <= 0 for size in shape):
        raise ValueError(f"{name} axes must be nonempty.")
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
