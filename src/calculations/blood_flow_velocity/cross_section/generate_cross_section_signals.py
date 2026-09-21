"""Profile measurements on segments prepared by calculations.topology."""

from __future__ import annotations

from collections import deque
from collections.abc import Mapping
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, replace

import numpy as np
from scipy import special

from calculations.compute_backend import optional_cupy_backend
from calculations.math import (
    nanmean_float32,
)
from calculations.topology import (
    BranchIdentityResult,
    PreparedSegmentChunks,
    PreparedTopology,
    SegmentRingSettings,
    extract_segment,
    prepare_topologies,
    prepare_segment_chunks,
    resolve_segment_rotations,
    resample_rotate_segment,
    dilate_segment_masks,
    longitudinal_profiles as calculate_longitudinal_profiles,
    transverse_profiles as calculate_transverse_profiles,
)

from runtime_limits import cap_parallel_jobs
@dataclass(frozen=True)
class CrossSectionSignalSettings:
    pixel_size_mm: float
    submask_size_percentile_kept: float = 0.95
    working_memory_mb: float = 512.0

    def __post_init__(self):
        if not np.isfinite(self.pixel_size_mm) or self.pixel_size_mm <= 0:
            raise ValueError("pixel_size_mm must be finite and positive.")
        if not np.isfinite(self.working_memory_mb) or self.working_memory_mb <= 0:
            raise ValueError("working_memory_mb must be finite and positive.")
        if not 0 < self.submask_size_percentile_kept <= 1:
            raise ValueError("submask_size_percentile_kept must be in (0, 1].")


@dataclass(frozen=True, kw_only=True)
class CrossSectionProfileOutputs:
    """Transverse and longitudinal cross-section profile outputs."""

    velocity_profiles: np.ndarray
    transverse_velocity_profiles_masked: np.ndarray
    longitudinal_velocity_profiles_unmasked: np.ndarray
    longitudinal_velocity_profiles_masked: np.ndarray
    profile_sample_count: np.ndarray
    profile_rotation_degrees: np.ndarray
    rotated_mean_images: np.ndarray
    rotated_mean_images_masked: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_window_side_pixels: int
    profile_pixel_size_mm: float
    profile_integration_limits_pixels: np.ndarray


@dataclass(frozen=True)
class CrossSectionTopology:
    spatial_shape: tuple[int, int]
    frame_count: int
    labels: np.ndarray
    branch_ids: np.ndarray
    section_masks: np.ndarray
    segment_masks: np.ndarray
    segment_center_xy: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_window_side_pixels: int
    profile_pixel_size_mm: float
    profile_rotation_degrees: np.ndarray
    profile_integration_limits_pixels: np.ndarray
    valid_segments: np.ndarray
    ring_settings: SegmentRingSettings
    branch_identity: BranchIdentityResult
    prepared_topology: PreparedTopology | None = None


@dataclass(frozen=True, kw_only=True)
class CrossSectionDisplacementResult:
    """Local waveforms, maps, sums, and cross-sectional motion metrics."""

    displacement: np.ndarray
    safe_displacement: np.ndarray
    displacement_maps_per_segment: np.ndarray | None
    transverse_displacement_profiles_unmasked: np.ndarray
    transverse_displacement_profiles_masked: np.ndarray
    longitudinal_displacement_profiles_unmasked: np.ndarray
    longitudinal_displacement_profiles_masked: np.ndarray
    x_sum_displacement_profile: np.ndarray
    y_sum_displacement_profile: np.ndarray
    cross_sectional_radial_movement_amplitude: np.ndarray
    cross_sectional_radial_asymmetry_index: np.ndarray


@dataclass(frozen=True)
class CrossSectionSignalResult(CrossSectionProfileOutputs):
    """Segment waveforms and their transverse profile measurements."""

    velocity: np.ndarray
    safe_velocity: np.ndarray
    velocity_maps_per_segment: np.ndarray | None
    velocity_map_segment_indexes: np.ndarray
    segment_masks: np.ndarray
    labels: np.ndarray
    branch_ids: np.ndarray
    segment_center_xy: np.ndarray
    branch_identity: BranchIdentityResult
    topology: CrossSectionTopology
    displacements: dict[str, CrossSectionDisplacementResult]
    transverse_velocity_fft_profiles_unmasked: np.ndarray | None = None
    transverse_velocity_fft_profiles_masked: np.ndarray | None = None

    @property
    def displacement_by_method(self) -> dict[str, CrossSectionDisplacementResult]:
        return self.displacements

    @property
    def projected_signal(self) -> np.ndarray:
        return self.velocity

    @property
    def full_profile_signal(self) -> np.ndarray:
        return self.safe_velocity

    @property
    def transverse_profiles(self) -> np.ndarray:
        return self.velocity_profiles


@dataclass(frozen=True)
class _CrossSectionVelocityMeasurement:
    raw: np.ndarray
    safe_velocity: np.ndarray
    transverse_profiles: np.ndarray
    longitudinal_profiles: np.ndarray
    rotated_stack: np.ndarray | None
    angle: float


@dataclass(frozen=True)
class _CrossSectionMeasurement:
    unmasked: _CrossSectionVelocityMeasurement
    masked: _CrossSectionVelocityMeasurement
    rotated_mean: np.ndarray
    rotated_mean_masked: np.ndarray
    rotated_mask: np.ndarray
    limits: tuple[int, int]
    sample_count: int


@dataclass(frozen=True)
class _DisplacementProfileMeasurement:
    raw: np.ndarray
    safe_displacement: np.ndarray


@dataclass(frozen=True)
class _CrossSectionDisplacementMeasurement:
    waveform: _DisplacementProfileMeasurement
    vectors: np.ndarray
    transverse_profiles_unmasked: np.ndarray
    transverse_profiles_masked: np.ndarray
    longitudinal_profiles_unmasked: np.ndarray
    longitudinal_profiles_masked: np.ndarray
    summed_displacement_xy: np.ndarray
    radial_movement_amplitude: np.ndarray
    radial_asymmetry_index: np.ndarray


@dataclass(frozen=True)
class _CrossSectionDisplacementWork:
    component_stacks: np.ndarray
    angle: float
    rotated_mask: np.ndarray
    limits: tuple[int, int]


_INTERPOLATED_SUBSTACK_SIDE = 128
_ROTATED_SUBSTACK_SIDE = int(_INTERPOLATED_SUBSTACK_SIDE * np.sqrt(2.0))
_MAX_PARALLEL_CROSS_SECTIONS = 8
_ARTERY_TRANSVERSE_MASK_DILATION_PIXELS = 10


@dataclass
class _CrossSectionBuffers:
    velocity: np.ndarray
    safe_velocity: np.ndarray
    velocity_maps_per_segment: np.ndarray | None
    velocity_map_segment_indexes: np.ndarray
    velocity_map_rows: np.ndarray
    segment_masks: np.ndarray
    segment_center_xy: np.ndarray
    velocity_profiles: np.ndarray
    transverse_velocity_profiles_masked: np.ndarray
    longitudinal_velocity_profiles_unmasked: np.ndarray
    longitudinal_velocity_profiles_masked: np.ndarray
    profile_sample_count: np.ndarray
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
        velocity_map_segment_indexes: np.ndarray,
        retain_velocity_maps: bool,
    ) -> _CrossSectionBuffers:
        signal_shape = (ring_count, branch_count, frame_count)
        profile_shape = (ring_count, branch_count)
        segment_indexes = np.asarray(velocity_map_segment_indexes, dtype=np.int32)
        if segment_indexes.ndim != 2 or segment_indexes.shape[1] != 2:
            raise ValueError("velocity-map segment indexes must have shape (segment, 2).")
        if segment_indexes.size and (
            np.any(segment_indexes[:, 0] < 0)
            or np.any(segment_indexes[:, 0] >= ring_count)
            or np.any(segment_indexes[:, 1] < 0)
            or np.any(segment_indexes[:, 1] >= branch_count)
        ):
            raise ValueError("velocity-map segment indexes are out of range.")
        if (
            segment_indexes.size
            and np.unique(segment_indexes, axis=0).shape[0]
            != segment_indexes.shape[0]
        ):
            raise ValueError("velocity-map segment indexes must be unique.")
        map_rows = np.full(profile_shape, -1, dtype=np.int32)
        if retain_velocity_maps and segment_indexes.size:
            map_rows[segment_indexes[:, 0], segment_indexes[:, 1]] = np.arange(
                segment_indexes.shape[0],
                dtype=np.int32,
            )
        return cls(
            velocity=np.full(signal_shape, np.nan, dtype=np.float32),
            safe_velocity=np.full(signal_shape, np.nan, dtype=np.float32),
            velocity_maps_per_segment=(
                np.full(
                    (
                        segment_indexes.shape[0],
                        frame_count,
                        _ROTATED_SUBSTACK_SIDE,
                        _ROTATED_SUBSTACK_SIDE,
                    ),
                    np.nan,
                    dtype=np.float32,
                )
                if retain_velocity_maps
                else None
            ),
            velocity_map_segment_indexes=(
                segment_indexes
                if retain_velocity_maps
                else np.empty((0, 2), dtype=np.int32)
            ),
            velocity_map_rows=map_rows,
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


@dataclass
class _CrossSectionDisplacementBuffers:
    displacement: np.ndarray
    safe_displacement: np.ndarray
    displacement_maps_per_segment: np.ndarray | None
    transverse_displacement_profiles_unmasked: np.ndarray
    transverse_displacement_profiles_masked: np.ndarray
    longitudinal_displacement_profiles_unmasked: np.ndarray
    longitudinal_displacement_profiles_masked: np.ndarray
    x_sum_displacement_profile: np.ndarray
    y_sum_displacement_profile: np.ndarray
    cross_sectional_radial_movement_amplitude: np.ndarray
    cross_sectional_radial_asymmetry_index: np.ndarray

    @classmethod
    def allocate(
        cls,
        *,
        frame_count: int,
        ring_count: int,
        branch_count: int,
        retain_maps: bool = True,
    ) -> _CrossSectionDisplacementBuffers:
        signal_shape = (ring_count, branch_count, frame_count)
        map_shape = (
            *signal_shape,
            _ROTATED_SUBSTACK_SIDE,
            _ROTATED_SUBSTACK_SIDE,
            2,
        )

        def filled(shape):
            return np.full(shape, np.nan, dtype=np.float32)

        profile_shape = (*signal_shape, _ROTATED_SUBSTACK_SIDE)
        return cls(
            displacement=filled(signal_shape),
            safe_displacement=filled(signal_shape),
            displacement_maps_per_segment=(filled(map_shape) if retain_maps else None),
            transverse_displacement_profiles_unmasked=filled(profile_shape),
            transverse_displacement_profiles_masked=filled(profile_shape),
            longitudinal_displacement_profiles_unmasked=filled(profile_shape),
            longitudinal_displacement_profiles_masked=filled(profile_shape),
            x_sum_displacement_profile=filled(signal_shape),
            y_sum_displacement_profile=filled(signal_shape),
            cross_sectional_radial_movement_amplitude=filled(signal_shape),
            cross_sectional_radial_asymmetry_index=filled(signal_shape),
        )


def _validate_velocity_map(velocity_map, vessel_mask: np.ndarray) -> None:
    if vessel_mask.ndim != 2:
        raise ValueError(
            f'vessel_mask must have shape (y, x), got {vessel_mask.shape!r}.'
        )
    shape = getattr(velocity_map, 'shape', None)
    if shape is None or len(shape) != 3:
        raise ValueError(
            f'velocity_map must have shape (frame, y, x), got {shape!r}.'
        )
    if tuple(shape[1:]) != tuple(vessel_mask.shape):
        raise ValueError(
            'velocity_map spatial shape must match vessel_mask: '
            f'{tuple(shape[1:])!r} != {tuple(vessel_mask.shape)!r}.'
        )


def _validate_displacement_maps(
    displacement_maps: Mapping[str, object] | None,
    velocity_map,
) -> dict[str, object]:
    if displacement_maps is None:
        return {}
    if not isinstance(displacement_maps, Mapping):
        raise ValueError('displacement_maps must be a mapping keyed by method.')
    expected_shape = (*tuple(velocity_map.shape), 2)
    normalized: dict[str, object] = {}
    for method, displacement_map in displacement_maps.items():
        if not isinstance(method, str) or not method:
            raise ValueError('displacement method names must be non-empty strings.')
        shape = getattr(displacement_map, 'shape', None)
        if shape is None or tuple(shape) != expected_shape:
            raise ValueError(
                f'displacement map {method!r} must have shape '
                f'{expected_shape!r}, got {shape!r}.'
            )
        normalized[method] = displacement_map
    return normalized


def generate_cross_section_signals(
    velocity_map,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    optic_disc_mask=None,
    displacement_maps: Mapping[str, object] | None = None,
    retain_displacement_maps: bool = True,
    retain_velocity_maps: bool = True,
) -> CrossSectionSignalResult:
    """Measure a scalar cube using the shared topology and streamed transforms."""

    vessel = np.asarray(vessel_mask, dtype=bool)
    _validate_velocity_map(velocity_map, vessel)
    prepared = prepare_topologies(
        {"vessel": vessel},
        optic_disc_mask,
        ring_settings,
        source_id="",
        optic_disc_center=optic_disc_center,
        window_size_percentile_kept=cross_section_settings.submask_size_percentile_kept,
    )["vessel"]
    prepared = resolve_segment_rotations(
        prepared,
        velocity_map,
        working_memory_mb=cross_section_settings.working_memory_mb,
    )
    workers = _cross_section_worker_count(
        int(np.count_nonzero(prepared.topology.valid_segments)),
        frame_count=int(velocity_map.shape[0]),
        working_memory_mb=cross_section_settings.working_memory_mb,
    )
    segments = prepare_segment_chunks(
        velocity_map,
        prepared,
        worker_count=workers,
        working_memory_mb=cross_section_settings.working_memory_mb,
        keep_on_device=optional_cupy_backend() is not None,
        include_masked_before_rotation=True,
    )
    return _generate_cross_section_signals_from_prepared(
        velocity_map, prepared, segments, ring_settings, cross_section_settings,
        displacement_maps=displacement_maps,
        retain_displacement_maps=retain_displacement_maps,
        retain_velocity_maps=retain_velocity_maps,
    )

def _generate_cross_section_signals_from_prepared(
    velocity_map,
    prepared_topology: PreparedTopology,
    prepared_segments: PreparedSegmentChunks,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    displacement_maps: Mapping[str, object] | None = None,
    retain_displacement_maps: bool = True,
    retain_velocity_maps: bool = True,
    segment_observer=None,
    transverse_mask_dilation_pixels: int = 0,
) -> CrossSectionSignalResult:
    segment_topology = prepared_topology.topology
    branches = segment_topology.branch_identity
    if branches is None:
        raise ValueError("prepared topology must retain its branch identity result.")

    vessel = np.asarray(branches.stages.vessel, dtype=bool)
    _validate_velocity_map(velocity_map, vessel)
    normalized_displacements = _validate_displacement_maps(
        displacement_maps,
        velocity_map,
    )
    substack_side_pixels = segment_topology.window_side_pixels
    profile_pixel_size_mm = _interpolated_pixel_size_mm(
        cross_section_settings.pixel_size_mm,
        substack_side_pixels,
    )
    if branches.branch_ids.size == 0:
        return _empty_result(
            velocity_map,
            vessel,
            ring_settings,
            branches,
            substack_side_pixels=substack_side_pixels,
            profile_pixel_size_mm=profile_pixel_size_mm,
            section_masks=segment_topology.annulus_masks,
            prepared_topology=prepared_topology,
            displacement_maps=normalized_displacements,
            retain_velocity_maps=retain_velocity_maps,
        )

    velocity_map_indexes = np.argwhere(
        segment_topology.valid_segments
        & np.isfinite(prepared_topology.rotation_degrees)
    )
    buffers = _CrossSectionBuffers.allocate(
        frame_count=velocity_map.shape[0],
        ring_count=ring_settings.ring_count,
        branch_count=branches.branch_ids.size,
        velocity_map_segment_indexes=velocity_map_indexes,
        retain_velocity_maps=retain_velocity_maps,
    )
    _fill_cross_section_buffers_from_prepared(
        buffers,
        prepared_topology,
        prepared_segments,
        cross_section_settings,
        substack_side_pixels,
        segment_observer=segment_observer,
        transverse_mask_dilation_pixels=transverse_mask_dilation_pixels,
    )
    topology = _legacy_topology_from_prepared(
        buffers,
        prepared_topology,
        ring_settings,
        frame_count=velocity_map.shape[0],
        profile_pixel_size_mm=profile_pixel_size_mm,
    )
    displacement_results = {
        method: _project_displacement_map(
            displacement_map,
            topology,
            retain_maps=retain_displacement_maps,
        )
        for method, displacement_map in normalized_displacements.items()
    }
    return _result_from_buffers(
        buffers,
        branches,
        cross_section_settings,
        substack_side_pixels,
        topology=topology,
        displacements=displacement_results,
    )


def _fill_cross_section_buffers_from_prepared(
    buffers: _CrossSectionBuffers,
    prepared_topology: PreparedTopology,
    prepared_segments: PreparedSegmentChunks,
    settings: CrossSectionSignalSettings,
    substack_side_pixels: int,
    *,
    segment_observer=None,
    transverse_mask_dilation_pixels: int = 0,
) -> None:
    """Reduce transformed chunks directly into retained output arrays."""

    topology = prepared_topology.topology
    mean_sum = np.zeros(buffers.rotated_mean_images.shape, dtype=np.float64)
    mean_count = np.zeros(buffers.rotated_mean_images.shape, dtype=np.int64)
    masked_sum = np.zeros(buffers.rotated_mean_images.shape, dtype=np.float64)
    masked_count = np.zeros(buffers.rotated_mean_images.shape, dtype=np.int64)
    initialized: set[tuple[int, int]] = set()
    for prepared_segment in prepared_segments:
        index = (prepared_segment.ring_index, prepared_segment.branch_index)
        frame_slice = prepared_segment.frame_slice
        angle = float(prepared_topology.rotation_degrees[index])
        rotated = prepared_segment.rotated
        rotated_masked = prepared_segment.rotated_masked
        if rotated_masked is None:
            raise ValueError(
                "velocity segments must include the pre-rotation masked companion."
            )
        rotated_mask = prepared_topology.rotated_masks[index]
        profile_mask = _dilate_profile_mask(
            rotated_mask,
            transverse_mask_dilation_pixels,
        )
        if segment_observer is not None:
            segment_observer(
                index[0],
                index[1],
                frame_slice,
                rotated,
                profile_mask,
            )
        _accumulate_rotated_means(
            rotated,
            rotated_masked,
            rotated_mask,
            mean_sum[index],
            mean_count[index],
            masked_sum[index],
            masked_count[index],
        )
        backend = optional_cupy_backend()
        if backend is not None and isinstance(rotated, backend.cupy.ndarray):
            measurement = _gpu_cross_section_measurement(
                rotated,
                rotated_masked,
                rotated_mask,
                profile_mask,
                angle,
                use_dilated_transverse_profile=(
                    transverse_mask_dilation_pixels > 0
                ),
                retain_rotated_stack=buffers.velocity_maps_per_segment is not None,
                cupy=backend.cupy,
            )
        else:
            rotated = np.asarray(rotated, dtype=np.float32)
            rotated_masked = np.asarray(rotated_masked, dtype=np.float32)
            transverse = calculate_transverse_profiles(rotated)
            longitudinal = calculate_longitudinal_profiles(rotated)
            transverse_masked = calculate_transverse_profiles(rotated_masked)
            longitudinal_masked = calculate_longitudinal_profiles(rotated_masked)
            profile_transverse = (
                _nanmean_float32_where(
                    rotated,
                    np.isfinite(rotated) & profile_mask[None, ...],
                    axis=1,
                )
                if transverse_mask_dilation_pixels > 0
                else transverse_masked
            )
            masked_measurement = _profile_measurement_from_profiles(
                transverse_masked,
                longitudinal_masked,
                angle,
                0,
                rotated.shape[-1] - 1,
                rotated_stack=None,
            )
            if transverse_mask_dilation_pixels > 0:
                masked_measurement = replace(
                    masked_measurement,
                    transverse_profiles=profile_transverse,
                )
            finite_masked_mean = (
                np.isfinite(rotated_masked) & rotated_mask[None, ...]
            )
            measurement = _CrossSectionMeasurement(
                unmasked=_profile_measurement_from_profiles(
                    transverse,
                    longitudinal,
                    angle,
                    0,
                    rotated.shape[-1] - 1,
                    rotated_stack=rotated,
                ),
                masked=masked_measurement,
                rotated_mean=nanmean_float32(rotated, axis=0),
                rotated_mean_masked=_nanmean_float32_where(
                    rotated_masked,
                    finite_masked_mean,
                    axis=0,
                ),
                rotated_mask=rotated_mask,
                limits=(0, rotated.shape[-1] - 1),
                sample_count=_rotated_profile_sample_count(angle),
            )
        buffers.segment_center_xy[index[1], index[0]] = (
            topology.segment_centers_xy[index]
        )
        _store_cross_section_measurement(
            buffers,
            index[0],
            index[1],
            measurement,
            tuple(int(value) for value in topology.window_bounds_xyxy[index]),
            frame_slice=frame_slice,
            initialize=index not in initialized,
        )
        initialized.add(index)

    np.divide(
        mean_sum,
        mean_count,
        out=buffers.rotated_mean_images,
        where=mean_count > 0,
    )
    np.divide(
        masked_sum,
        masked_count,
        out=buffers.rotated_mean_images_masked,
        where=masked_count > 0,
    )


def _accumulate_rotated_means(
    rotated,
    rotated_masked,
    rotated_mask: np.ndarray,
    total: np.ndarray,
    count: np.ndarray,
    masked_total: np.ndarray,
    masked_count: np.ndarray,
) -> None:
    backend = optional_cupy_backend()
    if backend is not None and isinstance(rotated, backend.cupy.ndarray):
        cupy = backend.cupy
        finite = cupy.isfinite(rotated)
        values = cupy.where(finite, rotated, cupy.float32(0.0))
        chunk_total = cupy.asnumpy(cupy.sum(values, axis=0, dtype=cupy.float64))
        chunk_count = cupy.asnumpy(cupy.sum(finite, axis=0, dtype=cupy.int64))
        masked_finite = cupy.isfinite(rotated_masked) & cupy.asarray(
            rotated_mask, dtype=cupy.bool_
        )[None]
        masked_values = cupy.where(
            masked_finite, rotated_masked, cupy.float32(0.0)
        )
        chunk_masked_total = cupy.asnumpy(
            cupy.sum(masked_values, axis=0, dtype=cupy.float64)
        )
        chunk_masked_count = cupy.asnumpy(
            cupy.sum(masked_finite, axis=0, dtype=cupy.int64)
        )
    else:
        values = np.asarray(rotated, dtype=np.float32)
        finite = np.isfinite(values)
        chunk_total = np.sum(
            values,
            axis=0,
            dtype=np.float64,
            where=finite,
        )
        chunk_count = np.sum(finite, axis=0, dtype=np.int64)
        masked_values = np.asarray(rotated_masked, dtype=np.float32)
        masked_finite = np.isfinite(masked_values) & rotated_mask[None]
        chunk_masked_total = np.sum(
            masked_values,
            axis=0,
            dtype=np.float64,
            where=masked_finite,
        )
        chunk_masked_count = np.sum(masked_finite, axis=0, dtype=np.int64)
    total += chunk_total
    count += chunk_count
    masked_total += chunk_masked_total
    masked_count += chunk_masked_count

def _result_from_buffers(
    buffers: _CrossSectionBuffers,
    branches: BranchIdentityResult,
    settings: CrossSectionSignalSettings,
    substack_side_pixels: int,
    *,
    topology: CrossSectionTopology,
    displacements: dict[str, CrossSectionDisplacementResult],
) -> CrossSectionSignalResult:
    profile_pixel_size_mm = _interpolated_pixel_size_mm(
        settings.pixel_size_mm,
        substack_side_pixels,
    )
    return CrossSectionSignalResult(
        velocity=buffers.velocity,
        safe_velocity=buffers.safe_velocity,
        velocity_maps_per_segment=buffers.velocity_maps_per_segment,
        velocity_map_segment_indexes=buffers.velocity_map_segment_indexes,
        segment_masks=buffers.segment_masks,
        labels=branches.labels,
        branch_ids=branches.branch_ids,
        segment_center_xy=buffers.segment_center_xy,
        branch_identity=branches,
        topology=topology,
        displacements=displacements,
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
        profile_sample_count=buffers.profile_sample_count,
        profile_rotation_degrees=buffers.profile_rotation_degrees,
        rotated_mean_images=buffers.rotated_mean_images,
        rotated_mean_images_masked=buffers.rotated_mean_images_masked,
        profile_window_bounds_xyxy=buffers.profile_window_bounds_xyxy,
        profile_window_side_pixels=int(substack_side_pixels),
        profile_pixel_size_mm=profile_pixel_size_mm,
        profile_integration_limits_pixels=(buffers.profile_integration_limits_pixels),
    )


def _empty_result(
    velocity_map: np.ndarray,
    vessel: np.ndarray,
    settings: SegmentRingSettings,
    branches: BranchIdentityResult,
    *,
    substack_side_pixels: int,
    profile_pixel_size_mm: float,
    section_masks: np.ndarray,
    prepared_topology: PreparedTopology | None = None,
    displacement_maps: Mapping[str, object],
    retain_velocity_maps: bool,
) -> CrossSectionSignalResult:
    shape = (settings.ring_count, 0, velocity_map.shape[0])
    empty_profiles = np.full(
        (
            settings.ring_count,
            0,
            velocity_map.shape[0],
            _ROTATED_SUBSTACK_SIDE,
        ),
        np.nan,
        dtype=np.float32,
    )
    empty_velocity_maps = (
        np.empty(
            (0, velocity_map.shape[0], _ROTATED_SUBSTACK_SIDE, _ROTATED_SUBSTACK_SIDE),
            dtype=np.float32,
        )
        if retain_velocity_maps
        else None
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
    topology = CrossSectionTopology(
        spatial_shape=tuple(vessel.shape),
        frame_count=int(velocity_map.shape[0]),
        labels=branches.labels.copy(),
        branch_ids=branches.branch_ids.copy(),
        section_masks=np.asarray(section_masks, dtype=bool).copy(),
        segment_masks=empty_segment_masks,
        segment_center_xy=np.full(
            (0, settings.ring_count, 2),
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
        profile_rotation_degrees=np.full(
            (settings.ring_count, 0),
            np.nan,
            dtype=np.float32,
        ),
        profile_integration_limits_pixels=np.full(
            (settings.ring_count, 0, 2),
            -1,
            dtype=np.int32,
        ),
        valid_segments=np.zeros((settings.ring_count, 0), dtype=bool),
        ring_settings=settings,
        branch_identity=branches,
        prepared_topology=prepared_topology,
    )
    displacement_results = {
        method: _empty_displacement_result(
            frame_count=velocity_map.shape[0],
            ring_count=settings.ring_count,
        )
        for method in displacement_maps
    }

    return CrossSectionSignalResult(
        velocity=np.full(shape, np.nan, dtype=np.float32),
        safe_velocity=np.full(shape, np.nan, dtype=np.float32),
        velocity_maps_per_segment=empty_velocity_maps,
        velocity_map_segment_indexes=np.empty((0, 2), dtype=np.int32),
        segment_masks=empty_segment_masks,
        labels=branches.labels,
        branch_ids=branches.branch_ids,
        segment_center_xy=np.full(
            (0, settings.ring_count, 2),
            np.nan,
            dtype=np.float32,
        ),
        branch_identity=branches,
        topology=topology,
        displacements=displacement_results,
        velocity_profiles=empty_profiles,
        transverse_velocity_profiles_masked=empty_profiles.copy(),
        longitudinal_velocity_profiles_unmasked=empty_profiles.copy(),
        longitudinal_velocity_profiles_masked=empty_profiles.copy(),
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
    )


def _cross_section_worker_count(
    work_count: int,
    *,
    frame_count: int = 1,
    working_memory_mb: float = 512.0,
) -> int:
    if work_count <= 1 or optional_cupy_backend() is not None:
        return 1
    # Bound each worker by one native and several transformed float32 stacks.
    bytes_per_worker = max(
        1,
        int(frame_count) * _ROTATED_SUBSTACK_SIDE**2 * 4 * 5,
    )
    memory_workers = max(
        1,
        int(float(working_memory_mb) * 1024**2) // bytes_per_worker,
    )
    return min(
        work_count,
        cap_parallel_jobs(_MAX_PARALLEL_CROSS_SECTIONS),
        memory_workers,
    )


def _store_cross_section_measurement(
    buffers: _CrossSectionBuffers,
    circle_index: int,
    branch_index: int,
    measurement: _CrossSectionMeasurement,
    bounds_xyxy: tuple[int, int, int, int],
    *,
    frame_slice: slice = slice(None),
    initialize: bool = True,
) -> None:
    masked = measurement.masked
    buffers.velocity[circle_index, branch_index, frame_slice] = masked.raw
    buffers.safe_velocity[circle_index, branch_index, frame_slice] = masked.safe_velocity
    if buffers.velocity_maps_per_segment is not None:
        map_row = int(buffers.velocity_map_rows[circle_index, branch_index])
        if map_row < 0:
            raise ValueError("Missing compact velocity-map row for valid segment.")
        if measurement.unmasked.rotated_stack is None:
            raise ValueError("Retained velocity maps require a rotated stack.")
        buffers.velocity_maps_per_segment[map_row, frame_slice] = measurement.unmasked.rotated_stack
    buffers.velocity_profiles[circle_index, branch_index, frame_slice] = (
        measurement.unmasked.transverse_profiles
    )
    buffers.transverse_velocity_profiles_masked[circle_index, branch_index, frame_slice] = (
        masked.transverse_profiles
    )
    buffers.longitudinal_velocity_profiles_unmasked[circle_index, branch_index, frame_slice] = (
        measurement.unmasked.longitudinal_profiles
    )
    buffers.longitudinal_velocity_profiles_masked[circle_index, branch_index, frame_slice] = (
        masked.longitudinal_profiles
    )
    if not initialize:
        return
    buffers.segment_masks[circle_index, branch_index] = measurement.rotated_mask
    buffers.profile_sample_count[circle_index, branch_index] = measurement.sample_count
    buffers.profile_rotation_degrees[circle_index, branch_index] = np.float32(
        masked.angle
    )
    buffers.profile_window_bounds_xyxy[circle_index, branch_index] = bounds_xyxy
    buffers.profile_integration_limits_pixels[circle_index, branch_index] = (
        measurement.limits
    )


def _project_displacement_map(
    displacement_map,
    topology: CrossSectionTopology,
    *,
    retain_maps: bool = True,
) -> CrossSectionDisplacementResult:
    buffers = _CrossSectionDisplacementBuffers.allocate(
        frame_count=topology.frame_count,
        ring_count=topology.section_masks.shape[0],
        branch_count=topology.branch_ids.size,
        retain_maps=retain_maps,
    )
    work_items = [
        (int(circle_index), int(branch_index))
        for circle_index, branch_index in np.argwhere(topology.valid_segments)
    ]

    def measure(indexes):
        circle_index, branch_index = indexes
        return _cross_section_displacement_from_topology(
            displacement_map,
            topology,
            circle_index,
            branch_index,
        )

    worker_count = (
        1
        if isinstance(displacement_map, np.memmap)
        else _cross_section_worker_count(len(work_items))
    )
    if worker_count == 1:
        for indexes in work_items:
            _store_displacement_measurement(
                buffers,
                indexes[0],
                indexes[1],
                measure(indexes),
            )
    elif isinstance(displacement_map, np.ndarray) and not isinstance(
        displacement_map, np.memmap
    ):
        with ThreadPoolExecutor(
            max_workers=worker_count,
            thread_name_prefix='cross-section-displacement',
        ) as executor:
            measurements = executor.map(measure, work_items)
            for indexes, measurement in zip(
                work_items,
                measurements,
                strict=True,
            ):
                _store_displacement_measurement(
                    buffers,
                    indexes[0],
                    indexes[1],
                    measurement,
                )
    else:
        _project_lazy_displacement_map(
            buffers,
            displacement_map,
            topology,
            work_items,
            worker_count,
        )
    return _displacement_result_from_buffers(buffers)


def _project_lazy_displacement_map(
    buffers: _CrossSectionDisplacementBuffers,
    displacement_map,
    topology: CrossSectionTopology,
    work_items: list[tuple[int, int]],
    worker_count: int,
) -> None:
    pending: deque[
        tuple[int, int, Future[_CrossSectionDisplacementMeasurement]]
    ] = deque()
    with ThreadPoolExecutor(
        max_workers=worker_count,
        thread_name_prefix='cross-section-displacement',
    ) as executor:
        for circle_index, branch_index in work_items:
            work = _prepare_cross_section_displacement_work(
                displacement_map,
                topology,
                circle_index,
                branch_index,
            )
            pending.append(
                (
                    circle_index,
                    branch_index,
                    executor.submit(_measure_cross_section_displacement, work),
                )
            )
            if len(pending) >= worker_count:
                _store_pending_displacement(buffers, pending.popleft())
        while pending:
            _store_pending_displacement(buffers, pending.popleft())


def _store_pending_displacement(
    buffers: _CrossSectionDisplacementBuffers,
    pending: tuple[
        int,
        int,
        Future[_CrossSectionDisplacementMeasurement],
    ],
) -> None:
    circle_index, branch_index, future = pending
    _store_displacement_measurement(
        buffers,
        circle_index,
        branch_index,
        future.result(),
    )


def _store_displacement_measurement(
    buffers: _CrossSectionDisplacementBuffers,
    circle_index: int,
    branch_index: int,
    measurement: _CrossSectionDisplacementMeasurement,
) -> None:
    index = (circle_index, branch_index)
    buffers.displacement[index] = measurement.waveform.raw
    buffers.safe_displacement[index] = measurement.waveform.safe_displacement
    if buffers.displacement_maps_per_segment is not None:
        buffers.displacement_maps_per_segment[index] = measurement.vectors
    buffers.transverse_displacement_profiles_unmasked[index] = (
        measurement.transverse_profiles_unmasked
    )
    buffers.transverse_displacement_profiles_masked[index] = (
        measurement.transverse_profiles_masked
    )
    buffers.longitudinal_displacement_profiles_unmasked[index] = (
        measurement.longitudinal_profiles_unmasked
    )
    buffers.longitudinal_displacement_profiles_masked[index] = (
        measurement.longitudinal_profiles_masked
    )
    buffers.x_sum_displacement_profile[index] = (
        measurement.summed_displacement_xy[:, 0]
    )
    buffers.y_sum_displacement_profile[index] = (
        measurement.summed_displacement_xy[:, 1]
    )
    buffers.cross_sectional_radial_movement_amplitude[index] = (
        measurement.radial_movement_amplitude
    )
    buffers.cross_sectional_radial_asymmetry_index[index] = (
        measurement.radial_asymmetry_index
    )


def _legacy_topology_from_prepared(
    buffers: _CrossSectionBuffers,
    prepared_topology: PreparedTopology,
    ring_settings: SegmentRingSettings,
    *,
    frame_count: int,
    profile_pixel_size_mm: float,
) -> CrossSectionTopology:
    topology = prepared_topology.topology
    branches = topology.branch_identity
    if branches is None:
        raise ValueError("prepared topology must retain its branch identity result.")

    bounds = topology.window_bounds_xyxy
    limits = buffers.profile_integration_limits_pixels
    valid_segments = (
        topology.valid_segments
        & np.isfinite(prepared_topology.rotation_degrees)
        & np.all(bounds >= 0, axis=-1)
        & (bounds[..., 0] < bounds[..., 1])
        & (bounds[..., 2] < bounds[..., 3])
        & (limits[..., 0] >= 0)
        & (limits[..., 0] <= limits[..., 1])
    )
    return CrossSectionTopology(
        spatial_shape=topology.spatial_shape,
        frame_count=int(frame_count),
        labels=topology.labels.copy(),
        branch_ids=topology.branch_ids.copy(),
        section_masks=topology.annulus_masks.copy(),
        segment_masks=buffers.segment_masks.copy(),
        segment_center_xy=np.transpose(
            topology.segment_centers_xy,
            (1, 0, 2),
        ).copy(),
        profile_window_bounds_xyxy=bounds.copy(),
        profile_window_side_pixels=int(topology.window_side_pixels),
        profile_pixel_size_mm=float(profile_pixel_size_mm),
        profile_rotation_degrees=prepared_topology.rotation_degrees.copy(),
        profile_integration_limits_pixels=limits.copy(),
        valid_segments=valid_segments,
        ring_settings=ring_settings,
        branch_identity=branches,
        prepared_topology=prepared_topology,
    )


def _displacement_result_from_buffers(
    buffers: _CrossSectionDisplacementBuffers,
) -> CrossSectionDisplacementResult:
    return CrossSectionDisplacementResult(
        displacement=buffers.displacement,
        safe_displacement=buffers.safe_displacement,
        displacement_maps_per_segment=buffers.displacement_maps_per_segment,
        transverse_displacement_profiles_unmasked=(
            buffers.transverse_displacement_profiles_unmasked
        ),
        transverse_displacement_profiles_masked=(
            buffers.transverse_displacement_profiles_masked
        ),
        longitudinal_displacement_profiles_unmasked=(
            buffers.longitudinal_displacement_profiles_unmasked
        ),
        longitudinal_displacement_profiles_masked=(
            buffers.longitudinal_displacement_profiles_masked
        ),
        x_sum_displacement_profile=buffers.x_sum_displacement_profile,
        y_sum_displacement_profile=buffers.y_sum_displacement_profile,
        cross_sectional_radial_movement_amplitude=(
            buffers.cross_sectional_radial_movement_amplitude
        ),
        cross_sectional_radial_asymmetry_index=(
            buffers.cross_sectional_radial_asymmetry_index
        ),
    )


def _empty_displacement_result(
    *,
    frame_count: int,
    ring_count: int,
) -> CrossSectionDisplacementResult:
    return _displacement_result_from_buffers(
        _CrossSectionDisplacementBuffers.allocate(
            frame_count=frame_count,
            ring_count=ring_count,
            branch_count=0,
        )
    )


def _interpolated_pixel_size_mm(
    native_pixel_size_mm: float,
    substack_side_pixels: int,
) -> float:
    if substack_side_pixels <= 0:
        return 0.0
    return float(
        native_pixel_size_mm * float(substack_side_pixels) / float(_INTERPOLATED_SUBSTACK_SIDE)
    )


def _cross_section_displacement_from_topology(
    displacement_map,
    topology: CrossSectionTopology,
    circle_index: int,
    branch_index: int,
) -> _CrossSectionDisplacementMeasurement:
    return _measure_cross_section_displacement(
        _prepare_cross_section_displacement_work(
            displacement_map,
            topology,
            circle_index,
            branch_index,
        )
    )


def _prepare_cross_section_displacement_work(
    displacement_map,
    topology: CrossSectionTopology,
    circle_index: int,
    branch_index: int,
) -> _CrossSectionDisplacementWork:
    prepared = topology.prepared_topology
    if prepared is None:
        raise ValueError("Displacement projection requires prepared topology.")
    index = (circle_index, branch_index)
    return _CrossSectionDisplacementWork(
        component_stacks=extract_segment(
            displacement_map, prepared.topology, *index, spatial_axes=(1, 2),
        ),
        angle=float(topology.profile_rotation_degrees[index]),
        rotated_mask=topology.segment_masks[index],
        limits=tuple(int(value) for value in topology.profile_integration_limits_pixels[index]),
    )

def _measure_cross_section_displacement(
    work: _CrossSectionDisplacementWork,
) -> _CrossSectionDisplacementMeasurement:
    frame_count = work.component_stacks.shape[0]
    component_frames = work.component_stacks.reshape(
        frame_count * 2,
        work.component_stacks.shape[-2],
        work.component_stacks.shape[-1],
    )
    rotated = resample_rotate_segment(
        component_frames, work.angle, _INTERPOLATED_SUBSTACK_SIDE,
    ).reshape(frame_count, 2, _ROTATED_SUBSTACK_SIDE, _ROTATED_SUBSTACK_SIDE)
    vectors = _correct_displacement_basis(
        rotated[:, 0],
        rotated[:, 1],
        work.angle,
    )
    profile_mask = dilate_segment_masks(
        work.rotated_mask,
        iterations=_ARTERY_TRANSVERSE_MASK_DILATION_PIXELS,
    )
    radial_amplitude, radial_asymmetry = _cross_sectional_radial_metrics(
        vectors,
        work.rotated_mask,
    )
    (
        transverse_profiles_unmasked,
        transverse_profiles_masked,
        longitudinal_profiles_unmasked,
        longitudinal_profiles_masked,
    ) = _displacement_profiles(vectors, profile_mask)
    c1, c2 = work.limits
    return _CrossSectionDisplacementMeasurement(
        waveform=_displacement_profile_measurement(
            vectors,
            c1,
            c2,
            mask=profile_mask,
        ),
        vectors=vectors,
        transverse_profiles_unmasked=transverse_profiles_unmasked,
        transverse_profiles_masked=transverse_profiles_masked,
        longitudinal_profiles_unmasked=longitudinal_profiles_unmasked,
        longitudinal_profiles_masked=longitudinal_profiles_masked,
        summed_displacement_xy=_nansum_float32(
            vectors,
            axis=(1, 2),
        ),
        radial_movement_amplitude=radial_amplitude,
        radial_asymmetry_index=radial_asymmetry,
    )


def _cross_sectional_radial_metrics(
    vectors: np.ndarray,
    vessel_mask: np.ndarray,
    *,
    epsilon: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Measure radial strength in wall bands extending outside the vessel mask."""

    if vectors.ndim != 4 or vectors.shape[-1] != 2:
        raise ValueError("vectors must have shape (frame, y, x, 2).")
    mask = np.asarray(vessel_mask, dtype=bool)
    if mask.shape != vectors.shape[1:3]:
        raise ValueError("vessel_mask must match the vector field spatial shape.")

    left_region, right_region = _wall_analysis_regions(mask)
    radial_strength = np.abs(vectors[..., 0])
    left = _mean_in_spatial_region(radial_strength, left_region)
    right = _mean_in_spatial_region(radial_strength, right_region)
    amplitude = np.float32(0.5) * (left + right)
    denominator = left + right + np.float32(epsilon)
    asymmetry = np.divide(
        left - right,
        denominator,
        out=np.full_like(amplitude, np.nan),
        where=np.isfinite(denominator),
    )
    return (
        amplitude.astype(np.float32, copy=False),
        asymmetry.astype(np.float32, copy=False),
    )


def _displacement_profiles(
    rotated_vectors: np.ndarray,
    profile_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Project displacement magnitudes using the velocity-profile reductions."""

    magnitude = np.hypot(
        rotated_vectors[..., 0],
        rotated_vectors[..., 1],
    ).astype(np.float32, copy=False)
    transverse_unmasked = calculate_transverse_profiles(magnitude)
    longitudinal_unmasked = calculate_longitudinal_profiles(magnitude)
    masked_magnitude = magnitude.copy()
    masked_magnitude[:, ~np.asarray(profile_mask, dtype=bool)] = np.nan
    transverse_masked = calculate_transverse_profiles(masked_magnitude)
    longitudinal_masked = calculate_longitudinal_profiles(masked_magnitude)
    return (
        transverse_unmasked,
        transverse_masked,
        longitudinal_unmasked,
        longitudinal_masked,
    )


def _wall_analysis_regions(vessel_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Use mask rows to locate walls; return equally scaled exterior bands."""

    height, width = vessel_mask.shape
    left_region = np.zeros((height, width), dtype=bool)
    right_region = np.zeros((height, width), dtype=bool)
    x = np.arange(width)
    for row in np.flatnonzero(np.any(vessel_mask, axis=1)):
        mask_x = np.flatnonzero(vessel_mask[row])
        left_wall = int(mask_x[0])
        right_wall = int(mask_x[-1])
        centerline = np.float32(0.5 * (left_wall + right_wall))
        band_width = max(
            int(np.ceil((right_wall - left_wall + 1) / 2.0)),
            1,
        )
        left_start = max(left_wall - band_width, 0)
        right_stop = min(right_wall + band_width + 1, width)
        left_region[row] = (
            (x >= left_start) & (x <= left_wall) & (x < centerline)
        )
        right_region[row] = (
            (x >= right_wall) & (x < right_stop) & (x > centerline)
        )
    return left_region, right_region


def _mean_in_spatial_region(
    values: np.ndarray,
    region: np.ndarray,
) -> np.ndarray:
    finite = np.isfinite(values) & region[None, ...]
    count = np.sum(finite, axis=(1, 2), dtype=np.int32)
    total = np.sum(
        values,
        axis=(1, 2),
        dtype=np.float32,
        where=finite,
    )
    return np.divide(
        total,
        count,
        out=np.full(values.shape[0], np.nan, dtype=np.float32),
        where=count > 0,
    )


def _correct_displacement_basis(
    rotated_dx: np.ndarray,
    rotated_dy: np.ndarray,
    angle_degrees: float,
) -> np.ndarray:
    cosine = np.float32(special.cosdg(angle_degrees))
    sine = np.float32(special.sindg(angle_degrees))
    local_x = cosine * rotated_dx + sine * rotated_dy
    local_y = -sine * rotated_dx + cosine * rotated_dy
    return np.stack((local_x, local_y), axis=-1).astype(np.float32, copy=False)


def _displacement_profile_measurement(
    rotated_vectors: np.ndarray,
    c1: int,
    c2: int,
    *,
    mask: np.ndarray | None = None,
) -> _DisplacementProfileMeasurement:
    axial_values = rotated_vectors[..., 1]
    finite = np.isfinite(axial_values)
    if mask is not None:
        if mask.shape != axial_values.shape[1:]:
            raise ValueError('displacement profile mask must match the vector field.')
        finite &= mask[None, ...]
    axial_profiles = _nanmean_float32_where(
        axial_values,
        finite,
        axis=1,
    )
    raw = nanmean_float32(axial_profiles[:, c1 : c2 + 1], axis=1)
    safe_displacement = nanmean_float32(axial_profiles, axis=1)
    return _DisplacementProfileMeasurement(
        raw=raw,
        safe_displacement=safe_displacement,
    )


def _nanmean_float32_where(
    values: np.ndarray,
    finite: np.ndarray,
    *,
    axis: int,
) -> np.ndarray:
    count = np.sum(finite, axis=axis, dtype=np.int32)
    total = np.sum(
        values,
        axis=axis,
        dtype=np.float32,
        where=finite,
    )
    return np.divide(
        total,
        count,
        out=np.full_like(total, np.nan, dtype=np.float32),
        where=count > 0,
    )


def _nansum_float32(
    values: np.ndarray,
    *,
    axis: tuple[int, ...],
) -> np.ndarray:
    finite = np.isfinite(values)
    count = np.sum(finite, axis=axis, dtype=np.int32)
    total = np.sum(
        values,
        axis=axis,
        dtype=np.float32,
        where=finite,
    )
    return np.where(
        count > 0,
        total,
        np.float32(np.nan),
    ).astype(np.float32, copy=False)


def _dilate_profile_mask(
    mask: np.ndarray,
    dilation_pixels: int = 0,
) -> np.ndarray:
    """Expand a rotated mask horizontally for profile export only."""

    return dilate_segment_masks(
        mask,
        iterations=dilation_pixels,
        horizontal_only=True,
    )


def _gpu_cross_section_measurement(
    rotated_stack,
    rotated_masked_stack,
    rotated_mask: np.ndarray,
    profile_mask: np.ndarray,
    angle: float,
    *,
    use_dilated_transverse_profile: bool,
    retain_rotated_stack: bool,
    cupy,
) -> _CrossSectionMeasurement:
    """Reduce a resident CUDA segment and transfer only profiles/summaries."""

    transverse = _gpu_nanmean(rotated_stack, axis=1, cupy=cupy)
    longitudinal = _gpu_nanmean(rotated_stack, axis=2, cupy=cupy)
    transverse_masked = _gpu_nanmean(rotated_masked_stack, axis=1, cupy=cupy)
    longitudinal_masked = _gpu_nanmean(
        rotated_masked_stack,
        axis=2,
        cupy=cupy,
    )
    transverse_profile = (
        _gpu_nanmean(
            rotated_stack,
            axis=1,
            spatial_mask=cupy.asarray(profile_mask, dtype=cupy.bool_),
            cupy=cupy,
        )
        if use_dilated_transverse_profile
        else transverse_masked
    )
    rotated_mean = _gpu_nanmean(rotated_stack, axis=0, cupy=cupy)
    rotated_mean_masked = _gpu_nanmean(
        rotated_masked_stack,
        axis=0,
        spatial_mask=cupy.asarray(rotated_mask, dtype=cupy.bool_),
        cupy=cupy,
    )
    (
        transverse_host,
        longitudinal_host,
        transverse_masked_host,
        transverse_profile_host,
        longitudinal_masked_host,
        rotated_mean_host,
        rotated_mean_masked_host,
    ) = (
        cupy.asnumpy(value)
        for value in (
            transverse,
            longitudinal,
            transverse_masked,
            transverse_profile,
            longitudinal_masked,
            rotated_mean,
            rotated_mean_masked,
        )
    )
    c1, c2 = _cross_section_limits(rotated_mean_masked_host)
    retained = cupy.asnumpy(rotated_stack) if retain_rotated_stack else None
    masked_measurement = _profile_measurement_from_profiles(
        transverse_masked_host,
        longitudinal_masked_host,
        angle,
        c1,
        c2,
        rotated_stack=None,
    )
    if use_dilated_transverse_profile:
        masked_measurement = replace(
            masked_measurement,
            transverse_profiles=transverse_profile_host,
        )
    return _CrossSectionMeasurement(
        unmasked=_profile_measurement_from_profiles(
            transverse_host,
            longitudinal_host,
            angle,
            c1,
            c2,
            rotated_stack=retained,
        ),
        masked=masked_measurement,
        rotated_mean=rotated_mean_host,
        rotated_mean_masked=rotated_mean_masked_host,
        rotated_mask=rotated_mask,
        limits=(c1, c2),
        sample_count=_rotated_profile_sample_count(angle),
    )


def _gpu_nanmean(values, *, axis: int, cupy, spatial_mask=None):
    finite = cupy.isfinite(values)
    if spatial_mask is not None:
        finite &= spatial_mask[None, ...]
    counts = cupy.sum(finite, axis=axis, dtype=cupy.int32)
    totals = cupy.sum(
        cupy.where(finite, values, cupy.float32(0.0)),
        axis=axis,
        dtype=cupy.float32,
    )
    nonempty = counts > 0
    safe_counts = cupy.where(nonempty, counts, cupy.int32(1))
    result = cupy.empty(totals.shape, dtype=cupy.float32)
    cupy.divide(totals, safe_counts, out=result)
    result[~nonempty] = cupy.nan
    return result


def _profile_measurement_from_profiles(
    transverse_profiles: np.ndarray,
    longitudinal_profiles: np.ndarray,
    angle: float,
    c1: int,
    c2: int,
    *,
    rotated_stack: np.ndarray | None,
) -> _CrossSectionVelocityMeasurement:
    transverse = np.asarray(transverse_profiles, dtype=np.float32)
    longitudinal = np.asarray(longitudinal_profiles, dtype=np.float32)
    raw = nanmean_float32(transverse[:, c1 : c2 + 1], axis=1)
    return _CrossSectionVelocityMeasurement(
        raw=raw,
        safe_velocity=nanmean_float32(transverse, axis=1),
        transverse_profiles=transverse,
        longitudinal_profiles=longitudinal,
        rotated_stack=rotated_stack,
        angle=float(angle),
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


def _cross_section_limits(image: np.ndarray) -> tuple[int, int]:
    """Use the complete transverse support; fitting belongs to analysis."""

    return 0, max(int(image.shape[1]) - 1, 0)
