"""Generic profile measurements on vessel-aligned signal-map segments."""

from __future__ import annotations

from collections.abc import Callable, Mapping, MutableMapping
from dataclasses import dataclass
from time import perf_counter

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32
from calculations.topology import (
    BranchIdentityResult,
    PreparedTopology,
    SegmentRingSettings,
    TopologyCacheKey,
    dilate_segment_masks,
    longitudinal_profiles,
    prepare_segment_chunks,
    prepare_topologies,
    resolve_segment_rotations,
    transverse_profiles,
)
from runtime_limits import cap_parallel_jobs
from utils.logger import Logger


SegmentObserver = Callable[[int, int, slice, object, np.ndarray], None]
SegmentObserverFactory = Callable[[str, PreparedTopology], SegmentObserver | None]


@dataclass(frozen=True)
class SegmentProfileSettings:
    """Resource and spatial-scale settings for generic profile measurement."""

    pixel_size_mm: float
    submask_size_percentile_kept: float = 0.95
    working_memory_mb: float = 512.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.pixel_size_mm) or self.pixel_size_mm <= 0:
            raise ValueError("pixel_size_mm must be finite and positive.")
        if not np.isfinite(self.working_memory_mb) or self.working_memory_mb <= 0:
            raise ValueError("working_memory_mb must be finite and positive.")
        if not 0 < self.submask_size_percentile_kept <= 1:
            raise ValueError("submask_size_percentile_kept must be in (0, 1].")

    @classmethod
    def from_value(cls, value) -> SegmentProfileSettings:
        if isinstance(value, cls):
            return value
        return cls(
            pixel_size_mm=float(value.pixel_size_mm),
            submask_size_percentile_kept=float(
                value.submask_size_percentile_kept
            ),
            working_memory_mb=float(value.working_memory_mb),
        )


@dataclass(frozen=True)
class SegmentProfileTopology:
    """Topology and profile geometry shared by every measured signal map."""

    spatial_shape: tuple[int, int]
    frame_count: int
    labels: np.ndarray
    branch_ids: np.ndarray
    section_masks: np.ndarray
    segment_masks: np.ndarray
    segment_centers_xy: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_window_side_pixels: int
    profile_pixel_size_mm: float
    profile_rotation_degrees: np.ndarray
    profile_integration_limits_pixels: np.ndarray
    valid_segments: np.ndarray
    ring_settings: SegmentRingSettings
    branch_identity: BranchIdentityResult
    prepared_topology: PreparedTopology


@dataclass(frozen=True, kw_only=True)
class SegmentProfileResult:
    """Signal-neutral segment waveforms, maps, and spatial profiles."""

    projected_signal: np.ndarray
    full_profile_signal: np.ndarray
    segment_maps: np.ndarray | None
    segment_map_indexes: np.ndarray
    segment_masks: np.ndarray
    labels: np.ndarray
    branch_ids: np.ndarray
    segment_centers_xy: np.ndarray
    branch_identity: BranchIdentityResult
    topology: SegmentProfileTopology
    transverse_profiles_unmasked: np.ndarray
    transverse_profiles_masked: np.ndarray
    longitudinal_profiles_unmasked: np.ndarray
    longitudinal_profiles_masked: np.ndarray
    profile_sample_count: np.ndarray
    profile_rotation_degrees: np.ndarray
    rotated_mean_images: np.ndarray
    rotated_mean_images_masked: np.ndarray
    profile_window_bounds_xyxy: np.ndarray
    profile_window_side_pixels: int
    profile_pixel_size_mm: float
    profile_integration_limits_pixels: np.ndarray


@dataclass
class _SegmentProfileBuffers:
    projected_signal: np.ndarray
    full_profile_signal: np.ndarray
    segment_maps: np.ndarray | None
    segment_map_indexes: np.ndarray
    segment_map_rows: np.ndarray
    segment_masks: np.ndarray
    transverse_profiles_unmasked: np.ndarray
    transverse_profiles_masked: np.ndarray
    longitudinal_profiles_unmasked: np.ndarray
    longitudinal_profiles_masked: np.ndarray
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
        canvas_side: int,
        segment_indexes: np.ndarray,
        retain_segment_maps: bool,
    ) -> _SegmentProfileBuffers:
        indexes = np.asarray(segment_indexes, dtype=np.int32).reshape((-1, 2))
        signal_shape = (ring_count, branch_count, frame_count)
        segment_shape = (ring_count, branch_count)
        profile_shape = (*signal_shape, canvas_side)
        image_shape = (*segment_shape, canvas_side, canvas_side)
        map_rows = np.full(segment_shape, -1, dtype=np.int32)
        if retain_segment_maps and indexes.size:
            map_rows[indexes[:, 0], indexes[:, 1]] = np.arange(
                indexes.shape[0], dtype=np.int32
            )

        def filled(shape):
            return np.full(shape, np.nan, dtype=np.float32)

        return cls(
            projected_signal=filled(signal_shape),
            full_profile_signal=filled(signal_shape),
            segment_maps=(
                filled((indexes.shape[0], frame_count, canvas_side, canvas_side))
                if retain_segment_maps
                else None
            ),
            segment_map_indexes=(
                indexes if retain_segment_maps else np.empty((0, 2), dtype=np.int32)
            ),
            segment_map_rows=map_rows,
            segment_masks=np.zeros(image_shape, dtype=bool),
            transverse_profiles_unmasked=filled(profile_shape),
            transverse_profiles_masked=filled(profile_shape),
            longitudinal_profiles_unmasked=filled(profile_shape),
            longitudinal_profiles_masked=filled(profile_shape),
            profile_sample_count=np.zeros(segment_shape, dtype=np.int32),
            profile_rotation_degrees=filled(segment_shape),
            rotated_mean_images=filled(image_shape),
            rotated_mean_images_masked=filled(image_shape),
            profile_window_bounds_xyxy=np.full(
                (*segment_shape, 4), -1, dtype=np.int32
            ),
            profile_integration_limits_pixels=np.full(
                (*segment_shape, 2), -1, dtype=np.int32
            ),
        )


def analyze_segment_profiles(
    signal_map,
    vessel_masks: Mapping[str, object],
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    profile_settings,
    *,
    optic_disc_mask=None,
    source_id: str = "",
    topology_cache: MutableMapping[TopologyCacheKey, object] | None = None,
    retain_segment_maps: bool = False,
    transverse_mask_dilation_pixels: int | Mapping[str, int] = 0,
    prepared_topologies: Mapping[str, PreparedTopology] | None = None,
    transform_mode: str = "fused",
    post_interpolation=None,
    temporal_halo: int = 0,
    scratch_array_count: int | None = None,
    segment_observer_factory: SegmentObserverFactory | None = None,
) -> dict[str, SegmentProfileResult]:
    """Measure profiles from a time-varying scalar map for each vessel mask."""

    settings = SegmentProfileSettings.from_value(profile_settings)
    masks = {
        str(name): np.asarray(mask, dtype=bool)
        for name, mask in vessel_masks.items()
    }
    if not masks:
        return {}
    for mask in masks.values():
        _validate_signal_map(signal_map, mask)

    backend = optional_cupy_backend()
    Logger.log(
        "Segment-profile compute backend: "
        + ("CuPy/CUDA" if backend is not None else "CPU/SciPy")
        + "."
    )
    Logger.log(
        "Segment-profile input: "
        f"map_shape={tuple(int(size) for size in signal_map.shape)}, "
        f"map_type={type(signal_map).__name__}, vessels={tuple(masks)}."
    )
    topology_started = perf_counter()
    if prepared_topologies is None:
        topologies = prepare_topologies(
            masks,
            optic_disc_mask,
            ring_settings,
            source_id=source_id,
            cache=topology_cache,
            optic_disc_center=optic_disc_center,
            window_size_percentile_kept=settings.submask_size_percentile_kept,
        )
    else:
        topologies = dict(prepared_topologies)
        if set(topologies) != set(masks):
            raise ValueError("prepared_topologies must match vessel mask names.")
    topologies = {
        name: resolve_segment_rotations(
            topology,
            signal_map,
            working_memory_mb=settings.working_memory_mb,
        )
        for name, topology in topologies.items()
    }
    Logger.log(
        f"Completed topology preparation in {perf_counter() - topology_started:.2f}s."
    )

    results: dict[str, SegmentProfileResult] = {}
    for name, topology in topologies.items():
        geometry = topology.topology
        Logger.log(
            f"Preparing {name} segments: radii={geometry.annulus_masks.shape[0]}, "
            f"branches={geometry.branch_ids.size}, "
            f"valid_segments={int(np.count_nonzero(geometry.valid_segments))}, "
            f"native_window={geometry.window_side_pixels}px."
        )
        worker_count = _profile_worker_count(
            int(np.count_nonzero(geometry.valid_segments)),
            frame_count=int(signal_map.shape[0]),
            canvas_side=int(topology.rotated_masks.shape[-1]),
            working_memory_mb=settings.working_memory_mb,
        )
        segments = prepare_segment_chunks(
            signal_map,
            topology,
            worker_count=worker_count,
            working_memory_mb=settings.working_memory_mb,
            keep_on_device=backend is not None,
            transform_mode=transform_mode,
            post_interpolation=post_interpolation,
            temporal_halo=temporal_halo,
            scratch_array_count=scratch_array_count,
            include_masked_before_rotation=True,
        )
        observer = (
            segment_observer_factory(name, topology)
            if segment_observer_factory is not None
            else None
        )
        Logger.log(f"Streaming {name} segments into profile measurement.")
        measurement_started = perf_counter()
        if backend is not None:
            backend.cupy.cuda.get_current_stream().synchronize()
        results[name] = _measure_segment_profiles_from_prepared(
            signal_map,
            topology,
            segments,
            ring_settings,
            settings,
            retain_segment_maps=retain_segment_maps,
            segment_observer=observer,
            transverse_mask_dilation_pixels=_dilation_pixels(
                transverse_mask_dilation_pixels,
                name,
            ),
        )
        if backend is not None:
            backend.cupy.cuda.get_current_stream().synchronize()
        Logger.log(
            f"Completed {name} segment profile measurements in "
            f"{perf_counter() - measurement_started:.2f}s."
        )
    return results


def _measure_segment_profiles_from_prepared(
    signal_map,
    prepared_topology: PreparedTopology,
    prepared_segments,
    ring_settings: SegmentRingSettings,
    settings: SegmentProfileSettings,
    *,
    retain_segment_maps: bool,
    segment_observer: SegmentObserver | None,
    transverse_mask_dilation_pixels: int,
) -> SegmentProfileResult:
    geometry = prepared_topology.topology
    branches = geometry.branch_identity
    if branches is None:
        raise ValueError("prepared topology must retain its branch identity result.")

    ring_count = int(geometry.annulus_masks.shape[0])
    branch_count = int(branches.branch_ids.size)
    frame_count = int(signal_map.shape[0])
    canvas_side = int(prepared_topology.rotated_masks.shape[-1])
    interpolated_side = int(prepared_topology.interpolated_masks.shape[-1])
    segment_indexes = np.argwhere(
        geometry.valid_segments
        & np.isfinite(prepared_topology.rotation_degrees)
    )
    buffers = _SegmentProfileBuffers.allocate(
        frame_count=frame_count,
        ring_count=ring_count,
        branch_count=branch_count,
        canvas_side=canvas_side,
        segment_indexes=segment_indexes,
        retain_segment_maps=retain_segment_maps,
    )
    mean_sum = np.zeros(buffers.rotated_mean_images.shape, dtype=np.float64)
    mean_count = np.zeros(buffers.rotated_mean_images.shape, dtype=np.int64)
    masked_sum = np.zeros(buffers.rotated_mean_images.shape, dtype=np.float64)
    masked_count = np.zeros(buffers.rotated_mean_images.shape, dtype=np.int64)
    initialized: set[tuple[int, int]] = set()

    for segment in prepared_segments:
        index = (segment.ring_index, segment.branch_index)
        rotated = segment.rotated
        rotated_masked = segment.rotated_masked
        if rotated_masked is None:
            raise ValueError(
                "segment profile measurement requires the masked segment companion."
            )
        rotated_mask = prepared_topology.rotated_masks[index]
        profile_mask = dilate_segment_masks(
            rotated_mask,
            iterations=transverse_mask_dilation_pixels,
            horizontal_only=True,
        )
        if segment_observer is not None:
            segment_observer(
                index[0],
                index[1],
                segment.frame_slice,
                rotated,
                profile_mask,
            )
        _accumulate_means(
            rotated,
            rotated_masked,
            rotated_mask,
            mean_sum[index],
            mean_count[index],
            masked_sum[index],
            masked_count[index],
        )
        transverse_unmasked = transverse_profiles(rotated)
        transverse_masked = transverse_profiles(rotated_masked)
        if transverse_mask_dilation_pixels > 0:
            transverse_masked_for_output = transverse_profiles(
                rotated,
                profile_mask,
            )
        else:
            transverse_masked_for_output = transverse_masked
        longitudinal_unmasked = longitudinal_profiles(rotated)
        longitudinal_masked = longitudinal_profiles(rotated_masked)
        frame_slice = segment.frame_slice
        masked_signal = _to_numpy(_profile_mean(transverse_masked))
        buffers.projected_signal[index][frame_slice] = masked_signal
        buffers.full_profile_signal[index][frame_slice] = masked_signal
        buffers.transverse_profiles_unmasked[index][frame_slice] = _to_numpy(
            transverse_unmasked
        )
        buffers.transverse_profiles_masked[index][frame_slice] = _to_numpy(
            transverse_masked_for_output
        )
        buffers.longitudinal_profiles_unmasked[index][frame_slice] = _to_numpy(
            longitudinal_unmasked
        )
        buffers.longitudinal_profiles_masked[index][frame_slice] = _to_numpy(
            longitudinal_masked
        )
        if buffers.segment_maps is not None:
            map_row = int(buffers.segment_map_rows[index])
            if map_row < 0:
                raise ValueError("Missing compact segment-map row for valid segment.")
            buffers.segment_maps[map_row, frame_slice] = _to_numpy(rotated)
        if index not in initialized:
            buffers.segment_masks[index] = rotated_mask
            buffers.profile_sample_count[index] = _rotated_profile_sample_count(
                float(prepared_topology.rotation_degrees[index]),
                interpolated_side,
                canvas_side,
            )
            buffers.profile_rotation_degrees[index] = np.float32(
                prepared_topology.rotation_degrees[index]
            )
            buffers.profile_window_bounds_xyxy[index] = geometry.window_bounds_xyxy[
                index
            ]
            buffers.profile_integration_limits_pixels[index] = (0, canvas_side - 1)
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
    profile_pixel_size_mm = _interpolated_pixel_size_mm(
        settings.pixel_size_mm,
        geometry.window_side_pixels,
        interpolated_side,
    )
    bounds = buffers.profile_window_bounds_xyxy
    limits = buffers.profile_integration_limits_pixels
    valid_segments = (
        geometry.valid_segments
        & np.isfinite(buffers.profile_rotation_degrees)
        & np.all(bounds >= 0, axis=-1)
        & (bounds[..., 0] < bounds[..., 1])
        & (bounds[..., 2] < bounds[..., 3])
        & (limits[..., 0] >= 0)
        & (limits[..., 0] <= limits[..., 1])
    )
    topology = SegmentProfileTopology(
        spatial_shape=geometry.spatial_shape,
        frame_count=frame_count,
        labels=geometry.labels.copy(),
        branch_ids=geometry.branch_ids.copy(),
        section_masks=geometry.annulus_masks.copy(),
        segment_masks=buffers.segment_masks.copy(),
        segment_centers_xy=geometry.segment_centers_xy.copy(),
        profile_window_bounds_xyxy=bounds.copy(),
        profile_window_side_pixels=int(geometry.window_side_pixels),
        profile_pixel_size_mm=profile_pixel_size_mm,
        profile_rotation_degrees=buffers.profile_rotation_degrees.copy(),
        profile_integration_limits_pixels=limits.copy(),
        valid_segments=valid_segments,
        ring_settings=ring_settings,
        branch_identity=branches,
        prepared_topology=prepared_topology,
    )
    return SegmentProfileResult(
        projected_signal=buffers.projected_signal,
        full_profile_signal=buffers.full_profile_signal,
        segment_maps=buffers.segment_maps,
        segment_map_indexes=buffers.segment_map_indexes,
        segment_masks=buffers.segment_masks,
        labels=branches.labels,
        branch_ids=branches.branch_ids,
        segment_centers_xy=geometry.segment_centers_xy.copy(),
        branch_identity=branches,
        topology=topology,
        transverse_profiles_unmasked=buffers.transverse_profiles_unmasked,
        transverse_profiles_masked=buffers.transverse_profiles_masked,
        longitudinal_profiles_unmasked=buffers.longitudinal_profiles_unmasked,
        longitudinal_profiles_masked=buffers.longitudinal_profiles_masked,
        profile_sample_count=buffers.profile_sample_count,
        profile_rotation_degrees=buffers.profile_rotation_degrees,
        rotated_mean_images=buffers.rotated_mean_images,
        rotated_mean_images_masked=buffers.rotated_mean_images_masked,
        profile_window_bounds_xyxy=buffers.profile_window_bounds_xyxy,
        profile_window_side_pixels=int(geometry.window_side_pixels),
        profile_pixel_size_mm=profile_pixel_size_mm,
        profile_integration_limits_pixels=buffers.profile_integration_limits_pixels,
    )


def _validate_signal_map(signal_map, vessel_mask: np.ndarray) -> None:
    if vessel_mask.ndim != 2:
        raise ValueError(
            f"vessel_mask must have shape (y, x), got {vessel_mask.shape!r}."
        )
    shape = getattr(signal_map, "shape", None)
    if shape is None or len(shape) != 3:
        raise ValueError(
            f"signal_map must have shape (frame, y, x), got {shape!r}."
        )
    if any(int(size) <= 0 for size in shape):
        raise ValueError("signal_map axes must be nonempty.")
    if tuple(shape[1:]) != tuple(vessel_mask.shape):
        raise ValueError(
            "signal_map spatial shape must match vessel_mask: "
            f"{tuple(shape[1:])!r} != {tuple(vessel_mask.shape)!r}."
        )


def _accumulate_means(
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
        chunk_total = cupy.asnumpy(
            cupy.sum(
                cupy.where(finite, rotated, cupy.float32(0.0)),
                axis=0,
                dtype=cupy.float64,
            )
        )
        chunk_count = cupy.asnumpy(cupy.sum(finite, axis=0, dtype=cupy.int64))
        masked_finite = cupy.isfinite(rotated_masked) & cupy.asarray(
            rotated_mask, dtype=cupy.bool_
        )[None]
        chunk_masked_total = cupy.asnumpy(
            cupy.sum(
                cupy.where(masked_finite, rotated_masked, cupy.float32(0.0)),
                axis=0,
                dtype=cupy.float64,
            )
        )
        chunk_masked_count = cupy.asnumpy(
            cupy.sum(masked_finite, axis=0, dtype=cupy.int64)
        )
    else:
        values = np.asarray(rotated, dtype=np.float32)
        finite = np.isfinite(values)
        chunk_total = np.sum(values, axis=0, dtype=np.float64, where=finite)
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


def _profile_mean(profiles):
    backend = optional_cupy_backend()
    if backend is not None and isinstance(profiles, backend.cupy.ndarray):
        finite = backend.cupy.isfinite(profiles)
        count = backend.cupy.sum(finite, axis=-1, dtype=backend.cupy.int32)
        total = backend.cupy.sum(
            backend.cupy.where(finite, profiles, backend.cupy.float32(0.0)),
            axis=-1,
            dtype=backend.cupy.float32,
        )
        return backend.cupy.where(
            count > 0,
            total / backend.cupy.maximum(count, 1),
            backend.cupy.float32(np.nan),
        ).astype(backend.cupy.float32, copy=False)
    return nanmean_float32(np.asarray(profiles, dtype=np.float32), axis=-1)


def _to_numpy(values) -> np.ndarray:
    backend = optional_cupy_backend()
    if backend is not None and isinstance(values, backend.cupy.ndarray):
        return backend.cupy.asnumpy(values)
    return np.asarray(values)


def _rotated_profile_sample_count(
    angle_degrees: float,
    interpolated_side: int,
    canvas_side: int,
) -> int:
    if not np.isfinite(angle_degrees):
        return 0
    radians = np.deg2rad(float(angle_degrees))
    scale = abs(float(np.cos(radians))) + abs(float(np.sin(radians)))
    count = int(np.floor(interpolated_side * scale + 0.5))
    return min(max(count, interpolated_side), canvas_side)


def _interpolated_pixel_size_mm(
    native_pixel_size_mm: float,
    native_side_pixels: int,
    interpolated_side_pixels: int,
) -> float:
    if native_side_pixels <= 0 or interpolated_side_pixels <= 0:
        return 0.0
    return float(
        native_pixel_size_mm
        * float(native_side_pixels)
        / float(interpolated_side_pixels)
    )


def _profile_worker_count(
    work_count: int,
    *,
    frame_count: int,
    canvas_side: int,
    working_memory_mb: float,
) -> int:
    if work_count <= 1 or optional_cupy_backend() is not None:
        return 1
    bytes_per_worker = max(1, frame_count * canvas_side**2 * 4 * 5)
    memory_workers = max(
        1,
        int(float(working_memory_mb) * 1024**2) // bytes_per_worker,
    )
    return min(work_count, cap_parallel_jobs(8), memory_workers)


def _dilation_pixels(value: int | Mapping[str, int], vessel_name: str) -> int:
    if isinstance(value, Mapping):
        return int(value.get(vessel_name, 0))
    return int(value)


__all__ = [
    "SegmentProfileResult",
    "SegmentProfileSettings",
    "SegmentProfileTopology",
    "analyze_segment_profiles",
]
