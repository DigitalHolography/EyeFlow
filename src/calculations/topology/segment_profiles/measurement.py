"""Generic profile measurements on vessel-aligned signal-map segments."""

from __future__ import annotations

from collections.abc import Callable, Mapping, MutableMapping
from time import perf_counter

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from runtime_limits import cap_parallel_jobs
from utils.logger import Logger

from ..cache import TopologyCacheKey
from ..geometry import AnnulusGeometry
from ..optic_disc import OpticDisc
from ..transforms import dilate_segment_masks
from ..workflow import (
    PreparedTopology,
    prepare_segment_chunks,
    prepare_topologies,
    resolve_segment_rotations,
)
from .accumulator import SegmentProfileAccumulator
from .models import SegmentProfileResult, SegmentProfileSettings

SegmentObserver = Callable[[int, int, slice, object, np.ndarray], None]
SegmentObserverFactory = Callable[[str, PreparedTopology], SegmentObserver | None]


def analyze_segment_profiles(
    signal_map,
    vessel_masks: Mapping[str, object],
    optic_disc: OpticDisc,
    ring_settings: AnnulusGeometry,
    profile_settings,
    *,
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
            optic_disc,
            ring_settings,
            source_id=source_id,
            cache=topology_cache,
            window_size_percentile_kept=settings.submask_size_percentile_kept,
        )
        topologies = {
            name: resolve_segment_rotations(
                topology,
                signal_map,
                working_memory_mb=settings.working_memory_mb,
            )
            for name, topology in topologies.items()
        }
    else:
        topologies = dict(prepared_topologies)
        if set(topologies) != set(masks):
            raise ValueError("prepared_topologies must match vessel mask names.")
    Logger.log(
        f"Completed topology preparation in {perf_counter() - topology_started:.2f}s."
    )

    results: dict[str, SegmentProfileResult] = {}
    for name, topology in topologies.items():
        geometry = topology.native
        Logger.log(
            f"Preparing {name} segments: radii={geometry.annulus_masks.shape[0]}, "
            f"branches={geometry.branch_ids.size}, "
            f"valid_segments={int(np.count_nonzero(geometry.valid_segments))}, "
            f"native_window={geometry.window_side_pixels}px."
        )
        worker_count = _profile_worker_count(
            int(np.count_nonzero(geometry.valid_segments)),
            keep_on_device=backend is not None,
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
    settings: SegmentProfileSettings,
    *,
    retain_segment_maps: bool,
    segment_observer: SegmentObserver | None,
    transverse_mask_dilation_pixels: int,
) -> SegmentProfileResult:
    geometry = prepared_topology.native
    frame_count = int(signal_map.shape[0])
    interpolated_side = int(prepared_topology.interpolated_masks.shape[-1])
    accumulator = SegmentProfileAccumulator(
        prepared_topology,
        frame_count=frame_count,
        retain_segment_maps=retain_segment_maps,
    )

    for segment in prepared_segments:
        index = (segment.ring_index, segment.branch_index)
        rotated = segment.rotated
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
        accumulator.add_chunk(segment, profile_mask)

    sample_spacing_mm = _interpolated_pixel_size_mm(
        settings.pixel_size_mm,
        geometry.window_side_pixels,
        interpolated_side,
    )
    return accumulator.finish(
        sample_spacing_mm=sample_spacing_mm,
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
    keep_on_device: bool,
) -> int:
    """Return the concurrency cap; chunk planning enforces the memory bound."""

    if work_count <= 1 or keep_on_device:
        return 1
    return min(work_count, cap_parallel_jobs(8))


def _dilation_pixels(value: int | Mapping[str, int], vessel_name: str) -> int:
    if isinstance(value, Mapping):
        return int(value.get(vessel_name, 0))
    return int(value)


__all__ = [
    "analyze_segment_profiles",
]
