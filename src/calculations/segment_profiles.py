"""Generic profile measurements on vessel-aligned signal-map segments."""

from __future__ import annotations

from collections.abc import Callable, Mapping, MutableMapping
from dataclasses import dataclass
from time import perf_counter
from typing import Literal

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32
from calculations.topology import (
    AnnulusGeometry,
    OpticDisc,
    PreparedTopology,
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


@dataclass(frozen=True, slots=True)
class MaskedArrays:
    """Matched unmasked and masked forms of one measured array."""

    unmasked: np.ndarray
    masked: np.ndarray

    def __post_init__(self) -> None:
        unmasked = np.asarray(self.unmasked)
        masked = np.asarray(self.masked)
        if unmasked.shape != masked.shape:
            raise ValueError(
                "masked and unmasked arrays must have the same shape: "
                f"{masked.shape!r} != {unmasked.shape!r}."
            )
        object.__setattr__(self, "unmasked", unmasked)
        object.__setattr__(self, "masked", masked)

    def select(self, *, masked: bool) -> np.ndarray:
        return self.masked if masked else self.unmasked


@dataclass(frozen=True, slots=True)
class CompactSegmentMaps:
    """Maps retained only for valid segments, plus their grid indexes."""

    values: np.ndarray
    indexes: np.ndarray

    def __post_init__(self) -> None:
        values = np.asarray(self.values)
        indexes = np.asarray(self.indexes)
        if indexes.ndim != 2 or indexes.shape[1] != 2:
            raise ValueError("segment map indexes must have shape (segment, 2).")
        if not np.issubdtype(indexes.dtype, np.integer):
            raise ValueError("segment map indexes must contain integers.")
        if values.ndim != 4:
            raise ValueError(
                "compact segment maps must have shape (segment, frame, y, x)."
            )
        if values.shape[0] != indexes.shape[0]:
            raise ValueError("segment map values and indexes must have equal row counts.")
        if indexes.shape[0] != np.unique(indexes, axis=0).shape[0]:
            raise ValueError("segment map indexes must be unique.")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "indexes", indexes.astype(np.int32, copy=False))

    def row_for(self, ring_index: int, branch_index: int) -> int | None:
        matches = np.flatnonzero(
            (self.indexes[:, 0] == int(ring_index))
            & (self.indexes[:, 1] == int(branch_index))
        )
        return None if matches.size == 0 else int(matches[0])

    def to_dense(self, segment_shape: tuple[int, int]) -> np.ndarray:
        """Expand compact maps to ``(ring, branch, frame, y, x)``."""

        shape = tuple(int(size) for size in segment_shape)
        if len(shape) != 2 or any(size < 0 for size in shape):
            raise ValueError("segment_shape must contain two non-negative sizes.")
        if self.indexes.size and (
            np.any(self.indexes < 0)
            or np.any(self.indexes[:, 0] >= shape[0])
            or np.any(self.indexes[:, 1] >= shape[1])
        ):
            raise ValueError("segment map indexes fall outside segment_shape.")
        dense = np.full((*shape, *self.values.shape[1:]), np.nan, dtype=self.values.dtype)
        if self.indexes.size:
            dense[self.indexes[:, 0], self.indexes[:, 1]] = self.values
        return dense


@dataclass(frozen=True, slots=True, kw_only=True)
class SegmentProfileResult:
    """Measurements for every segment in one prepared vessel topology."""

    topology: PreparedTopology
    segment_signal: np.ndarray
    transverse: MaskedArrays
    longitudinal: MaskedArrays
    mean_images: MaskedArrays
    sample_spacing_mm: float
    maps: CompactSegmentMaps | None = None

    def __post_init__(self) -> None:
        segment_shape = self.topology.segment_shape
        profile_side = self.topology.profile_side_pixels
        if self.segment_signal.ndim != 3 or self.segment_signal.shape[:2] != segment_shape:
            raise ValueError(
                "segment_signal must have shape (ring, branch, frame) matching topology."
            )
        expected_prefix = self.segment_signal.shape
        for name, pair in (
            ("transverse", self.transverse),
            ("longitudinal", self.longitudinal),
        ):
            if pair.unmasked.ndim != 4 or pair.unmasked.shape[:3] != expected_prefix:
                raise ValueError(
                    f"{name} profiles must have shape (ring, branch, frame, sample)."
                )
            if pair.unmasked.shape[-1] != profile_side:
                raise ValueError(f"{name} profiles must match the prepared profile size.")
        if self.mean_images.unmasked.ndim != 4:
            raise ValueError("mean images must have shape (ring, branch, y, x).")
        if self.mean_images.unmasked.shape[:2] != segment_shape:
            raise ValueError("mean images must match the topology segment grid.")
        if self.mean_images.unmasked.shape[-2:] != (profile_side, profile_side):
            raise ValueError("mean images must match the prepared profile size.")
        if not np.isfinite(self.sample_spacing_mm) or self.sample_spacing_mm < 0:
            raise ValueError("sample_spacing_mm must be finite and non-negative.")
        if self.maps is not None:
            if self.maps.values.shape[1:] != (
                self.frame_count,
                profile_side,
                profile_side,
            ):
                raise ValueError("retained maps must match the frame and profile sizes.")
            if not np.array_equal(self.maps.indexes, self.topology.valid_indexes()):
                raise ValueError("retained map indexes must match valid topology segments.")

    def profile(
        self,
        direction: Literal["transverse", "longitudinal"],
        *,
        masked: bool,
    ) -> np.ndarray:
        if direction == "transverse":
            return self.transverse.select(masked=masked)
        if direction == "longitudinal":
            return self.longitudinal.select(masked=masked)
        raise ValueError("direction must be 'transverse' or 'longitudinal'.")

    def require_maps(self) -> CompactSegmentMaps:
        if self.maps is None:
            raise RuntimeError(
                "Per-segment maps were not retained; request them during profile analysis."
            )
        return self.maps

    @property
    def frame_count(self) -> int:
        return int(self.segment_signal.shape[-1])

    @property
    def segment_shape(self) -> tuple[int, int]:
        return self.topology.segment_shape


class _SegmentProfileAccumulator:
    """Allocate, populate, and finalize one vessel's profile measurements."""

    def __init__(
        self,
        topology: PreparedTopology,
        *,
        frame_count: int,
        retain_segment_maps: bool,
    ) -> None:
        self.topology = topology
        self.frame_count = int(frame_count)
        segment_shape = topology.segment_shape
        canvas_side = topology.profile_side_pixels
        signal_shape = (*segment_shape, self.frame_count)
        profile_shape = (*signal_shape, canvas_side)
        image_shape = (*segment_shape, canvas_side, canvas_side)

        def filled(shape):
            return np.full(shape, np.nan, dtype=np.float32)

        self.segment_signal = filled(signal_shape)
        self.transverse = MaskedArrays(filled(profile_shape), filled(profile_shape))
        self.longitudinal = MaskedArrays(filled(profile_shape), filled(profile_shape))
        self.mean_images = MaskedArrays(filled(image_shape), filled(image_shape))
        self._mean_sum = np.zeros(image_shape, dtype=np.float64)
        self._mean_count = np.zeros(image_shape, dtype=np.int64)
        self._masked_sum = np.zeros(image_shape, dtype=np.float64)
        self._masked_count = np.zeros(image_shape, dtype=np.int64)

        indexes = topology.valid_indexes().reshape((-1, 2))
        self.maps = (
            CompactSegmentMaps(
                values=filled(
                    (indexes.shape[0], self.frame_count, canvas_side, canvas_side)
                ),
                indexes=indexes,
            )
            if retain_segment_maps
            else None
        )
        self._map_rows = np.full(segment_shape, -1, dtype=np.int32)
        if self.maps is not None and indexes.size:
            self._map_rows[indexes[:, 0], indexes[:, 1]] = np.arange(
                indexes.shape[0], dtype=np.int32
            )

    def add_chunk(self, segment, profile_mask: np.ndarray) -> None:
        index = (segment.ring_index, segment.branch_index)
        rotated = segment.rotated
        rotated_masked = segment.rotated_masked
        if rotated_masked is None:
            raise ValueError(
                "segment profile measurement requires the masked segment companion."
            )
        rotated_mask = self.topology.rotated_masks[index]
        _accumulate_means(
            rotated,
            rotated_masked,
            rotated_mask,
            self._mean_sum[index],
            self._mean_count[index],
            self._masked_sum[index],
            self._masked_count[index],
        )

        transverse_unmasked = transverse_profiles(rotated)
        transverse_masked = transverse_profiles(rotated_masked)
        transverse_masked_for_output = (
            transverse_profiles(rotated, profile_mask)
            if not np.array_equal(profile_mask, rotated_mask)
            else transverse_masked
        )
        frame_slice = segment.frame_slice
        self.segment_signal[index][frame_slice] = _to_numpy(
            _profile_mean(transverse_masked)
        )
        self.transverse.unmasked[index][frame_slice] = _to_numpy(
            transverse_unmasked
        )
        self.transverse.masked[index][frame_slice] = _to_numpy(
            transverse_masked_for_output
        )
        self.longitudinal.unmasked[index][frame_slice] = _to_numpy(
            longitudinal_profiles(rotated)
        )
        self.longitudinal.masked[index][frame_slice] = _to_numpy(
            longitudinal_profiles(rotated_masked)
        )
        if self.maps is not None:
            map_row = int(self._map_rows[index])
            if map_row < 0:
                raise ValueError("Missing compact segment-map row for valid segment.")
            self.maps.values[map_row, frame_slice] = _to_numpy(rotated)

    def finish(self, *, sample_spacing_mm: float) -> SegmentProfileResult:
        np.divide(
            self._mean_sum,
            self._mean_count,
            out=self.mean_images.unmasked,
            where=self._mean_count > 0,
        )
        np.divide(
            self._masked_sum,
            self._masked_count,
            out=self.mean_images.masked,
            where=self._masked_count > 0,
        )
        return SegmentProfileResult(
            topology=self.topology,
            segment_signal=self.segment_signal,
            transverse=self.transverse,
            longitudinal=self.longitudinal,
            mean_images=self.mean_images,
            sample_spacing_mm=sample_spacing_mm,
            maps=self.maps,
        )


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
    accumulator = _SegmentProfileAccumulator(
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
    "CompactSegmentMaps",
    "MaskedArrays",
    "SegmentProfileResult",
    "SegmentProfileSettings",
    "analyze_segment_profiles",
]
