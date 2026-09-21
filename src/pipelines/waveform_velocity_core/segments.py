"""Prepare and analyze velocity segments for configured vessel classes."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from dataclasses import replace
from time import perf_counter

import numpy as np
from scipy.signal import resample

from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    _ARTERY_TRANSVERSE_MASK_DILATION_PIXELS,
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    _cross_section_worker_count,
    _generate_cross_section_signals_from_prepared,
    _validate_velocity_map,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32, next_power_of_two
from calculations.topology import (
    SegmentRingSettings,
    TopologyCacheKey,
    prepare_segment_chunks,
    prepare_topologies,
    resolve_segment_rotations,
)
from pipelines.displacement_map.segments import analyze_displacement_segments
from utils.logger import Logger

_FFT_PROFILE_X_BATCH = 32


class _VelocityProfileFftAccumulator:
    """Accumulate FFT profiles while each rotated segment is still resident."""

    def __init__(
        self,
        *,
        frame_count: int,
        ring_count: int,
        branch_count: int,
        canvas_side: int,
        cycle_boundary_indexes,
        index_base: int,
    ) -> None:
        self.boundaries = normalize_cycle_boundaries(
            cycle_boundary_indexes,
            frame_count,
            index_base=index_base,
        )
        self.time_count = next_power_of_two(
            int(np.max(np.diff(self.boundaries)))
        )
        beat_count = self.boundaries.size - 1
        output_shape = (
            canvas_side,
            self.time_count,
            beat_count,
            branch_count,
            ring_count,
        )
        self.unmasked = np.full(output_shape, np.nan, dtype=np.float32)
        self.masked = np.full(output_shape, np.nan, dtype=np.float32)
        self.elapsed_seconds = 0.0
        self._pending_key = None
        self._pending_parts = []

    def observe(
        self,
        ring_index: int,
        branch_index: int,
        frame_slice,
        rotated_chunk=None,
        profile_mask: np.ndarray | None = None,
    ) -> None:
        """Consume one ordered chunk while retaining at most one beat."""

        started = perf_counter()
        if profile_mask is None:
            profile_mask = rotated_chunk
            rotated_chunk = frame_slice
            frame_slice = slice(0, int(rotated_chunk.shape[0]))
        chunk_start = int(frame_slice.start or 0)
        chunk_stop = int(frame_slice.stop)
        mask = np.asarray(profile_mask, dtype=bool)
        if rotated_chunk.ndim != 3 or mask.shape != rotated_chunk.shape[1:]:
            raise ValueError(
                "FFT profile input must contain a (frame, y, x) chunk and "
                "matching (y, x) mask."
            )
        for beat_index in range(self.boundaries.size - 1):
            beat_start = int(self.boundaries[beat_index])
            beat_stop = int(self.boundaries[beat_index + 1]) + 1
            overlap_start = max(chunk_start, beat_start)
            overlap_stop = min(chunk_stop, beat_stop)
            if overlap_start >= overlap_stop:
                continue
            part = rotated_chunk[
                overlap_start - chunk_start : overlap_stop - chunk_start
            ].copy()
            key = (ring_index, branch_index, beat_index)
            if self._pending_key not in (None, key):
                raise RuntimeError("FFT chunks arrived out of segment/frame order.")
            self._pending_key = key
            self._pending_parts.append(part)
            buffered = sum(int(value.shape[0]) for value in self._pending_parts)
            expected = beat_stop - beat_start
            if buffered == expected:
                backend = optional_cupy_backend()
                if (
                    backend is not None
                    and isinstance(self._pending_parts[0], backend.cupy.ndarray)
                ):
                    beat = backend.cupy.concatenate(self._pending_parts, axis=0)
                else:
                    beat = np.concatenate(self._pending_parts, axis=0)
                self._write_beat(
                    ring_index,
                    branch_index,
                    beat_index,
                    beat,
                    mask,
                )
                self._pending_key = None
                self._pending_parts = []
            elif buffered > expected:
                raise RuntimeError("FFT chunk buffering exceeded the current beat.")
        self.elapsed_seconds += perf_counter() - started

    def _write_beat(
        self,
        ring_index: int,
        branch_index: int,
        beat_index: int,
        beat,
        mask: np.ndarray,
    ) -> None:
        backend = optional_cupy_backend()
        if backend is not None and isinstance(beat, backend.cupy.ndarray):
            try:
                self._write_beat_gpu(
                    ring_index,
                    branch_index,
                    beat_index,
                    beat,
                    mask,
                    backend.cupy,
                )
                backend.cupy.cuda.get_current_stream().synchronize()
                return
            except Exception as exc:
                Logger.log_debug(
                    "CuPy velocity FFT profiles failed; using CPU fallback: "
                    f"{type(exc).__name__}: {exc}"
                )
                beat = backend.cupy.asnumpy(beat)
        self._write_beat_cpu(
            ring_index,
            branch_index,
            beat_index,
            np.asarray(beat, dtype=np.float32),
            mask,
        )

    def _write_beat_cpu(
        self,
        ring_index: int,
        branch_index: int,
        beat_index: int,
        beat: np.ndarray,
        mask: np.ndarray,
    ) -> None:
        for x_start in range(0, beat.shape[2], _FFT_PROFILE_X_BATCH):
            x_stop = min(x_start + _FFT_PROFILE_X_BATCH, beat.shape[2])
            interpolated = _interpft_stack_axis0(
                beat[:, :, x_start:x_stop],
                self.time_count + 1,
            )[:-1]
            magnitude = np.abs(
                np.fft.fft(interpolated, axis=0)
            ).astype(np.float32, copy=False)
            output_slice = (
                slice(x_start, x_stop),
                slice(None),
                beat_index,
                branch_index,
                ring_index,
            )
            self.unmasked[output_slice] = nanmean_float32(
                magnitude,
                axis=1,
            ).T
            self.masked[output_slice] = nanmean_float32(
                np.where(
                    mask[None, :, x_start:x_stop],
                    magnitude,
                    np.float32(np.nan),
                ),
                axis=1,
            ).T

    def _write_beat_gpu(
        self,
        ring_index: int,
        branch_index: int,
        beat_index: int,
        beat,
        mask: np.ndarray,
        cupy,
    ) -> None:
        pixel_count = int(beat.shape[1] * beat.shape[2])
        flattened = beat.reshape((beat.shape[0], pixel_count))
        active_pixels = cupy.any(cupy.isfinite(flattened), axis=0)
        if int(cupy.count_nonzero(active_pixels)) == 0:
            return
        interpolated = _gpu_fourier_resample_axis0(
            flattened[:, active_pixels],
            self.time_count + 1,
            cupy,
        )[:-1]
        magnitude_active = cupy.abs(
            cupy.fft.fft(interpolated, axis=0)
        ).astype(cupy.float32, copy=False)
        magnitude = cupy.full(
            (self.time_count, pixel_count),
            cupy.nan,
            dtype=cupy.float32,
        )
        magnitude[:, active_pixels] = magnitude_active
        magnitude = magnitude.reshape(
            (self.time_count, beat.shape[1], beat.shape[2])
        )
        gpu_mask = cupy.asarray(mask, dtype=cupy.bool_)
        unmasked = _gpu_nanmean_axis1(magnitude, cupy)
        masked = _gpu_nanmean_axis1(magnitude, cupy, mask=gpu_mask)
        output_slice = (
            slice(None),
            slice(None),
            beat_index,
            branch_index,
            ring_index,
        )
        self.unmasked[output_slice] = cupy.asnumpy(unmasked).T
        self.masked[output_slice] = cupy.asnumpy(masked).T

def _interpft_stack_axis0(values: np.ndarray, target_length: int) -> np.ndarray:
    """Vectorized Fourier interpolation of one segment stack's frame axis."""

    source = np.asarray(values, dtype=np.float32)
    source_length = int(source.shape[0])
    if source_length == 0:
        raise ValueError("interpft requires a non-empty frame axis.")
    if target_length <= 0:
        raise ValueError("interpft target_length must be positive.")
    if target_length == source_length:
        return source.copy()

    spatial_shape = source.shape[1:]
    flattened = source.reshape(source_length, -1)
    active_pixels = np.any(np.isfinite(flattened), axis=0)
    interpolated = np.full(
        (int(target_length), flattened.shape[1]),
        np.nan,
        dtype=np.float32,
    )
    if np.any(active_pixels):
        interpolated[:, active_pixels] = resample(
            flattened[:, active_pixels],
            int(target_length),
            axis=0,
        ).astype(np.float32, copy=False)
    return interpolated.reshape(int(target_length), *spatial_shape)


def _gpu_fourier_resample_axis0(values, target_length: int, cupy):
    """CuPy equivalent of SciPy's real Fourier resampling along axis zero."""

    source_length = int(values.shape[0])
    if source_length == 0:
        raise ValueError("interpft requires a non-empty frame axis.")
    if target_length <= 0:
        raise ValueError("interpft target_length must be positive.")
    if target_length == source_length:
        return values.copy()

    spectrum = cupy.fft.rfft(values, axis=0)
    output_spectrum = cupy.zeros(
        (target_length // 2 + 1, *spectrum.shape[1:]),
        dtype=spectrum.dtype,
    )
    common_length = min(source_length, int(target_length))
    nyquist_stop = common_length // 2 + 1
    output_spectrum[:nyquist_stop] = spectrum[:nyquist_stop]
    if common_length % 2 == 0:
        nyquist_index = common_length // 2
        if target_length < source_length:
            output_spectrum[nyquist_index] *= cupy.float32(2.0)
        else:
            output_spectrum[nyquist_index] *= cupy.float32(0.5)
    result = cupy.fft.irfft(output_spectrum, n=target_length, axis=0)
    result *= cupy.float32(target_length / source_length)
    return result


def _gpu_nanmean_axis1(values, cupy, *, mask=None):
    finite = cupy.isfinite(values)
    if mask is not None:
        finite &= mask[None, ...]
    counts = cupy.sum(finite, axis=1)
    totals = cupy.sum(
        cupy.where(finite, values, cupy.float32(0.0)),
        axis=1,
        dtype=cupy.float32,
    )
    safe_counts = counts.copy()
    safe_counts[safe_counts == 0] = 1
    result = cupy.divide(totals, safe_counts)
    result[counts == 0] = cupy.nan
    return result.astype(cupy.float32, copy=False)


def analyze_velocity_segments(
    velocity_map,
    vessel_masks: Mapping[str, object],
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    optic_disc_mask=None,
    source_id: str = "",
    topology_cache: MutableMapping[TopologyCacheKey, object] | None = None,
    retain_velocity_maps: bool = False,
    cycle_boundary_indexes=None,
    velocity_profile_fft: bool = False,
    index_base: int = 0,
    transverse_mask_dilation_pixels: int | None = None,
    displacement_maps_by_vessel: Mapping[str, Mapping[str, object]] | None = None,
    retain_displacement_maps: bool = False,
    prepared_topologies: Mapping[str, object] | None = None,
    transform_mode: str = "fused",
    post_interpolation=None,
    temporal_halo: int = 0,
    scratch_array_count: int | None = None,
) -> dict[str, CrossSectionSignalResult]:
    """Analyze velocity-map segments for every named vessel mask.

    Vessel names are retained as result keys. All selected vessels contribute
    to one shared segment-window size before their maps are prepared.
    """

    masks = {
        str(name): np.asarray(mask, dtype=bool)
        for name, mask in vessel_masks.items()
    }
    if not masks:
        return {}
    if velocity_profile_fft and cycle_boundary_indexes is None:
        raise ValueError(
            "cycle_boundary_indexes are required for velocity FFT profiles."
        )
    for mask in masks.values():
        _validate_velocity_map(velocity_map, mask)

    backend = optional_cupy_backend()
    Logger.log(
        "Cross-section compute backend: "
        + ("CuPy/CUDA" if backend is not None else "CPU/SciPy")
        + "."
    )
    Logger.log(
        "Segment input: "
        f"map_shape={tuple(int(size) for size in velocity_map.shape)}, "
        f"map_type={type(velocity_map).__name__}, vessels={tuple(masks)}."
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
            window_size_percentile_kept=(
                cross_section_settings.submask_size_percentile_kept
            ),
        )
    else:
        topologies = dict(prepared_topologies)
        if set(topologies) != set(masks):
            raise ValueError("prepared_topologies must match vessel mask names.")
    topologies = {
        name: resolve_segment_rotations(
            topology,
            velocity_map,
            working_memory_mb=cross_section_settings.working_memory_mb,
        )
        for name, topology in topologies.items()
    }
    Logger.log(
        f"Completed topology preparation in {perf_counter() - topology_started:.2f}s."
    )

    results: dict[str, CrossSectionSignalResult] = {}
    for name, topology in topologies.items():
        geometry = topology.topology
        Logger.log(
            f"Preparing {name} segments: radii={geometry.annulus_masks.shape[0]}, "
            f"branches={geometry.branch_ids.size}, "
            f"valid_segments={int(np.count_nonzero(geometry.valid_segments))}, "
            f"native_window={geometry.window_side_pixels}px."
        )
        worker_count = _cross_section_worker_count(
            int(np.count_nonzero(geometry.valid_segments)),
            frame_count=int(velocity_map.shape[0]),
            working_memory_mb=float(cross_section_settings.working_memory_mb),
        )
        segments = prepare_segment_chunks(
            velocity_map,
            topology,
            worker_count=worker_count,
            working_memory_mb=cross_section_settings.working_memory_mb,
            keep_on_device=backend is not None,
            transform_mode=transform_mode,
            post_interpolation=post_interpolation,
            temporal_halo=temporal_halo,
            scratch_array_count=scratch_array_count,
            include_masked_before_rotation=True,
        )
        Logger.log(f"Streaming {name} segments into profile measurement.")
        fft_profiles = (
            _VelocityProfileFftAccumulator(
                frame_count=int(velocity_map.shape[0]),
                ring_count=int(geometry.annulus_masks.shape[0]),
                branch_count=int(geometry.branch_ids.size),
                canvas_side=int(topology.rotated_masks.shape[-1]),
                cycle_boundary_indexes=cycle_boundary_indexes,
                index_base=index_base,
            )
            if velocity_profile_fft
            else None
        )
        measurement_started = perf_counter()
        if backend is not None:
            backend.cupy.cuda.get_current_stream().synchronize()
        result = _generate_cross_section_signals_from_prepared(
            velocity_map,
            topology,
            segments,
            ring_settings,
            cross_section_settings,
            retain_velocity_maps=retain_velocity_maps,
            segment_observer=(fft_profiles.observe if fft_profiles else None),
            transverse_mask_dilation_pixels=(
                _legacy_profile_dilation_pixels(name)
                if transverse_mask_dilation_pixels is None
                else int(transverse_mask_dilation_pixels)
            ),
            displacement_maps=None,
            retain_displacement_maps=retain_displacement_maps,
        )
        displacement_results = analyze_displacement_segments(
            (displacement_maps_by_vessel or {}).get(name, {}),
            topology,
            retain_maps=retain_displacement_maps,
            working_memory_mb=cross_section_settings.working_memory_mb,
        )
        results[name] = replace(
            result,
            displacements=displacement_results,
            transverse_velocity_fft_profiles_unmasked=(
                fft_profiles.unmasked if fft_profiles else None
            ),
            transverse_velocity_fft_profiles_masked=(
                fft_profiles.masked if fft_profiles else None
            ),
        )
        if backend is not None:
            backend.cupy.cuda.get_current_stream().synchronize()
        if fft_profiles is not None:
            Logger.log(
                f"Completed {name} optional velocity-profile FFT in "
                f"{fft_profiles.elapsed_seconds:.2f}s."
            )
        Logger.log(
            f"Completed {name} segment profile measurements in "
            f"{perf_counter() - measurement_started:.2f}s."
        )
    return results


def _legacy_profile_dilation_pixels(vessel_name: str) -> int:
    """Retain the historical artery-only transverse profile expansion."""

    return (
        _ARTERY_TRANSVERSE_MASK_DILATION_PIXELS
        if str(vessel_name).lower() == "artery"
        else 0
    )
