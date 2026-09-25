"""Prepare and analyze velocity segments for configured vessel classes."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from time import perf_counter

import numpy as np
from scipy.signal import resample

from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    _ARTERY_TRANSVERSE_MASK_DILATION_PIXELS,
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    CrossSectionTopology,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32, next_power_of_two
from calculations.segment_profiles import SegmentProfileResult, analyze_segment_profiles
from calculations.topology import (
    AnnulusGeometry,
    OpticDisc,
    TopologyCacheKey,
)
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


def analyze_velocity_segment_profiles(
    velocity_map,
    vessel_masks: Mapping[str, object],
    optic_disc: OpticDisc,
    ring_settings: AnnulusGeometry,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    source_id: str = "",
    topology_cache: MutableMapping[TopologyCacheKey, object] | None = None,
    retain_velocity_maps: bool = False,
    cycle_boundary_indexes=None,
    velocity_profile_fft: bool = False,
    index_base: int = 0,
    transverse_mask_dilation_pixels: int | None = None,
    prepared_topologies: Mapping[str, object] | None = None,
    transform_mode: str = "fused",
    post_interpolation=None,
    temporal_halo: int = 0,
    scratch_array_count: int | None = None,
) -> dict[str, CrossSectionSignalResult]:
    """Measure velocity profiles and attach velocity-only optional products."""

    if not vessel_masks:
        return {}
    if velocity_profile_fft and cycle_boundary_indexes is None:
        raise ValueError(
            "cycle_boundary_indexes are required for velocity FFT profiles."
    )
    fft_profiles: dict[str, _VelocityProfileFftAccumulator] = {}

    def segment_observer_factory(name, topology):
        geometry = topology.topology
        if not velocity_profile_fft:
            return None
        accumulator = _VelocityProfileFftAccumulator(
            frame_count=int(velocity_map.shape[0]),
            ring_count=int(geometry.annulus_masks.shape[0]),
            branch_count=int(geometry.branch_ids.size),
            canvas_side=int(topology.rotated_masks.shape[-1]),
            cycle_boundary_indexes=cycle_boundary_indexes,
            index_base=index_base,
        )
        fft_profiles[name] = accumulator
        return accumulator.observe

    dilation_pixels = (
        {
            str(name): _legacy_profile_dilation_pixels(str(name))
            for name in vessel_masks
        }
        if transverse_mask_dilation_pixels is None
        else int(transverse_mask_dilation_pixels)
    )
    profile_results = analyze_segment_profiles(
        velocity_map,
        vessel_masks,
        optic_disc,
        ring_settings,
        cross_section_settings,
        source_id=source_id,
        topology_cache=topology_cache,
        retain_segment_maps=retain_velocity_maps,
        transverse_mask_dilation_pixels=dilation_pixels,
        prepared_topologies=prepared_topologies,
        transform_mode=transform_mode,
        post_interpolation=post_interpolation,
        temporal_halo=temporal_halo,
        scratch_array_count=scratch_array_count,
        segment_observer_factory=segment_observer_factory,
    )

    results: dict[str, CrossSectionSignalResult] = {}
    for name, profile_result in profile_results.items():
        accumulator = fft_profiles.get(name)
        results[name] = _velocity_result(
            profile_result,
            fft_profiles=accumulator,
        )
        if accumulator is not None:
            Logger.log(
                f"Completed {name} optional velocity-profile FFT in "
                f"{accumulator.elapsed_seconds:.2f}s."
            )
    return results


def _velocity_result(
    profiles: SegmentProfileResult,
    *,
    fft_profiles: _VelocityProfileFftAccumulator | None,
) -> CrossSectionSignalResult:
    """Adapt a neutral profile result to the established velocity result schema."""

    profile_topology = profiles.topology
    legacy_centers = np.transpose(
        profile_topology.segment_centers_xy,
        (1, 0, 2),
    ).copy()
    topology = CrossSectionTopology(
        spatial_shape=profile_topology.spatial_shape,
        optic_disc_center_xy=profile_topology.optic_disc_center_xy,
        frame_count=profile_topology.frame_count,
        labels=profile_topology.labels,
        branch_ids=profile_topology.branch_ids,
        segment_masks=profile_topology.segment_masks,
        segment_center_xy=legacy_centers,
        profile_window_bounds_xyxy=profile_topology.profile_window_bounds_xyxy,
        profile_window_side_pixels=profile_topology.profile_window_side_pixels,
        profile_pixel_size_mm=profile_topology.profile_pixel_size_mm,
        profile_rotation_degrees=profile_topology.profile_rotation_degrees,
        profile_integration_limits_pixels=(
            profile_topology.profile_integration_limits_pixels
        ),
        valid_segments=profile_topology.valid_segments,
        ring_settings=profile_topology.ring_settings,
        branch_identity=profile_topology.branch_identity,
        prepared_topology=profile_topology.prepared_topology,
    )
    return CrossSectionSignalResult(
        velocity=profiles.projected_signal,
        safe_velocity=profiles.full_profile_signal,
        velocity_maps_per_segment=profiles.segment_maps,
        velocity_map_segment_indexes=profiles.segment_map_indexes,
        segment_masks=profiles.segment_masks,
        labels=profiles.labels,
        branch_ids=profiles.branch_ids,
        segment_center_xy=legacy_centers,
        branch_identity=profiles.branch_identity,
        topology=topology,
        displacements={},
        velocity_profiles=profiles.transverse_profiles_unmasked,
        transverse_velocity_profiles_masked=profiles.transverse_profiles_masked,
        longitudinal_velocity_profiles_unmasked=(
            profiles.longitudinal_profiles_unmasked
        ),
        longitudinal_velocity_profiles_masked=profiles.longitudinal_profiles_masked,
        profile_sample_count=profiles.profile_sample_count,
        profile_rotation_degrees=profiles.profile_rotation_degrees,
        rotated_mean_images=profiles.rotated_mean_images,
        rotated_mean_images_masked=profiles.rotated_mean_images_masked,
        profile_window_bounds_xyxy=profiles.profile_window_bounds_xyxy,
        profile_window_side_pixels=profiles.profile_window_side_pixels,
        profile_pixel_size_mm=profiles.profile_pixel_size_mm,
        profile_integration_limits_pixels=profiles.profile_integration_limits_pixels,
        transverse_velocity_fft_profiles_unmasked=(
            fft_profiles.unmasked if fft_profiles else None
        ),
        transverse_velocity_fft_profiles_masked=(
            fft_profiles.masked if fft_profiles else None
        ),
    )


def _legacy_profile_dilation_pixels(vessel_name: str) -> int:
    """Retain the historical artery-only transverse profile expansion."""

    return (
        _ARTERY_TRANSVERSE_MASK_DILATION_PIXELS
        if str(vessel_name).lower() == "artery"
        else 0
    )
