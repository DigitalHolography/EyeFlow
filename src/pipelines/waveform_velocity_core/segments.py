"""Prepare and analyze velocity segments for configured vessel classes."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from dataclasses import replace
from time import perf_counter

import numpy as np
from scipy.signal import resample

from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32, next_power_of_two
from calculations.topology import (
    SegmentRingSettings,
    TopologyCacheKey,
    prepare_segments,
    prepare_topologies,
)
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    _generate_cross_section_signals_from_prepared,
    _validate_velocity_map,
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

    def observe(
        self,
        ring_index: int,
        branch_index: int,
        rotated_stack: np.ndarray,
        profile_mask: np.ndarray,
    ) -> None:
        stack = np.asarray(rotated_stack, dtype=np.float32)
        mask = np.asarray(profile_mask, dtype=bool)
        if stack.ndim != 3 or mask.shape != stack.shape[1:]:
            raise ValueError(
                "FFT profile input must contain a (frame, y, x) stack and "
                "matching (y, x) mask."
            )
        if stack.shape[0] <= int(self.boundaries[-1]):
            raise ValueError("FFT profile boundaries exceed the segment stack.")

        for beat_index in range(self.boundaries.size - 1):
            start = int(self.boundaries[beat_index])
            stop = int(self.boundaries[beat_index + 1]) + 1
            for x_start in range(0, stack.shape[2], _FFT_PROFILE_X_BATCH):
                x_stop = min(x_start + _FFT_PROFILE_X_BATCH, stack.shape[2])
                interpolated = _interpft_stack_axis0(
                    stack[start:stop, :, x_start:x_stop],
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
                del interpolated, magnitude


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
        segments = prepare_segments(velocity_map, topology)
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
        result = _generate_cross_section_signals_from_prepared(
            velocity_map,
            topology,
            segments,
            ring_settings,
            cross_section_settings,
            retain_velocity_maps=retain_velocity_maps,
            segment_observer=(fft_profiles.observe if fft_profiles else None),
        )
        results[name] = replace(
            result,
            transverse_velocity_fft_profiles_unmasked=(
                fft_profiles.unmasked if fft_profiles else None
            ),
            transverse_velocity_fft_profiles_masked=(
                fft_profiles.masked if fft_profiles else None
            ),
        )
        Logger.log(
            f"Completed {name} segment profile measurements in "
            f"{perf_counter() - measurement_started:.2f}s."
        )
    return results
