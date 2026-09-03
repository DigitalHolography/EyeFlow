"""Per-beat rotated velocity maps and segment-mask output packing."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.signal import resample

from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.math import next_power_of_two
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue
from runtime_limits import cap_parallel_jobs


_MAX_PARALLEL_SEGMENT_INTERPOLATIONS = 8
_DISPLACEMENT_MAP_ROOT = "Processing/DisplacementMapPerSegment"


def pack_segment_map_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    """Pack rotated maps and masks for artery and vein segments."""
    schema = _resolve_output_paths(output_paths)
    outputs = _pack_vessel_segment_maps(
        artery_segments,
        schema.artery_segments,
        cycle_boundary_indexes,
        index_base=index_base,
    )
    outputs.update(
        _pack_vessel_segment_maps(
            vein_segments,
            schema.vein_segments,
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return outputs


def pack_displacement_segment_map_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    outputs = _pack_vessel_displacement_maps(
        artery_segments,
        "Artery",
        cycle_boundary_indexes,
        index_base=index_base,
    )
    outputs.update(
        _pack_vessel_displacement_maps(
            vein_segments,
            "Vein",
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return outputs


def _pack_vessel_displacement_maps(
    segments,
    vessel_name: str,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> dict[str, object]:
    if segments is None:
        return {}

    outputs: dict[str, object] = {}
    displacement_results = getattr(segments, "displacements", {})
    for raw_method, displacement in sorted(displacement_results.items()):
        method = _hdf_method_name(raw_method)
        displacement_maps_per_beat = np.stack(
            [
                interpolate_velocity_maps_per_beat(
                    displacement.displacement_maps_per_segment[
                        ..., component_index
                    ],
                    cycle_boundary_indexes,
                    index_base=index_base,
                )
                for component_index in range(2)
            ],
            axis=-1,
        )
        outputs[f"{_DISPLACEMENT_MAP_ROOT}/{method}/{vessel_name}"] = (
            DatasetValue(
                data=displacement_maps_per_beat,
                attrs={
                    "unit": "pixels",
                    "dimDesc": [
                        "x",
                        "y",
                        "time",
                        "beat",
                        "branch",
                        "radius",
                        "displacement_orientation",
                    ],
                    "coordinate_system": "rotated_segment_pixel",
                    "components": ["local_x", "local_y"],
                    "component_basis": "rotated_segment_local",
                },
                h5_options=_velocity_map_h5_options(
                    displacement_maps_per_beat.shape
                ),
            )
        )
    return outputs


def interpolate_velocity_maps_per_beat(
    velocity_maps: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> np.ndarray:
    """Interpolate maps to ``(x, y, time, beat, branch, radius)``."""
    maps = np.asarray(velocity_maps, dtype=np.float32)
    if maps.ndim != 5:
        raise ValueError(
            "velocity maps must have shape (radius, branch, frame, y, x)."
        )

    radius_count, branch_count, frame_count, y_count, x_count = maps.shape
    boundaries = normalize_cycle_boundaries(
        cycle_boundary_indexes,
        frame_count,
        index_base=index_base,
    )
    beat_count = boundaries.size - 1
    time_count = next_power_of_two(int(np.max(np.diff(boundaries))))
    output = np.full(
        (
            x_count,
            y_count,
            time_count,
            beat_count,
            branch_count,
            radius_count,
        ),
        np.nan,
        dtype=np.float32,
    )

    segment_indexes = [
        (radius_index, branch_index)
        for radius_index in range(radius_count)
        for branch_index in range(branch_count)
    ]

    def interpolate_segment(segment_index: tuple[int, int]) -> None:
        radius_index, branch_index = segment_index
        for beat_index in range(beat_count):
            start = int(boundaries[beat_index])
            stop = int(boundaries[beat_index + 1]) + 1
            interpolated = _interpft_maps_axis0(
                maps[radius_index, branch_index, start:stop],
                time_count + 1,
            )[:-1]
            output[
                :, :, :, beat_index, branch_index, radius_index
            ] = interpolated.transpose(2, 1, 0)

    worker_count = _segment_map_worker_count(len(segment_indexes))
    if worker_count == 1:
        for segment_index in segment_indexes:
            interpolate_segment(segment_index)
    else:
        with ThreadPoolExecutor(
            max_workers=worker_count,
            thread_name_prefix="segment-map",
        ) as executor:
            for _ in executor.map(interpolate_segment, segment_indexes):
                pass
    return output


def _segment_map_worker_count(segment_count: int) -> int:
    if segment_count <= 1:
        return 1
    return min(
        segment_count,
        cap_parallel_jobs(_MAX_PARALLEL_SEGMENT_INTERPOLATIONS),
    )


def _pack_vessel_segment_maps(
    segments,
    paths,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> dict[str, object]:
    if segments is None:
        return {}

    outputs: dict[str, object] = {}
    if paths.velocity_map_per_segment is not None:
        maps_per_beat = interpolate_velocity_maps_per_beat(
            segments.velocity_maps_per_segment,
            cycle_boundary_indexes,
            index_base=index_base,
        )
        outputs[paths.velocity_map_per_segment] = DatasetValue(
            data=maps_per_beat,
            attrs={
                "unit": "mm/s",
                "dimDesc": ["x", "y", "time", "beat", "branch", "radius"],
                "coordinate_system": "rotated_segment_pixel",
                "mask_output": paths.segments,
            },
            h5_options=_velocity_map_h5_options(maps_per_beat.shape),
        )

    if paths.segments is not None:
        masks = np.asarray(segments.segment_masks, dtype=bool)
        if masks.ndim != 4:
            raise ValueError(
                "segment masks must have shape (radius, branch, y, x)."
            )
        serialized_masks = masks.transpose(3, 2, 1, 0)
        outputs[paths.segments] = DatasetValue(
            data=serialized_masks,
            attrs={
                "dimDesc": ["x", "y", "branch", "radius"],
                "coordinate_system": "rotated_segment_pixel",
            },
            h5_options=_segment_mask_h5_options(serialized_masks.shape),
        )
    return outputs


def _interpft_maps_axis0(values: np.ndarray, target_length: int) -> np.ndarray:
    """Vectorized equivalent of ``interpft_real`` along the frame axis."""
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


def _velocity_map_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    options: dict[str, object] = {
        "compression": "lzf",
        "shuffle": True,
    }
    if len(shape) in (6, 7) and all(shape):
        options["chunks"] = (
            (shape[0], shape[1], 1, 1, 1, 1)
            if len(shape) == 6
            else (shape[0], shape[1], 1, 1, 1, 1, shape[-1])
        )
    return options


def _hdf_method_name(value: object) -> str:
    method = str(value).strip()
    if not method or "/" in method:
        raise ValueError(
            "Displacement registration method names must be non-empty HDF5 path segments."
        )
    return method


def _segment_mask_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    options: dict[str, object] = {
        "dtype": np.bool_,
        "compression": "gzip",
        "compression_opts": 4,
        "shuffle": True,
    }
    if len(shape) == 4 and all(shape):
        options["chunks"] = (shape[0], shape[1], 1, 1)
    return options


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
