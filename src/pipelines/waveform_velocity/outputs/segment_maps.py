"""Serialization of per-beat segment velocity maps and masks."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue
from pipelines.retinal_velocity.models import RetinalVelocity
from pipelines.retinal_velocity.semantics import velocity_dataset_attrs

from .paths import resolve_output_paths


def pack_segment_map_outputs(
    artery_segments,
    vein_segments,
    artery_velocity_maps_per_beat: np.ndarray | None,
    vein_velocity_maps_per_beat: np.ndarray | None,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    velocity_analysis: RetinalVelocity | None = None,
) -> dict[str, object]:
    """Pack prepared per-beat maps and masks for artery and vein segments."""

    schema = resolve_output_paths(output_paths)
    velocity_attrs = velocity_dataset_attrs(velocity_analysis)
    outputs = _pack_vessel_segment_maps(
        artery_segments,
        artery_velocity_maps_per_beat,
        schema.artery_segments,
        velocity_attrs,
    )
    outputs.update(
        _pack_vessel_segment_maps(
            vein_segments,
            vein_velocity_maps_per_beat,
            schema.vein_segments,
            velocity_attrs,
        )
    )
    return outputs


def _pack_vessel_segment_maps(
    segments,
    velocity_maps_per_beat: np.ndarray | None,
    paths,
    velocity_attrs: dict[str, object],
) -> dict[str, object]:
    if segments is None:
        if velocity_maps_per_beat is not None:
            raise ValueError(
                "velocity_maps_per_beat must be None when segments is None."
            )
        return {}

    outputs: dict[str, object] = {}
    if paths.velocity_map_per_segment is not None:
        if velocity_maps_per_beat is None:
            raise ValueError(
                "velocity_maps_per_beat is required when map output is enabled."
            )
        maps_per_beat = np.asarray(velocity_maps_per_beat, dtype=np.float32)
        if maps_per_beat.ndim != 6:
            raise ValueError(
                "velocity_maps_per_beat must have shape "
                "(x, y, time, beat, branch, radius)."
            )
        masks = np.asarray(segments.profile.topology.rotated_masks, dtype=bool)
        expected_mask_shape = (
            maps_per_beat.shape[5],
            maps_per_beat.shape[4],
            maps_per_beat.shape[1],
            maps_per_beat.shape[0],
        )
        if masks.shape != expected_mask_shape:
            raise ValueError(
                "segment_masks must have shape (radius, branch, y, x) "
                "matching velocity_maps_per_beat."
            )
        outputs[paths.velocity_map_per_segment] = DatasetValue(
            data=maps_per_beat,
            attrs={
                **velocity_attrs,
                "dimDesc": ["x", "y", "time", "beat", "branch", "radius"],
                "coordinate_system": "rotated_segment_pixel",
                "mask_output": paths.segments,
            },
            h5_options=_velocity_map_h5_options(maps_per_beat.shape),
        )

    if paths.segments is not None:
        masks = np.asarray(segments.profile.topology.rotated_masks, dtype=bool)
        if masks.ndim != 4:
            raise ValueError("segment masks must have shape (radius, branch, y, x).")
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


def _velocity_map_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    options: dict[str, object] = {"compression": "lzf", "shuffle": True}
    if len(shape) in (6, 7) and all(shape):
        options["chunks"] = (
            (shape[0], shape[1], 1, 1, 1, 1)
            if len(shape) == 6
            else (shape[0], shape[1], 1, 1, 1, 1, shape[-1])
        )
    return options


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


__all__ = ["pack_segment_map_outputs"]
