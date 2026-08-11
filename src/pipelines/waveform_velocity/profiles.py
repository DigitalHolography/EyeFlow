"""Minimal HDF5 export for AngioEye transverse velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from input_output.schema import CrossSectionProfileOutputPaths, EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue


def pack_cross_section_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    schema = _resolve_output_paths(output_paths)
    metrics = _pack_vessel_profiles(
        schema.artery_cross_section_profiles,
        artery_segments,
        cycle_boundary_indexes,
        index_base=index_base,
    )
    metrics.update(
        _pack_vessel_profiles(
            schema.vein_cross_section_profiles,
            vein_segments,
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return metrics


def _pack_vessel_profiles(
    paths: CrossSectionProfileOutputPaths,
    segments,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> dict[str, object]:
    return {
        paths.velocity_profile: _profile_dataset(
            np.asarray(segments.velocity_profiles, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
        ),
        paths.velocity_profile_masked: _profile_dataset(
            np.asarray(segments.velocity_profiles_masked, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
        ),
    }


def _profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> DatasetValue:
    if profiles.ndim != 4:
        raise ValueError(
            "profile arrays must have shape "
            "(radius, branch, frame, transverse_sample)."
        )
    profiles_per_beat = interpolate_velocity_profiles_per_beat(
        profiles,
        cycle_boundary_indexes,
        index_base=index_base,
    )
    return DatasetValue(
        data=profiles_per_beat,
        attrs={
            "unit": "mm/s",
            "dimDesc": ["x", "time", "beat", "branch", "radius"],
        },
        h5_options=_profile_h5_options(profiles_per_beat.shape),
    )


def _profile_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    """Use lossless compression with chunks aligned to one segment profile."""
    options: dict[str, object] = {
        "compression": "gzip",
        "compression_opts": 4,
        "shuffle": True,
    }
    if len(shape) != 5 or not all(shape):
        return options

    x_count, time_count = shape[:2]
    target_elements = (1024 * 1024) // np.dtype(np.float32).itemsize
    time_chunk = min(time_count, max(target_elements // x_count, 1))
    options["chunks"] = (x_count, time_chunk, 1, 1, 1)
    return options


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
