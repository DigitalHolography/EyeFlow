"""Minimal HDF5 export for AngioEye transverse velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from input_output.schema import CrossSectionProfileOutputPaths, EyeFlowOutputPaths

from .outputs import metric_value


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
    profiles = np.asarray(segments.centered_velocity_profiles, dtype=np.float32)
    coordinates = np.asarray(
        segments.centered_profile_x_micrometers,
        dtype=np.float32,
    )
    if profiles.ndim != 4:
        raise ValueError(
            "centered profiles must have shape "
            "(radius, branch, frame, transverse_sample)."
        )
    expected_coordinates = (
        profiles.shape[0],
        profiles.shape[1],
        profiles.shape[3],
    )
    if coordinates.shape != expected_coordinates:
        raise ValueError(
            "centered profile coordinates must have shape "
            f"{expected_coordinates}, got {coordinates.shape}."
        )
    profiles_per_beat = interpolate_velocity_profiles_per_beat(
        profiles,
        cycle_boundary_indexes,
        index_base=index_base,
    )
    coordinates_by_x = np.transpose(coordinates, (2, 1, 0))
    return {
        paths.velocity_profile: metric_value(
            profiles_per_beat,
            unit="mm/s",
            dim_desc=("x", "time", "beat", "branch", "radius"),
        ),
        paths.transverse_coordinate_micrometers: metric_value(
            coordinates_by_x,
            unit="um",
            dim_desc=("x", "branch", "radius"),
        ),
    }


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
