"""Minimal HDF5 export for AngioEye transverse velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from input_output.schema import CrossSectionProfileOutputPaths, EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue

from pipelines.waveform_velocity_core.per_beat_outputs import metric_value


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
    raw_profiles, raw_coordinates = _profile_arrays(segments, "raw_profile")
    interpolated_profiles, interpolated_coordinates = _profile_arrays(
        segments,
        "interpolated_profile",
    )
    return {
        **_pack_profile_variant(
            paths.raw_profile,
            raw_profiles,
            raw_coordinates,
            cycle_boundary_indexes,
            index_base=index_base,
        ),
        **_pack_profile_variant(
            paths.interpolated_profile,
            interpolated_profiles,
            interpolated_coordinates,
            cycle_boundary_indexes,
            index_base=index_base,
        ),
    }


def _profile_arrays(segments, variant: str) -> tuple[np.ndarray, np.ndarray]:
    profile_data = getattr(segments, variant, None)
    if profile_data is not None:
        return (
            np.asarray(profile_data.velocity, dtype=np.float32),
            np.asarray(profile_data.x_micrometers, dtype=np.float32),
        )

    if variant == "raw_profile":
        return (
            np.asarray(segments.velocity_profiles, dtype=np.float32),
            np.asarray(segments.profile_x_micrometers, dtype=np.float32),
        )
    return (
        np.asarray(segments.centered_velocity_profiles, dtype=np.float32),
        np.asarray(segments.centered_profile_x_micrometers, dtype=np.float32),
    )


def _pack_profile_variant(
    paths,
    profiles: np.ndarray,
    coordinates: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> dict[str, object]:
    if profiles.ndim != 4:
        raise ValueError(
            "profile arrays must have shape "
            "(radius, branch, frame, transverse_sample)."
        )
    if coordinates.ndim == 1:
        if coordinates.shape != (profiles.shape[3],):
            raise ValueError(
                "raw profile coordinates must have shape "
                f"({profiles.shape[3]},), got {coordinates.shape}."
            )
        coordinates_by_x = coordinates
        coordinate_dim_desc = ("x",)
    else:
        expected_coordinates = (
            profiles.shape[0],
            profiles.shape[1],
            profiles.shape[3],
        )
        if coordinates.shape != expected_coordinates:
            raise ValueError(
                "interpolated profile coordinates must have shape "
                f"{expected_coordinates}, got {coordinates.shape}."
            )
        coordinates_by_x = np.transpose(coordinates, (2, 1, 0))
        coordinate_dim_desc = ("x", "branch", "radius")
    profiles_per_beat = interpolate_velocity_profiles_per_beat(
        profiles,
        cycle_boundary_indexes,
        index_base=index_base,
    )
    return {
        paths.velocity_profile: DatasetValue(
            data=profiles_per_beat,
            attrs={
                "unit": "mm/s",
                "dimDesc": ["x", "time", "beat", "branch", "radius"],
            },
            h5_options=_profile_h5_options(profiles_per_beat.shape),
        ),
        paths.transverse_coordinate_micrometers: metric_value(
            coordinates_by_x,
            unit="um",
            dim_desc=coordinate_dim_desc,
        ),
    }


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
