"""Minimal HDF5 export for transverse and longitudinal velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from input_output.schema import CrossSectionProfileOutputPaths, EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue
from pipelines.waveform_velocity_core.dopplerview.outputs import metric_value


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
    raw_coordinates = np.asarray(
        segments.profile_x_micrometers,
        dtype=np.float32,
    )
    interpolated_coordinates = np.asarray(
        segments.centered_profile_x_micrometers,
        dtype=np.float32,
    )
    return {
        paths.raw_profile.velocity_profile: _profile_dataset(
            np.asarray(segments.velocity_profiles, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
        ),
        paths.raw_profile.transverse_coordinate_micrometers: _coordinate_dataset(
            raw_coordinates,
            spatial_axis="x",
        ),
        paths.raw_profile.transverse_velocity_profile_masked: _profile_dataset(
            np.asarray(
                segments.transverse_velocity_profiles_masked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
        ),
        paths.raw_profile.longitudinal_velocity_profile_masked: _profile_dataset(
            np.asarray(
                segments.longitudinal_velocity_profiles_masked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
        ),
        paths.raw_profile.transverse_velocity_profile_masked_centroid: (
            _profile_dataset(
                np.asarray(
                    segments.transverse_velocity_profile_masked_centroids,
                    dtype=np.float32,
                ),
                cycle_boundary_indexes,
                index_base=index_base,
                spatial_axis="y",
                unit="pixel",
            )
        ),
        paths.raw_profile.longitudinal_coordinate_micrometers: _coordinate_dataset(
            raw_coordinates,
            spatial_axis="y",
        ),
        paths.interpolated_profile.velocity_profile: _profile_dataset(
            np.asarray(segments.centered_velocity_profiles, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
        ),
        paths.interpolated_profile.transverse_coordinate_micrometers: (
            _coordinate_dataset(interpolated_coordinates, spatial_axis="x")
        ),
    }


def _coordinate_dataset(coordinates: np.ndarray, *, spatial_axis: str):
    if coordinates.ndim == 1:
        return metric_value(
            coordinates,
            unit="um",
            dim_desc=(spatial_axis,),
        )
    if coordinates.ndim != 3:
        raise ValueError(
            "profile coordinates must have shape (sample,) or "
            "(radius, branch, sample)."
        )
    return metric_value(
        np.transpose(coordinates, (2, 1, 0)),
        unit="um",
        dim_desc=(spatial_axis, "branch", "radius"),
    )


def _profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    spatial_axis: str = "x",
    unit: str = "mm/s",
) -> DatasetValue:
    if profiles.ndim != 4:
        raise ValueError(
            "profile arrays must have shape "
            "(radius, branch, frame, spatial_sample)."
        )
    profiles_per_beat = interpolate_velocity_profiles_per_beat(
        profiles,
        cycle_boundary_indexes,
        index_base=index_base,
    )
    return DatasetValue(
        data=profiles_per_beat,
        attrs={
            "unit": unit,
            "dimDesc": [spatial_axis, "time", "beat", "branch", "radius"],
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
