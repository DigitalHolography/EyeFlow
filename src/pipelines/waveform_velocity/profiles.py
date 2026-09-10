"""Minimal HDF5 export for transverse and longitudinal velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from calculations.math import nanmean_float32
from input_output.schema import EyeFlowOutputPaths, VelocityProfileOutputPaths
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
        schema.artery_velocity_profiles,
        artery_segments,
        cycle_boundary_indexes,
        index_base=index_base,
        include_temporal_means=True,
    )
    metrics.update(
        _pack_vessel_profiles(
            schema.vein_velocity_profiles,
            vein_segments,
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return metrics


def _pack_vessel_profiles(
    paths: VelocityProfileOutputPaths,
    segments,
    cycle_boundary_indexes,
    *,
    index_base: int,
    include_temporal_means: bool = False,
) -> dict[str, object]:
    transverse_unmasked = _profile_dataset(
        np.asarray(segments.velocity_profiles, dtype=np.float32),
        cycle_boundary_indexes,
        index_base=index_base,
    )
    transverse_masked = _profile_dataset(
        np.asarray(
            segments.transverse_velocity_profiles_masked,
            dtype=np.float32,
        ),
        cycle_boundary_indexes,
        index_base=index_base,
        spatial_axis="x",
    )
    longitudinal_unmasked = _profile_dataset(
        np.asarray(
            segments.longitudinal_velocity_profiles_unmasked,
            dtype=np.float32,
        ),
        cycle_boundary_indexes,
        index_base=index_base,
        spatial_axis="y",
    )
    longitudinal_masked = _profile_dataset(
        np.asarray(
            segments.longitudinal_velocity_profiles_masked,
            dtype=np.float32,
        ),
        cycle_boundary_indexes,
        index_base=index_base,
        spatial_axis="y",
    )
    outputs = {
        paths.transverse_velocity_profile_unmasked: transverse_unmasked,
        paths.transverse_velocity_profile_masked: transverse_masked,
        paths.longitudinal_velocity_profile_unmasked: longitudinal_unmasked,
        paths.longitudinal_velocity_profile_masked: longitudinal_masked,
    }
    if include_temporal_means:
        outputs.update(
            {
                paths.transverse_velocity_profile_unmasked_meaned: (
                    _temporally_meaned_profile_dataset(transverse_unmasked)
                ),
                paths.transverse_velocity_profile_masked_meaned: (
                    _temporally_meaned_profile_dataset(transverse_masked)
                ),
                paths.longitudinal_velocity_profile_unmasked_meaned: (
                    _temporally_meaned_profile_dataset(longitudinal_unmasked)
                ),
                paths.longitudinal_velocity_profile_masked_meaned: (
                    _temporally_meaned_profile_dataset(longitudinal_masked)
                ),
            }
        )
    return outputs


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


def _temporally_meaned_profile_dataset(profile: DatasetValue) -> DatasetValue:
    """Average an interpolated profile over time within each beat."""
    data = nanmean_float32(np.asarray(profile.data), axis=1)
    attrs = dict(profile.attrs or {})
    dim_desc = list(attrs.get("dimDesc", ()))
    if len(dim_desc) < 2 or dim_desc[1] != "time":
        raise ValueError("profile dataset must have time as its second dimension.")
    del dim_desc[1]
    attrs["dimDesc"] = dim_desc
    attrs["temporal_reduction"] = "mean_over_interpolated_beat_time"
    return DatasetValue(
        data=data,
        attrs=attrs,
        h5_options=_profile_h5_options(data.shape),
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
