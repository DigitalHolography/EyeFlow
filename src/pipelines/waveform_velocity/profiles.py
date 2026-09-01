"""Minimal HDF5 export for transverse and longitudinal velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from input_output.schema import EyeFlowOutputPaths, VelocityProfileOutputPaths
from pipeline_engine.base import DatasetValue


_DISPLACEMENT_PROFILE_ROOT = "Processing/DisplacementProfiles"
_DEBUG_DISPLACEMENT_PROFILE_ROOT = "Processing/Debug/DisplacementProfiles"
_DISPLACEMENT_PROFILE_FIELDS = (
    (
        "TransverseDisplacementProfileUnmasked",
        "displacement_profiles",
        "displacement_vector_profiles",
        "x",
    ),
    (
        "TransverseDisplacementProfileMasked",
        "transverse_displacement_profiles_masked",
        "transverse_displacement_vector_profiles_masked",
        "x",
    ),
    (
        "LongitudinalDisplacementProfileUnmasked",
        "longitudinal_displacement_profiles_unmasked",
        "longitudinal_displacement_vector_profiles_unmasked",
        "y",
    ),
    (
        "LongitudinalDisplacementProfileMasked",
        "longitudinal_displacement_profiles_masked",
        "longitudinal_displacement_vector_profiles_masked",
        "y",
    ),
)


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


def pack_displacement_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    outputs = _pack_vessel_displacement_profiles(
        artery_segments,
        "Artery",
        cycle_boundary_indexes,
        index_base=index_base,
    )
    outputs.update(
        _pack_vessel_displacement_profiles(
            vein_segments,
            "Vein",
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return outputs


def _pack_vessel_displacement_profiles(
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
        root = f"{_DISPLACEMENT_PROFILE_ROOT}/{method}/{vessel_name}"
        debug_root = (
            f"{_DEBUG_DISPLACEMENT_PROFILE_ROOT}/{method}/{vessel_name}"
        )
        for leaf, scalar_field, vector_field, spatial_axis in (
            _DISPLACEMENT_PROFILE_FIELDS
        ):
            outputs[f"{root}/{leaf}/value"] = _profile_dataset(
                np.asarray(getattr(displacement, scalar_field), dtype=np.float32),
                cycle_boundary_indexes,
                index_base=index_base,
                spatial_axis=spatial_axis,
                unit="pixels",
            )
            outputs[f"{debug_root}/{leaf}/value"] = _vector_profile_dataset(
                np.asarray(getattr(displacement, vector_field), dtype=np.float32),
                cycle_boundary_indexes,
                index_base=index_base,
                spatial_axis=spatial_axis,
            )
    return outputs


def _pack_vessel_profiles(
    paths: VelocityProfileOutputPaths,
    segments,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> dict[str, object]:
    return {
        paths.transverse_velocity_profile_unmasked: _profile_dataset(
            np.asarray(segments.velocity_profiles, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
        ),
        paths.transverse_velocity_profile_masked: _profile_dataset(
            np.asarray(
                segments.transverse_velocity_profiles_masked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
        ),
        paths.longitudinal_velocity_profile_unmasked: _profile_dataset(
            np.asarray(
                segments.longitudinal_velocity_profiles_unmasked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
        ),
        paths.longitudinal_velocity_profile_masked: _profile_dataset(
            np.asarray(
                segments.longitudinal_velocity_profiles_masked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
        ),
    }


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


def _vector_profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    spatial_axis: str,
) -> DatasetValue:
    if profiles.ndim != 5 or profiles.shape[-1] != 2:
        raise ValueError(
            "displacement vector profiles must have shape "
            "(radius, branch, frame, spatial_sample, 2)."
        )
    profiles_per_beat = np.stack(
        [
            interpolate_velocity_profiles_per_beat(
                profiles[..., component_index],
                cycle_boundary_indexes,
                index_base=index_base,
            )
            for component_index in range(2)
        ],
        axis=-1,
    )
    return DatasetValue(
        data=profiles_per_beat,
        attrs={
            "unit": "pixels",
            "dimDesc": [
                spatial_axis,
                "time",
                "beat",
                "branch",
                "radius",
                "displacement_orientation",
            ],
            "components": ["local_x", "local_y"],
            "temporary_debug_output": True,
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
    if len(shape) not in (5, 6) or not all(shape):
        return options

    x_count, time_count = shape[:2]
    target_elements = (1024 * 1024) // np.dtype(np.float32).itemsize
    time_chunk = min(time_count, max(target_elements // x_count, 1))
    options["chunks"] = (
        (x_count, time_chunk, 1, 1, 1)
        if len(shape) == 5
        else (x_count, time_chunk, 1, 1, 1, shape[-1])
    )
    return options


def _hdf_method_name(value: object) -> str:
    method = str(value).strip()
    if not method or "/" in method:
        raise ValueError(
            "Displacement registration method names must be non-empty HDF5 path segments."
        )
    return method


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
