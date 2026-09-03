"""HDF5 export for transverse and longitudinal velocity/displacement profiles."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.cross_section.profile_processing import (
    interpolate_velocity_profiles_per_beat,
)
from input_output.schema import EyeFlowOutputPaths, VelocityProfileOutputPaths
from pipeline_engine.base import DatasetValue


_DISPLACEMENT_PROFILE_ROOT = "Processing/DisplacementProfiles"
_DISPLACEMENT_PROFILE_FIELDS = (
    ("X", "x_sum_displacement_profile", "local_x"),
    ("Y", "y_sum_displacement_profile", "local_y"),
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
    """Pack signed local-component sums for each vessel segment."""

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
        for axis_name, profile_field, component in _DISPLACEMENT_PROFILE_FIELDS:
            outputs[f"{root}/{axis_name}_sum_displacement_profile/value"] = (
                _summed_displacement_profile_dataset(
                    np.asarray(
                        getattr(displacement, profile_field),
                        dtype=np.float32,
                    ),
                    cycle_boundary_indexes,
                    index_base=index_base,
                    component=component,
                )
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


def _summed_displacement_profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    component: str,
) -> DatasetValue:
    """Interpolate one signed unmasked-subimage component sum per beat."""

    if profiles.ndim != 3:
        raise ValueError(
            "summed displacement profiles must have shape "
            "(radius, branch, frame)."
        )
    profiles_per_beat = interpolate_velocity_profiles_per_beat(
        profiles[..., None],
        cycle_boundary_indexes,
        index_base=index_base,
    )
    return DatasetValue(
        data=profiles_per_beat[0],
        attrs={
            "unit": "pixels",
            "dimDesc": ["time", "beat", "branch", "radius"],
            "coordinate_system": "rotated_segment_pixel",
            "component_basis": "rotated_segment_local",
            "component": component,
            "spatial_region": "full_unmasked_subimage",
            "spatial_reduction": "sum_over_valid_subimage_pixels",
            "displacement_reference": "temporal_mean_image",
        },
        h5_options=_profile_h5_options(profiles_per_beat.shape[1:]),
    )


def _profile_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    """Use lossless compression with chunks aligned to one segment profile."""
    options: dict[str, object] = {
        "compression": "gzip",
        "compression_opts": 4,
        "shuffle": True,
    }
    if len(shape) not in (4, 5) or not all(shape):
        return options

    sample_count, time_count = shape[:2]
    target_elements = (1024 * 1024) // np.dtype(np.float32).itemsize
    if len(shape) == 4:
        options["chunks"] = (sample_count, 1, 1, 1)
        return options
    time_chunk = min(time_count, max(target_elements // sample_count, 1))
    options["chunks"] = (sample_count, time_chunk, 1, 1, 1)
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
