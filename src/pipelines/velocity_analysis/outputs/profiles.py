"""HDF5 export for transverse and longitudinal velocity profiles."""

from __future__ import annotations

import numpy as np

from input_output.profile_datasets import (
    _profile_dataset,
    _profile_h5_options,
    _temporally_meaned_profile_dataset,
)
from input_output.schema import EyeFlowOutputPaths, VelocityProfileOutputPaths
from pipeline_engine.base import DatasetValue
from pipelines.velocity.models import RetinalVelocity
from pipelines.velocity.semantics import velocity_dataset_attrs

from ..analysis.profiles.profiles import DEFAULT_PROFILE_MASK_DILATION_ITERATIONS
from .paths import resolve_output_paths


def pack_cross_section_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    index_base: int = 0,
    velocity: RetinalVelocity | None = None,
) -> dict[str, object]:
    schema = resolve_output_paths(output_paths)
    velocity_attrs = velocity_dataset_attrs(velocity)
    metrics = _pack_vessel_profiles(
        schema.artery_velocity_profiles,
        artery_segments,
        cycle_boundary_indexes,
        index_base=index_base,
        include_temporal_means=True,
        velocity_attrs=velocity_attrs,
    )
    metrics.update(
        _pack_vessel_profiles(
            schema.vein_velocity_profiles,
            vein_segments,
            cycle_boundary_indexes,
            index_base=index_base,
            velocity_attrs=velocity_attrs,
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
    velocity_attrs: dict[str, object] | None = None,
) -> dict[str, object]:
    profile = segments.profile
    valid_segments = np.asarray(profile.topology.valid_segments, dtype=bool)
    transverse_masked = np.asarray(
        profile.transverse.masked,
        dtype=np.float32,
    )
    outputs = {
        paths.transverse_velocity_profile_unmasked: _profile_dataset(
            np.asarray(profile.transverse.unmasked, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
            unit=str((velocity_attrs or {}).get("unit", "mm/s")),
            attrs=velocity_attrs,
            valid_segments=valid_segments,
        ),
        paths.transverse_velocity_profile_masked: _profile_dataset(
            transverse_masked,
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
            unit=str((velocity_attrs or {}).get("unit", "mm/s")),
            attrs=velocity_attrs,
            valid_segments=valid_segments,
        ),
        paths.longitudinal_velocity_profile_unmasked: _profile_dataset(
            np.asarray(
                profile.longitudinal.unmasked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
            unit=str((velocity_attrs or {}).get("unit", "mm/s")),
            attrs=velocity_attrs,
            valid_segments=valid_segments,
        ),
        paths.longitudinal_velocity_profile_masked: _profile_dataset(
            np.asarray(
                profile.longitudinal.masked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
            unit=str((velocity_attrs or {}).get("unit", "mm/s")),
            attrs=velocity_attrs,
            valid_segments=valid_segments,
        ),
    }
    if include_temporal_means:
        for field in (
            "transverse_velocity_profile_unmasked",
            "transverse_velocity_profile_masked",
            "longitudinal_velocity_profile_unmasked",
            "longitudinal_velocity_profile_masked",
        ):
            outputs[getattr(paths, field + "_meaned")] = (
                _temporally_meaned_profile_dataset(outputs[getattr(paths, field)])
            )
    return outputs


def pack_velocity_profile_fft_outputs(
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack FFT profiles accumulated during streamed segment processing."""

    schema = resolve_output_paths(output_paths)
    outputs = _pack_vessel_velocity_fft_profiles(
        schema.artery_velocity_profiles,
        artery_segments,
        mask_dilation_pixels=DEFAULT_PROFILE_MASK_DILATION_ITERATIONS,
    )
    outputs.update(
        _pack_vessel_velocity_fft_profiles(
            schema.vein_velocity_profiles,
            vein_segments,
            mask_dilation_pixels=0,
        )
    )
    return outputs



def _pack_vessel_velocity_fft_profiles(
    paths: VelocityProfileOutputPaths,
    segments,
    *,
    mask_dilation_pixels: int,
) -> dict[str, object]:
    unmasked_path = paths.transverse_velocity_profile_fft_unmasked
    masked_path = paths.transverse_velocity_profile_fft_masked
    if segments is None or unmasked_path is None or masked_path is None:
        return {}
    fft = segments.transverse_fft
    if fft is None:
        raise RuntimeError(
            "Velocity FFT profiles were not accumulated during streamed "
            "segment processing."
        )
    unmasked = np.asarray(fft.unmasked, dtype=np.float32)
    masked = np.asarray(fft.masked, dtype=np.float32)
    if unmasked.ndim != 5 or masked.shape != unmasked.shape:
        raise ValueError(
            "Streamed FFT profiles must have matching "
            "(x, frequency, beat, branch, radius) shapes."
        )
    shared_attrs = {
        "unit": "a.u.",
        "dimDesc": ["x", "frequency", "beat", "branch", "radius"],
        "source_dimensions": [
            "x",
            "y",
            "time",
            "beat",
            "branch",
            "radius",
        ],
        "temporal_transform": "absolute_value_of_fft_over_time",
        "fft_spectrum": "full",
        "fft_normalization": "none",
        "spatial_reduction": "nanmean_over_y",
    }
    return {
        unmasked_path: DatasetValue(
            data=unmasked,
            attrs={
                **shared_attrs,
                "mask_applied": False,
                "mask_dilation_iterations": 0,
            },
            h5_options=_profile_h5_options(unmasked.shape),
        ),
        masked_path: DatasetValue(
            data=masked,
            attrs={
                **shared_attrs,
                "mask_applied": True,
                "mask_dilation_iterations": int(mask_dilation_pixels),
                "mask_dilation_axis": "transverse_x",
            },
            h5_options=_profile_h5_options(masked.shape),
        ),
    }
