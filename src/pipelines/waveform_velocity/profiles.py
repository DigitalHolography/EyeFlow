"""HDF5 export for transverse and longitudinal velocity profiles."""

from __future__ import annotations

import numpy as np

from calculations.math import nanmean_float32
from calculations.topology import dilate_segment_masks
from input_output.profile_datasets import (
    _profile_dataset,
    _profile_h5_options,
    _temporally_meaned_profile_dataset,
)
from input_output.schema import EyeFlowOutputPaths, VelocityProfileOutputPaths
from pipeline_engine.base import DatasetValue
from pipelines.retinal_velocity.models import RetinalVelocity
from pipelines.retinal_velocity.semantics import (
    resolve_velocity_semantics,
)

_PROFILE_MASK_DILATION_ITERATIONS = 10

def pack_cross_section_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    index_base: int = 0,
    velocity_analysis: RetinalVelocity | None = None,
) -> dict[str, object]:
    schema = _resolve_output_paths(output_paths)
    velocity_unit = resolve_velocity_semantics(velocity_analysis).unit
    metrics = _pack_vessel_profiles(
        schema.artery_velocity_profiles,
        artery_segments,
        cycle_boundary_indexes,
        index_base=index_base,
        include_temporal_means=True,
        velocity_unit=velocity_unit,
    )
    metrics.update(
        _pack_vessel_profiles(
            schema.vein_velocity_profiles,
            vein_segments,
            cycle_boundary_indexes,
            index_base=index_base,
            velocity_unit=velocity_unit,
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
    velocity_unit: str = "mm/s",
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
            unit=velocity_unit,
            valid_segments=valid_segments,
        ),
        paths.transverse_velocity_profile_masked: _profile_dataset(
            transverse_masked,
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
            unit=velocity_unit,
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
            unit=velocity_unit,
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
            unit=velocity_unit,
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


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)

def pack_velocity_profile_fft_outputs(
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack FFT profiles accumulated during streamed segment processing."""

    schema = _resolve_output_paths(output_paths)
    outputs = _pack_vessel_velocity_fft_profiles(
        schema.artery_velocity_profiles,
        artery_segments,
        mask_dilation_pixels=_PROFILE_MASK_DILATION_ITERATIONS,
    )
    outputs.update(
        _pack_vessel_velocity_fft_profiles(
            schema.vein_velocity_profiles,
            vein_segments,
            mask_dilation_pixels=0,
        )
    )
    return outputs



def velocity_fft_transverse_profiles(
    velocity_maps_per_beat: np.ndarray,
    segment_masks: np.ndarray,
    *,
    mask_dilation_pixels: int = _PROFILE_MASK_DILATION_ITERATIONS,
) -> tuple[np.ndarray, np.ndarray]:
    """Return FFT profiles shaped ``(x, frequency, beat, branch, radius)``.

    ``velocity_maps_per_beat`` must already have shape
    ``(x, y, time, beat, branch, radius)``. The FFT is applied along its time
    axis independently for every pixel, beat, branch, and radius.

    The two returned arrays contain the unmasked and horizontally expanded
    mask projections used by the historical artery profile workflow.
    """

    maps = np.asarray(velocity_maps_per_beat, dtype=np.float32)
    masks = np.asarray(segment_masks, dtype=bool)
    if maps.ndim != 6:
        raise ValueError(
            "velocity_maps_per_beat must have shape "
            "(x, y, time, beat, branch, radius)."
        )
    expected_mask_shape = (
        maps.shape[5],
        maps.shape[4],
        maps.shape[1],
        maps.shape[0],
    )
    if masks.shape != expected_mask_shape:
        raise ValueError(
            "segment_masks must have shape (radius, branch, y, x) matching "
            "velocity_maps_per_beat."
        )

    dilated_masks = dilate_segment_masks(
        masks,
        iterations=mask_dilation_pixels,
        horizontal_only=True,
    )
    output_shape = (
        maps.shape[0],
        maps.shape[2],
        maps.shape[3],
        maps.shape[4],
        maps.shape[5],
    )
    unmasked = np.full(output_shape, np.nan, dtype=np.float32)
    masked = np.full(output_shape, np.nan, dtype=np.float32)
    for radius_index in range(maps.shape[5]):
        for branch_index in range(maps.shape[4]):
            magnitude = np.abs(
                np.fft.fft(
                    maps[..., branch_index, radius_index],
                    axis=2,
                )
            ).astype(np.float32, copy=False)
            unmasked[..., branch_index, radius_index] = nanmean_float32(
                magnitude,
                axis=1,
            )
            xy_mask = dilated_masks[radius_index, branch_index].T
            masked[..., branch_index, radius_index] = nanmean_float32(
                np.where(
                    xy_mask[:, :, None, None],
                    magnitude,
                    np.float32(np.nan),
                ),
                axis=1,
            )
    return unmasked, masked



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
