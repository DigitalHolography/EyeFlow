"""Cross-section profiles extracted from the quantitative moment0ff gradient."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity import segment_velocity_results
from calculations.math import nanmean_float32
from pipeline_engine.base import DatasetValue
from pipelines.spatial_gradient_moment0.runner import (
    STATE_KEY,
    SpatialGradientMoment0Artifacts,
)
from pipelines.waveform_velocity_core.runner import _segment_ring_settings

from .profiles import (
    _profile_dataset,
    _temporally_meaned_profile_dataset,
)


SPATIAL_GRADIENT_PROFILE_ROOT = "Processing/SpatialGradientProfiles"
SPATIAL_GRADIENT_METRICS_ROOT = "Processing/SpatialGradientMetrics"
SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES = 5


def extract_spatial_gradient_segments(ctx, waveform_context):
    """Prepare annular branch topology and project the gradient video onto it."""

    artifacts = ctx.state.get(STATE_KEY)
    if not isinstance(artifacts, SpatialGradientMoment0Artifacts):
        raise RuntimeError(
            "The spatial_gradient_moment0 prerequisite did not publish its "
            "quantitative gradient video."
        )
    gradient_map = np.load(artifacts.gradient_path, mmap_mode="r")
    source = waveform_context.source_data
    try:
        ring_settings = _segment_ring_settings(
            source.optic_disc_width,
            source.optic_disc_height,
            image_shape=gradient_map.shape[-2:],
            optic_disc_center=source.optic_disc_center,
            number_of_radii_in_FOV=int(
                waveform_context.attrs["number_of_radii_in_FOV"]
            ),
        )
        return segment_velocity_results(
            gradient_map,
            source.retinal_artery_mask,
            source.retinal_vein_mask,
            source.optic_disc_center,
            ring_settings,
            source.cross_section_settings,
            retain_displacement_maps=False,
        )
    finally:
        mmap = getattr(gradient_map, "_mmap", None)
        if mmap is not None:
            mmap.close()
        artifacts.cleanup()


def pack_spatial_gradient_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    """Pack masked, unmasked, and temporal-mean transverse gradient profiles."""

    outputs = _pack_vessel_spatial_gradient_profiles(
        artery_segments,
        "Artery",
        cycle_boundary_indexes,
        index_base=index_base,
    )
    outputs.update(
        _pack_vessel_spatial_gradient_profiles(
            vein_segments,
            "Vein",
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return outputs


def cleanup_spatial_gradient_artifacts(ctx) -> None:
    artifacts = ctx.state.get(STATE_KEY)
    if isinstance(artifacts, SpatialGradientMoment0Artifacts):
        artifacts.cleanup()


def _pack_vessel_spatial_gradient_profiles(
    segments,
    vessel_name: str,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> dict[str, object]:
    if segments is None:
        return {}
    root = f"{SPATIAL_GRADIENT_PROFILE_ROOT}/{vessel_name}/Transverse"
    unmasked = _gradient_profile_dataset(
        np.asarray(segments.velocity_profiles, dtype=np.float32),
        cycle_boundary_indexes,
        index_base=index_base,
        mask="unmasked",
    )
    masked = _gradient_profile_dataset(
        np.asarray(
            segments.transverse_velocity_profiles_masked,
            dtype=np.float32,
        ),
        cycle_boundary_indexes,
        index_base=index_base,
        mask="vessel_segment",
    )
    meaned = _temporally_meaned_profile_dataset(masked)
    peak_metrics = _spatial_gradient_peak_metrics(
        masked,
        meaned,
        minimum_gap=SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES,
    )
    metrics_root = f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel_name}/Transverse"
    return {
        f"{root}/TransverseSpatialGradientProfileMasked": masked,
        f"{root}/TransverseSpatialGradientProfileMaskedMeaned": meaned,
        f"{root}/TransverseSpatialGradientProfileUnmasked": unmasked,
        **{
            f"{metrics_root}/{name}": value
            for name, value in peak_metrics.items()
        },
    }


def _spatial_gradient_peak_metrics(
    masked_profile: DatasetValue,
    meaned_profile: DatasetValue,
    *,
    minimum_gap: int,
) -> dict[str, DatasetValue]:
    """Find time-resolved peak positions and derive lumen diameters."""

    if minimum_gap < 1:
        raise ValueError("minimum peak gap must be at least one sample.")

    left_indexes, right_indexes, _, _ = _spatial_gradient_peak_arrays(
        masked_profile,
        ["x", "time", "beat", "branch", "radius"],
        minimum_gap=minimum_gap,
    )
    _, _, left_values, right_values = _spatial_gradient_peak_arrays(
        meaned_profile,
        ["x", "beat", "branch", "radius"],
        minimum_gap=minimum_gap,
    )
    meaned_attrs = dict(meaned_profile.attrs or {})
    common_attrs = {
        "minimum_peak_gap_samples": np.int32(minimum_gap),
        "peak_selection": "two_highest_finite_values_with_minimum_index_gap",
    }
    index_attrs = {
        **common_attrs,
        "dimDesc": ["time", "beat", "branch", "radius"],
        "source_profile": "TransverseSpatialGradientProfileMasked",
        "unit": "pixels",
        "index_base": np.int32(0),
        "position_axis": "x",
    }
    value_attrs = {
        **common_attrs,
        "dimDesc": ["beat", "branch", "radius"],
        "source_profile": "TransverseSpatialGradientProfileMaskedMeaned",
        "unit": meaned_attrs.get("unit", "a.u."),
    }
    time_mean_left = nanmean_float32(left_indexes, axis=0)
    time_mean_right = nanmean_float32(right_indexes, axis=0)
    radius_mean_left = nanmean_float32(left_indexes, axis=-1)
    radius_mean_right = nanmean_float32(right_indexes, axis=-1)
    diameter_attrs = {
        **common_attrs,
        "source_profile": "TransverseSpatialGradientProfileMasked",
        "unit": "pixels",
        "definition": "index_right_max - index_left_max",
    }
    temporal_reduction = "mean_over_interpolated_beat_time"
    radius_reduction = "mean_over_valid_radii"
    return {
        "index_left_max": DatasetValue(
            left_indexes,
            {**index_attrs, "peak_side": "left"},
        ),
        "index_right_max": DatasetValue(
            right_indexes,
            {**index_attrs, "peak_side": "right"},
        ),
        "peak_value_left_max": DatasetValue(
            left_values,
            {**value_attrs, "peak_side": "left"},
        ),
        "peak_value_right_max": DatasetValue(
            right_values,
            {**value_attrs, "peak_side": "right"},
        ),
        "lumen_diameter_branch": DatasetValue(
            nanmean_float32(time_mean_right, axis=-1)
            - nanmean_float32(time_mean_left, axis=-1),
            {
                **diameter_attrs,
                "dimDesc": ["beat", "branch"],
                "temporal_reduction": temporal_reduction,
                "radius_reduction": radius_reduction,
            },
        ),
        "lumen_diameter_radius": DatasetValue(
            time_mean_right - time_mean_left,
            {
                **diameter_attrs,
                "dimDesc": ["beat", "branch", "radius"],
                "temporal_reduction": temporal_reduction,
            },
        ),
        "lumen_diameter_time": DatasetValue(
            radius_mean_right - radius_mean_left,
            {
                **diameter_attrs,
                "dimDesc": ["time", "beat", "branch"],
                "radius_reduction": radius_reduction,
            },
        ),
        "lumen_diameter": DatasetValue(
            right_indexes - left_indexes,
            {**diameter_attrs, "dimDesc": ["time", "beat", "branch", "radius"]},
        ),
    }


def _spatial_gradient_peak_arrays(
    profile: DatasetValue,
    expected_dimensions: list[str],
    *,
    minimum_gap: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Find the two strongest separated peaks along each profile's x axis."""

    values = np.asarray(profile.data, dtype=np.float32)
    dim_desc = list((profile.attrs or {}).get("dimDesc", ()))
    if values.ndim != len(expected_dimensions) or dim_desc != expected_dimensions:
        raise ValueError(
            "spatial-gradient profile dimensions must be "
            f"{tuple(expected_dimensions)}."
        )
    metric_shape = values.shape[1:]
    left_indexes = np.full(metric_shape, np.nan, dtype=np.float32)
    right_indexes = np.full(metric_shape, np.nan, dtype=np.float32)
    left_values = np.full(metric_shape, np.nan, dtype=np.float32)
    right_values = np.full(metric_shape, np.nan, dtype=np.float32)
    for indexes in np.ndindex(metric_shape):
        curve = values[(slice(None), *indexes)]
        peaks = _two_highest_separated_indexes(curve, minimum_gap)
        if len(peaks) >= 1:
            left_indexes[indexes] = peaks[0]
            left_values[indexes] = curve[peaks[0]]
        if len(peaks) == 2:
            right_indexes[indexes] = peaks[1]
            right_values[indexes] = curve[peaks[1]]

    return left_indexes, right_indexes, left_values, right_values


def _two_highest_separated_indexes(
    profile: np.ndarray,
    minimum_gap: int,
) -> tuple[int, ...]:
    finite_indexes = np.flatnonzero(np.isfinite(profile))
    if finite_indexes.size == 0:
        return ()
    ranked = finite_indexes[
        np.lexsort((finite_indexes, -profile[finite_indexes]))
    ]
    strongest = int(ranked[0])
    second = next(
        (
            int(index)
            for index in ranked[1:]
            if abs(int(index) - strongest) >= minimum_gap
        ),
        None,
    )
    if second is None:
        return (strongest,)
    return tuple(sorted((strongest, second)))


def _gradient_profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    mask: str,
) -> DatasetValue:
    value = _profile_dataset(
        profiles,
        cycle_boundary_indexes,
        index_base=index_base,
        unit="a.u.",
    )
    attrs = dict(value.attrs or {})
    attrs.update(
        {
            "measurement": "spatial_gradient_magnitude",
            "source_dataset": "/moment0ff",
            "spatial_operator": "3x3 Sobel magnitude",
            "spatial_region": mask,
        }
    )
    value.attrs = attrs
    return value


__all__ = [
    "SPATIAL_GRADIENT_PROFILE_ROOT",
    "SPATIAL_GRADIENT_METRICS_ROOT",
    "SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES",
    "cleanup_spatial_gradient_artifacts",
    "extract_spatial_gradient_segments",
    "pack_spatial_gradient_profile_outputs",
]
