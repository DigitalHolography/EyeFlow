"""Cross-section profiles extracted from the quantitative moment0ff gradient."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity import segment_velocity_results
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
    meaned_profile: DatasetValue,
    *,
    minimum_gap: int,
) -> dict[str, DatasetValue]:
    """Find the two highest samples separated by ``minimum_gap`` samples."""

    values = np.asarray(meaned_profile.data, dtype=np.float32)
    attrs = dict(meaned_profile.attrs or {})
    dim_desc = list(attrs.get("dimDesc", ()))
    if values.ndim != 4 or dim_desc != ["x", "beat", "branch", "radius"]:
        raise ValueError(
            "meaned spatial-gradient profiles must have dimensions "
            "(x, beat, branch, radius)."
        )
    if minimum_gap < 1:
        raise ValueError("minimum peak gap must be at least one sample.")

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

    common_attrs = {
        "dimDesc": ["beat", "branch", "radius"],
        "minimum_peak_gap_samples": np.int32(minimum_gap),
        "peak_selection": "two_highest_finite_values_with_minimum_index_gap",
        "source_profile": "TransverseSpatialGradientProfileMaskedMeaned",
    }
    index_attrs = {
        **common_attrs,
        "unit": "pixels",
        "index_base": np.int32(0),
        "position_axis": "x",
    }
    value_attrs = {
        **common_attrs,
        "unit": attrs.get("unit", "a.u."),
    }
    lumen_diameter = right_indexes - left_indexes
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
        "lumen_diameter": DatasetValue(
            lumen_diameter,
            {
                **common_attrs,
                "unit": "pixels",
                "definition": "index_right_max - index_left_max",
            },
        ),
    }


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
