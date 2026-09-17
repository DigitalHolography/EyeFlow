"""Cross-section gradients computed after interpolating raw moment0ff."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from calculations.blood_flow_velocity import segment_velocity_results
from calculations.math import nanmedian
from pipeline_engine.base import DatasetValue
from pipelines.spatial_gradient_moment0.runner import (
    STATE_KEY,
    TBKR_LUMEN_SIZE_QC_THRESHOLD,
    TEMPORAL_MEDIAN_WINDOW,
    SpatialGradientMoment0Artifacts,
)
from pipelines.waveform_velocity_core.runner import _segment_ring_settings
from pipelines.waveform_velocity_core.segmentation import _optic_disc_mask

from .profiles import (
    _profile_dataset,
    _temporally_meaned_profile_dataset,
)

SPATIAL_GRADIENT_PROFILE_ROOT = "Processing/SpatialGradientProfiles"
SPATIAL_GRADIENT_METRICS_ROOT = "Processing/SpatialGradientMetrics"
SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES = 5
_SPATIAL_GRADIENT_MASK_DILATION_PIXELS = 5


def extract_spatial_gradient_segments(ctx, waveform_context):
    """Interpolate raw moment0ff, filter its gradient, then rotate each segment."""

    artifacts = ctx.state.get(STATE_KEY)
    if not isinstance(artifacts, SpatialGradientMoment0Artifacts):
        raise RuntimeError(
            "The spatial_gradient_moment0 prerequisite did not publish its "
            "quantitative gradient video."
        )
    source = waveform_context.source_data
    try:
        moment0ff = ctx.inputs.hd.as_holodoppler().moment0_flat_field_dataset()
        if moment0ff is None:
            raise KeyError("Missing flat-field HoloDoppler moment0 dataset: moment0ff/M0FF.")
        ring_settings = _segment_ring_settings(
            source.optic_disc_width,
            source.optic_disc_height,
            image_shape=moment0ff.shape[-2:],
            optic_disc_center=source.optic_disc_center,
            number_of_radii_in_FOV=int(
                waveform_context.attrs["number_of_radii_in_FOV"]
            ),
        )
        optic_disc_mask, _ = _optic_disc_mask(
            source, source.retinal_artery_mask.shape
        )
        return segment_velocity_results(
            moment0ff,
            source.retinal_artery_mask,
            source.retinal_vein_mask,
            source.optic_disc_center,
            ring_settings,
            replace(source.cross_section_settings, spatial_gradient=True),
            optic_disc_mask=optic_disc_mask if np.any(optic_disc_mask) else None,
            artery_transverse_mask_dilation_pixels=(
                _SPATIAL_GRADIENT_MASK_DILATION_PIXELS
            ),
            vein_transverse_mask_dilation_pixels=(
                _SPATIAL_GRADIENT_MASK_DILATION_PIXELS
            ),
            retain_displacement_maps=False,
        )
    finally:
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
        unmasked_profile=unmasked,
        minimum_gap=SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES,
    )
    metrics_root = f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel_name}/Transverse"
    return {
        f"{root}/Masked/SpatialGradientProfile/value": masked,
        f"{root}/Masked/SpatialGradientProfileMeaned/value": meaned,
        f"{root}/Unmasked/SpatialGradientProfile/value": unmasked,
        **{
            f"{metrics_root}/{name}": value
            for name, value in peak_metrics.items()
        },
    }


def _spatial_gradient_peak_metrics(
    masked_profile: DatasetValue,
    meaned_profile: DatasetValue,
    *,
    unmasked_profile: DatasetValue | None = None,
    minimum_gap: int,
) -> dict[str, DatasetValue]:
    """Find time-resolved peak positions and derive lumen diameters."""

    if minimum_gap < 1:
        raise ValueError("minimum peak gap must be at least one sample.")

    masked_left_indexes, masked_right_indexes, _, _ = (
        _spatial_gradient_peak_arrays(
            masked_profile,
            ["x", "time", "beat", "branch", "radius"],
            minimum_gap=minimum_gap,
        )
    )
    common_attrs = {
        "minimum_peak_gap_samples": np.int32(minimum_gap),
        "peak_selection": "two_highest_finite_values_with_minimum_index_gap",
    }
    if unmasked_profile is not None:
        metrics = _spatial_gradient_edge_index_metrics(
            masked_left_indexes,
            masked_right_indexes,
            mask_name="Masked",
            source_profile="Masked/SpatialGradientProfile/value",
            common_attrs=common_attrs,
        )
    else:
        index_attrs = {
            **common_attrs,
            "dimDesc": ["time", "beat", "branch", "radius"],
            "source_profile": "Masked/SpatialGradientProfile/value",
            "unit": "pixels",
            "index_base": np.int32(0),
            "position_axis": "x",
            "index_refinement": "three_point_quadratic_vertex_with_integer_fallback",
        }
        metrics = {
            "index_left_max": DatasetValue(
                masked_left_indexes,
                {**index_attrs, "peak_side": "left"},
            ),
            "index_right_max": DatasetValue(
                masked_right_indexes,
                {**index_attrs, "peak_side": "right"},
            ),
        }
    if unmasked_profile is None:
        _, _, left_values, right_values = _spatial_gradient_peak_arrays(
            meaned_profile,
            ["x", "beat", "branch", "radius"],
            minimum_gap=minimum_gap,
        )
        meaned_attrs = dict(meaned_profile.attrs or {})
        value_attrs = {
            **common_attrs,
            "dimDesc": ["beat", "branch", "radius"],
            "source_profile": "Masked/SpatialGradientProfileMeaned/value",
            "unit": meaned_attrs.get("unit", "a.u."),
        }
        metrics.update(
            {
                "peak_value_left_max": DatasetValue(
                    left_values,
                    {**value_attrs, "peak_side": "left"},
                ),
                "peak_value_right_max": DatasetValue(
                    right_values,
                    {**value_attrs, "peak_side": "right"},
                ),
            }
        )
    if unmasked_profile is not None:
        unmasked_left_indexes, unmasked_right_indexes, _, _ = (
            _spatial_gradient_peak_arrays(
                unmasked_profile,
                ["x", "time", "beat", "branch", "radius"],
                minimum_gap=minimum_gap,
            )
        )
        metrics.update(
            _spatial_gradient_edge_index_metrics(
                unmasked_left_indexes,
                unmasked_right_indexes,
                mask_name="Unmasked",
                source_profile="Unmasked/SpatialGradientProfile/value",
                common_attrs=common_attrs,
            )
        )
    return metrics


def _spatial_gradient_edge_index_metrics(
    left_indexes: np.ndarray,
    right_indexes: np.ndarray,
    *,
    mask_name: str,
    source_profile: str,
    common_attrs: dict[str, object],
) -> dict[str, DatasetValue]:
    """Build the requested hierarchy of NaN-median edge-index reductions."""

    left_bkr = nanmedian(left_indexes, axis=0)
    right_bkr = nanmedian(right_indexes, axis=0)
    left_kr = nanmedian(left_bkr, axis=0)
    right_kr = nanmedian(right_bkr, axis=0)
    hierarchy = {
        "tbkr": (
            left_indexes,
            right_indexes,
            ["time", "beat", "branch", "radius"],
        ),
        "bkr": (
            left_bkr,
            right_bkr,
            ["beat", "branch", "radius"],
        ),
        "bk": (
            nanmedian(left_bkr, axis=-1),
            nanmedian(right_bkr, axis=-1),
            ["beat", "branch"],
        ),
        "kr": (
            left_kr,
            right_kr,
            ["branch", "radius"],
        ),
        "k": (
            nanmedian(left_kr, axis=-1),
            nanmedian(right_kr, axis=-1),
            ["branch"],
        ),
    }
    tbkr_lumen_size = right_indexes - left_indexes
    tbkr_distribution_qc, tbkr_median, tbkr_standard_deviation = (
        _lumen_size_standard_deviation_quality_control(tbkr_lumen_size)
    )
    kr_lumen_size_qc = _kr_lumen_size_quality_control(tbkr_distribution_qc)
    tbkr_lumen_size_qc = _tbkr_lumen_size_quality_control(
        kr_lumen_size_qc,
        tbkr_lumen_size.shape,
        threshold=TBKR_LUMEN_SIZE_QC_THRESHOLD,
    )
    metrics: dict[str, DatasetValue] = {}
    for dimension_tag, (left, right, dim_desc) in hierarchy.items():
        attrs = {
            **common_attrs,
            "dimDesc": dim_desc,
            "source_profile": source_profile,
            "unit": "pixels",
            "index_base": np.int32(0),
            "position_axis": "x",
            "index_refinement": "three_point_quadratic_vertex_with_integer_fallback",
        }
        metrics[f"{mask_name}/{dimension_tag}/left_edge_index"] = DatasetValue(
            left,
            {**attrs, "peak_side": "left"},
        )
        metrics[f"{mask_name}/{dimension_tag}/right_edge_index"] = DatasetValue(
            right,
            {**attrs, "peak_side": "right"},
        )
        lumen_size = right - left
        if dimension_tag == "tbkr":
            lumen_size_qc = tbkr_lumen_size_qc
            median = tbkr_median
            standard_deviation = tbkr_standard_deviation
            qc_unit = "binary"
            qc_attrs = {
                "definition": (
                    "1 when the corresponding kr lumen_size QC is strictly "
                    "greater than the threshold; 0 otherwise. Time and beat "
                    "dimensions are broadcast and do not affect the result"
                ),
                "source_metric": f"{mask_name}/kr/lumen/size_qc",
                "threshold": np.float32(TBKR_LUMEN_SIZE_QC_THRESHOLD),
                "comparison": "greater_than",
                "broadcast_dimensions": ["time", "beat"],
            }
        elif dimension_tag == "kr":
            lumen_size_qc = kr_lumen_size_qc
            qc_unit = "fraction"
            qc_attrs = {
                "definition": (
                    "fraction of time-and-beat lumen_size samples within the "
                    "inclusive median plus or minus one population standard "
                    "deviation"
                ),
                "source_metric": f"{mask_name}/tbkr/lumen/size",
                "reduction": "mean_over_time_and_beat",
                "median": tbkr_median,
                "standard_deviation": tbkr_standard_deviation,
                "ddof": np.int32(0),
                "lower_limit": np.float32(
                    tbkr_median - tbkr_standard_deviation
                ),
                "upper_limit": np.float32(
                    tbkr_median + tbkr_standard_deviation
                ),
            }
        else:
            lumen_size_qc, lower_limit, upper_limit = (
                _lumen_size_quality_control(lumen_size)
            )
            qc_unit = "binary"
            qc_attrs = {
                "definition": (
                    "1 for finite lumen_size values within the inclusive "
                    "0.5th-to-99.5th percentile interval; 0 otherwise"
                ),
                "lower_percentile": np.float32(0.5),
                "upper_percentile": np.float32(99.5),
                "lower_limit": lower_limit,
                "upper_limit": upper_limit,
            }
        lumen_path = f"{mask_name}/{dimension_tag}/lumen/size"
        metrics[lumen_path] = DatasetValue(
            lumen_size,
            {
                **attrs,
                "definition": "right_edge_index - left_edge_index",
            },
        )
        metrics[f"{lumen_path}_qc"] = DatasetValue(
            lumen_size_qc,
            {
                **attrs,
                "unit": qc_unit,
                **qc_attrs,
            },
        )
        if mask_name == "Masked" and dimension_tag == "tbkr":
            statistic_attrs = {
                **attrs,
                "dimDesc": [],
                "source_metric": "Masked/tbkr/lumen/size",
                "distribution": "all_finite_values",
            }
            metrics[f"{lumen_path}_median"] = DatasetValue(
                median,
                {
                    **statistic_attrs,
                    "statistic": "median",
                },
            )
            metrics[f"{lumen_path}_std"] = DatasetValue(
                standard_deviation,
                {
                    **statistic_attrs,
                    "statistic": "standard_deviation",
                    "ddof": np.int32(0),
                },
            )
    return metrics


def _lumen_size_quality_control(
    lumen_size: np.ndarray,
) -> tuple[np.ndarray, np.float32, np.float32]:
    """Flag the inclusive central 99% of a finite lumen-size distribution."""

    values = np.asarray(lumen_size, dtype=np.float32)
    qc = np.zeros(values.shape, dtype=np.uint8)
    finite = np.isfinite(values)
    if not np.any(finite):
        nan = np.float32(np.nan)
        return qc, nan, nan

    lower_limit, upper_limit = np.percentile(values[finite], [0.5, 99.5])
    qc[finite & (values >= lower_limit) & (values <= upper_limit)] = 1
    return qc, np.float32(lower_limit), np.float32(upper_limit)


def _lumen_size_distribution_statistics(
    lumen_size: np.ndarray,
) -> tuple[np.float32, np.float32]:
    """Return the median and population standard deviation of finite values."""

    values = np.asarray(lumen_size, dtype=np.float32)
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        nan = np.float32(np.nan)
        return nan, nan
    return (
        np.float32(np.median(finite_values)),
        np.float32(np.std(finite_values, dtype=np.float64, ddof=0)),
    )


def _lumen_size_standard_deviation_quality_control(
    lumen_size: np.ndarray,
) -> tuple[np.ndarray, np.float32, np.float32]:
    """Flag finite values within the inclusive median-plus/minus-std interval."""

    values = np.asarray(lumen_size, dtype=np.float32)
    qc = np.zeros(values.shape, dtype=np.uint8)
    median, standard_deviation = _lumen_size_distribution_statistics(values)
    if not np.isfinite(median) or not np.isfinite(standard_deviation):
        return qc, median, standard_deviation

    finite = np.isfinite(values)
    lower_limit = median - standard_deviation
    upper_limit = median + standard_deviation
    qc[finite & (values >= lower_limit) & (values <= upper_limit)] = 1
    return qc, median, standard_deviation


def _kr_lumen_size_quality_control(
    tbkr_distribution_qc: np.ndarray,
) -> np.ndarray:
    """Average preliminary time/beat QC into a branch/radius fraction."""

    tbkr_qc = np.asarray(tbkr_distribution_qc, dtype=np.uint8)
    if tbkr_qc.ndim != 4:
        raise ValueError(
            "tbkr lumen-size distribution QC must have dimensions (t, b, k, r)."
        )
    if tbkr_qc.shape[0] == 0 or tbkr_qc.shape[1] == 0:
        return np.zeros(tbkr_qc.shape[2:], dtype=np.float32)

    kr_qc = np.mean(tbkr_qc, axis=(0, 1), dtype=np.float32)
    kr_qc = np.asarray(kr_qc, dtype=np.float32)
    return np.clip(kr_qc, 0.0, 1.0).astype(np.float32, copy=False)


def _tbkr_lumen_size_quality_control(
    kr_lumen_size_qc: np.ndarray,
    tbkr_shape: tuple[int, ...],
    *,
    threshold: float,
) -> np.ndarray:
    """Threshold kr QC and broadcast it across unchanged time/beat axes."""

    kr_qc = np.asarray(kr_lumen_size_qc, dtype=np.float32)
    if kr_qc.ndim != 2:
        raise ValueError("kr lumen-size QC must have dimensions (k, r).")
    if len(tbkr_shape) != 4 or tuple(tbkr_shape[2:]) != kr_qc.shape:
        raise ValueError(
            "tbkr shape must have dimensions (t, b, k, r) matching kr QC."
        )

    binary_kr_qc = np.asarray(kr_qc > threshold, dtype=np.uint8)
    return np.broadcast_to(binary_kr_qc, tbkr_shape).astype(
        np.uint8,
        copy=True,
    )


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
            left_indexes[indexes] = _fractional_peak_index(curve, peaks[0])
            left_values[indexes] = curve[peaks[0]]
        if len(peaks) == 2:
            right_indexes[indexes] = _fractional_peak_index(curve, peaks[1])
            right_values[indexes] = curve[peaks[1]]

    return left_indexes, right_indexes, left_values, right_values


def _fractional_peak_index(profile: np.ndarray, peak_index: int) -> np.float32:
    """Refine an integer peak with a three-sample quadratic vertex."""

    integer_position = np.float32(peak_index)
    if peak_index <= 0 or peak_index >= profile.size - 1:
        return integer_position

    left, center, right = np.asarray(
        profile[peak_index - 1 : peak_index + 2],
        dtype=np.float64,
    )
    if not np.all(np.isfinite((left, center, right))):
        return integer_position
    if center < left or center < right:
        return integer_position

    curvature = left - 2.0 * center + right
    if curvature >= 0.0:
        return integer_position
    offset = 0.5 * (left - right) / curvature
    if not np.isfinite(offset) or abs(offset) > 0.5:
        return integer_position
    return np.float32(peak_index + offset)


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
            "temporal_filter": "centered_pixelwise_median",
            "temporal_median_window": np.int32(TEMPORAL_MEDIAN_WINDOW),
            "temporal_median_passes": np.int32(2),
            "temporal_boundary_mode": "replicated_edges",
            "spatial_operator": "3x3 Sobel magnitude",
            "processing_order": (
                "spatial_interpolation, temporal_median, sobel_magnitude, "
                "temporal_median, segment_rotation"
            ),
            "spatial_region": mask,
        }
    )
    value.attrs = attrs
    return value


__all__ = [
    "SPATIAL_GRADIENT_PROFILE_ROOT",
    "SPATIAL_GRADIENT_METRICS_ROOT",
    "SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES",
    "TBKR_LUMEN_SIZE_QC_THRESHOLD",
    "cleanup_spatial_gradient_artifacts",
    "extract_spatial_gradient_segments",
    "pack_spatial_gradient_profile_outputs",
]
