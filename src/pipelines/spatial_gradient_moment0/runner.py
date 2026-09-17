"""Spatial gradient workflow: velocity geometry -> M0ff segments -> filter -> Sobel -> filter.

All spatial gradient settings, processing, exports and metrics live here.
The disabled full-frame legacy workflow is preserved at the bottom.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from types import SimpleNamespace

import numpy as np
from scipy import ndimage

from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    _center_pad_for_rotation,
    _resize_subimage_stack,
    _rotate_stack_with_nan,
    _subimage_values_from_bounds,
)
from calculations.blood_flow_velocity.cross_section.segment_array import SegmentArray
from calculations.math import nanmedian
from pipeline_engine.base import DatasetValue

from input_output.output_manager import OutputType
from input_output.writers.avi import MjpegAviWriter
from pipelines.displacement_map.filtering import CenteredMedianBuffer

# Shared settings for processing and output metadata.
# Each stage independently accepts "none", "median", or "gaussian".
PRE_SOBEL_TEMPORAL_FILTER = "median"
POST_SOBEL_TEMPORAL_FILTER = "gaussian"
TEMPORAL_MEDIAN_WINDOW = 17
TEMPORAL_GAUSSIAN_SIGMA = 2.0  # Frames; no spatial smoothing.
TEMPORAL_GAUSSIAN_TRUNCATE = 4.0
TRANSVERSE_MASK_DILATION_PIXELS = 5
TBKR_LUMEN_SIZE_QC_THRESHOLD = 0.5
SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES = 5
SPATIAL_GRADIENT_PROFILE_ROOT = "Processing/SpatialGradientProfiles"
SPATIAL_GRADIENT_METRICS_ROOT = "Processing/SpatialGradientMetrics"
STATE_KEY = "spatial_gradient_moment0_segments"
EXPORT_DEBUG_SEGMENT_VIDEO = True  # Temporary inspection of the first valid segment.
AVI_FILENAME = "spatial_gradient_moment0_segment_debug.avi"
DEFAULT_FPS = 25.0
CONTRAST_HIGH_PERCENTILE = 99.9
CONTRAST_GAMMA = 1.0


def run_spatial_gradient_moment0(ctx):
    """Process segments from the waveform core and cache them for output packing."""
    ctx.require_inputs("hd")
    from pipelines.waveform_velocity_core.runner import WAVEFORM_CONTEXT_STATE
    context = ctx.state.get(WAVEFORM_CONTEXT_STATE)
    if context is None:
        raise RuntimeError("Spatial gradient requires waveform_velocity_core geometry.")
    segments = extract_spatial_gradient_segments(ctx, context)
    ctx.state.set(STATE_KEY, segments)
    return segments


def spatial_gradient(frame) -> np.ndarray:
    """Return the magnitude of the horizontal and vertical 3x3 Sobel filters."""

    image = np.asarray(frame, dtype=np.float32)
    if image.ndim != 2:
        raise ValueError(f"A moment0 frame must be 2-D, got shape {image.shape}.")
    # ImageJ's Find Edges command combines the two Sobel responses as
    # sqrt(Gx**2 + Gy**2). Nearest-edge extension matches its border handling.
    finite_image = np.nan_to_num(image, nan=0.0, posinf=0.0, neginf=0.0)
    horizontal = ndimage.sobel(finite_image, axis=1, mode="nearest")
    vertical = ndimage.sobel(finite_image, axis=0, mode="nearest")
    return np.hypot(horizontal, vertical, dtype=np.float32)


def temporal_median_filter(
    frames: Iterable[np.ndarray],
    *,
    window: int | None = None,
) -> Iterator[np.ndarray]:
    """Yield centered temporal medians using replicated boundary frames."""

    if window is None:
        window = TEMPORAL_MEDIAN_WINDOW
    if window < 1 or window % 2 == 0:
        raise ValueError("temporal median window must be a positive odd integer.")
    median_buffer = CenteredMedianBuffer(window)
    for frame in frames:
        image = np.asarray(frame, dtype=np.float32)
        if image.ndim != 2:
            raise ValueError(f"A moment0 frame must be 2-D, got shape {image.shape}.")
        yield from median_buffer.push(image)
    yield from median_buffer.finish()


def temporal_gaussian_filter(frames: Iterable[np.ndarray]) -> Iterator[np.ndarray]:
    """Smooth time only, normalizing finite samples and preserving NaN padding."""
    if not math.isfinite(TEMPORAL_GAUSSIAN_SIGMA) or TEMPORAL_GAUSSIAN_SIGMA <= 0:
        raise ValueError("temporal Gaussian sigma must be finite and positive.")
    if not math.isfinite(TEMPORAL_GAUSSIAN_TRUNCATE) or TEMPORAL_GAUSSIAN_TRUNCATE <= 0:
        raise ValueError("temporal Gaussian truncate must be finite and positive.")
    frames = list(frames)
    if not frames:
        return
    stack = np.asarray(frames, dtype=np.float32)
    if stack.ndim != 3:
        raise ValueError("temporal Gaussian input must be a stack of 2-D frames.")
    finite = np.isfinite(stack)
    kwargs = dict(sigma=TEMPORAL_GAUSSIAN_SIGMA, axis=0, mode="nearest",
                  truncate=TEMPORAL_GAUSSIAN_TRUNCATE)
    values = ndimage.gaussian_filter1d(np.where(finite, stack, 0), **kwargs)
    weights = ndimage.gaussian_filter1d(finite.astype(np.float32), **kwargs)
    result = np.full(stack.shape, np.nan, dtype=np.float32)
    np.divide(values, weights, out=result, where=finite & (weights > 0))
    yield from result


# Add future filter implementations here; both stages use the same dispatcher.
TEMPORAL_FILTERS = {
    "none": lambda frames: iter(frames),
    "median": temporal_median_filter,
    "gaussian": temporal_gaussian_filter,
}


def apply_temporal_filter(frames: Iterable[np.ndarray], filter_type: str):
    try:
        implementation = TEMPORAL_FILTERS[filter_type]
    except KeyError:
        raise ValueError(
            f"Unknown temporal filter {filter_type!r}; choose from {tuple(TEMPORAL_FILTERS)}."
        ) from None
    return implementation(frames)


def filtered_segment_gradients(stack):
    """Apply the selected pre-filter, Sobel, and selected post-filter in order."""
    def gradients():
        for frame in apply_temporal_filter(stack, PRE_SOBEL_TEMPORAL_FILTER):
            gradient = spatial_gradient(frame)
            # Padding must not create artificial vessel edges at NaN boundaries.
            valid = ndimage.minimum_filter(np.isfinite(frame), size=3, mode="nearest")
            gradient[~valid] = np.nan
            yield gradient
    yield from apply_temporal_filter(gradients(), POST_SOBEL_TEMPORAL_FILTER)


def _temporal_filter_metadata():
    attrs = {
        "pre_sobel_temporal_filter": PRE_SOBEL_TEMPORAL_FILTER,
        "post_sobel_temporal_filter": POST_SOBEL_TEMPORAL_FILTER,
        "temporal_boundary_mode": "replicated_edges",
        "processing_order": (
            f"extract_resize_rotate_{PRE_SOBEL_TEMPORAL_FILTER}_sobel_"
            f"{POST_SOBEL_TEMPORAL_FILTER}_profile"
        ),
    }
    for stage, filter_type in (("pre_sobel", PRE_SOBEL_TEMPORAL_FILTER),
                               ("post_sobel", POST_SOBEL_TEMPORAL_FILTER)):
        if filter_type == "median":
            attrs[f"{stage}_temporal_median_window"] = TEMPORAL_MEDIAN_WINDOW
        elif filter_type == "gaussian":
            attrs[f"{stage}_temporal_gaussian_sigma_frames"] = TEMPORAL_GAUSSIAN_SIGMA
            attrs[f"{stage}_temporal_gaussian_truncate"] = TEMPORAL_GAUSSIAN_TRUNCATE
    return attrs


def extract_spatial_gradient_segments(ctx, waveform_context):
    """Project M0ff using velocity geometry, then filter and find segment edges."""
    cached = getattr(ctx, "state", None)
    cached = cached.get(STATE_KEY) if cached is not None else None
    if cached is not None:
        return cached
    debug_exported = False
    def export_debug(stack, ring, branch, vessel):
        nonlocal debug_exported
        if (EXPORT_DEBUG_SEGMENT_VIDEO and not debug_exported
                and getattr(getattr(ctx, "output", None), "available", False)):
            _export_debug_segment(ctx, stack, ring, branch, vessel, source_dataset=source_dataset)
            debug_exported = True
    hd = ctx.inputs.hd.as_holodoppler()
    moment0ff = hd.moment0_flat_field_dataset()
    source_dataset = "/moment0ff"
    if moment0ff is None:
        ctx.log_error("Missing flat-field HoloDoppler moment0 dataset (moment0ff/M0FF); falling back to M0.")
        moment0ff = hd.moment0_dataset()
        source_dataset = "/moment0"
        if moment0ff is None:
            raise KeyError("Spatial gradient requires M0ff or fallback M0; neither dataset is available.")
    results = []
    for vessel, segments in (
        ("Artery", waveform_context.artery_segment_result),
        ("Vein", waveform_context.vein_segment_result),
    ):
        debug_callback = None
        if (not debug_exported and EXPORT_DEBUG_SEGMENT_VIDEO
                and getattr(getattr(ctx, "output", None), "available", False)):
            def debug_callback(stack, ring, branch):
                export_debug(stack, ring, branch, vessel)
        results.append(_gradient_segments_from_velocity_geometry(
            moment0ff, segments, debug_callback=debug_callback, source_dataset=source_dataset,
        ))
    return tuple(results)


def _gradient_segments_from_velocity_geometry(
    moment0ff, segments, *, debug_callback=None, source_dataset="/moment0ff",
):
    """Reuse fitted windows, angles and rotated masks without refitting topology."""
    if segments is None:
        return None
    shape = tuple(segments.velocity_profiles.shape)
    if shape[2] != moment0ff.shape[0]:
        raise ValueError("M0ff and velocity segment frame counts must match.")
    unmasked = SegmentArray(shape, np.nan, dtype=np.float32)
    masked = SegmentArray(shape, np.nan, dtype=np.float32)
    topology = segments.topology
    for ring, branch in np.argwhere(topology.valid_segments):
        center = tuple(int(v) for v in segments.segment_center_xy[branch, ring])
        bounds = tuple(int(v) for v in segments.profile_window_bounds_xyxy[ring, branch])
        stack = _subimage_values_from_bounds(
            moment0ff, bounds, loc_xy=center,
            side_pixels=segments.profile_window_side_pixels,
        )
        stack = _rotate_stack_with_nan(
            _center_pad_for_rotation(_resize_subimage_stack(stack), np.nan),
            float(segments.profile_rotation_degrees[ring, branch]),
        )
        mask = ndimage.binary_dilation(
            np.asarray(segments.segment_masks[ring, branch], dtype=bool),
            structure=np.ones((1, 2 * TRANSVERSE_MASK_DILATION_PIXELS + 1), bool),
        )
        debug_stack = np.empty_like(stack) if debug_callback is not None else None
        for time, gradient in enumerate(filtered_segment_gradients(stack)):
            if debug_stack is not None:
                debug_stack[time] = gradient
            unmasked[ring, branch, time] = _mean_transverse(gradient)
            masked[ring, branch, time] = _mean_transverse(
                np.where(mask, gradient, np.nan)
            )
        if debug_stack is not None:
            debug_callback(debug_stack, ring, branch)
            debug_callback = None  # Only retain the first valid segment movie.
    return SimpleNamespace(
        source_dataset=source_dataset,
        velocity_profiles=unmasked,
        transverse_velocity_profiles_masked=masked,
        labels=segments.labels,
        branch_ids=segments.branch_ids,
        segment_center_xy=segments.segment_center_xy,
    )


def _mean_transverse(frame):
    finite = np.isfinite(frame)
    count = finite.sum(axis=0)
    result = np.full(frame.shape[1], np.nan, dtype=np.float32)
    np.divide(np.where(finite, frame, 0).sum(axis=0), count,
              out=result, where=count > 0)
    return result


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
        source_dataset=getattr(segments, "source_dataset", "/moment0ff"),
    )
    masked = _gradient_profile_dataset(
        np.asarray(
            segments.transverse_velocity_profiles_masked,
            dtype=np.float32,
        ),
        cycle_boundary_indexes,
        index_base=index_base,
        mask="vessel_segment",
        source_dataset=getattr(segments, "source_dataset", "/moment0ff"),
    )
    from pipelines.waveform_velocity.profiles import _temporally_meaned_profile_dataset
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
    source_dataset: str = "/moment0ff",
) -> DatasetValue:
    from pipelines.waveform_velocity.profiles import _profile_dataset
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
            "source_dataset": source_dataset,
            **_temporal_filter_metadata(),
            "spatial_operator": "3x3 Sobel magnitude",
            "spatial_region": mask,
            "transverse_mask_dilation_pixels": np.int32(TRANSVERSE_MASK_DILATION_PIXELS),
        }
    )
    value.attrs = attrs
    return value


def _export_debug_segment(ctx, stack, ring, branch, vessel, *, source_dataset="/moment0ff"):
    finite_stack = np.nan_to_num(stack, nan=0.0, posinf=0.0, neginf=0.0)
    maximum = _contrast_maximum(finite_stack.mean(axis=0), float(finite_stack.max()))
    fps = resolve_frame_rate(ctx)
    path = ctx.output.path_for(OutputType.AVI, AVI_FILENAME)
    metadata = {
        "source_dataset": source_dataset, "vessel": vessel,
        "ring_index": int(ring), "branch_index": int(branch),
        **_temporal_filter_metadata(),
        "display_range": [0.0, maximum], "fps": fps,
    }
    with MjpegAviWriter(path, width=stack.shape[2], height=stack.shape[1],
                        fps=fps, metadata=metadata) as video:
        for frame in finite_stack:
            video.write_frame(_display_frame(frame, maximum))
    ctx.log(f"Exported temporary spatial gradient segment video: {path}.")


def _display_frame(gradient: np.ndarray, maximum: float) -> np.ndarray:
    if not math.isfinite(maximum) or maximum <= 0.0:
        return np.zeros(gradient.shape, dtype=np.uint8)
    scaled = np.clip(np.asarray(gradient, dtype=np.float32) / maximum, 0.0, 1.0)
    scaled = np.power(scaled, CONTRAST_GAMMA, dtype=np.float32)
    return np.rint(scaled * 255.0).astype(np.uint8)


def _contrast_maximum(mean_gradient: np.ndarray, fallback: float) -> float:
    values = np.asarray(mean_gradient)
    positive_finite = values[np.isfinite(values) & (values > 0.0)]
    if positive_finite.size:
        maximum = float(
            np.percentile(positive_finite, CONTRAST_HIGH_PERCENTILE)
        )
        if math.isfinite(maximum) and maximum > 0.0:
            return maximum
    return float(fallback)


def resolve_frame_rate(ctx, fallback: float = DEFAULT_FPS) -> float:
    """Resolve the HoloDoppler acquisition rate, with an older-file fallback."""

    try:
        dt_seconds = float(ctx.inputs.hd.as_holodoppler().timing().dt_seconds)
        fps = 1.0 / dt_seconds
        if math.isfinite(fps) and fps > 0.0:
            return fps
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        pass
    fallback = float(fallback)
    if not math.isfinite(fallback) or fallback <= 0.0:
        raise ValueError("fallback fps must be finite and greater than zero.")
    ctx.log_warning(f"Could not resolve HoloDoppler timing; using {fallback:g} fps.")
    return fallback


# LEGACY WORKFLOW (disabled): full-frame median -> Sobel -> segment extraction.
# Preserved for reference; no full-frame AVI or mean PNG is exported.
# @dataclass(frozen=True, slots=True)
# class SpatialGradientMoment0Artifacts:
#     avi_path: Path
#     mean_png_path: Path
#     gradient_path: Path
#     frame_count: int
#     display_maximum: float
#     temporary_directory: object
#
#     def cleanup(self) -> None:
#         cleanup = getattr(self.temporary_directory, "cleanup", None)
#         if cleanup is not None:
#             cleanup()
#
#
# def run_spatial_gradient_moment0(ctx) -> SpatialGradientMoment0Artifacts:
#     """Export a globally scaled gradient AVI and its temporal-mean PNG."""
#
#     ctx.require_inputs("hd")
#     if not ctx.output.available:
#         raise ValueError("An output manager is required to export gradient artifacts.")
#
#     moment0ff = ctx.inputs.hd.as_holodoppler().moment0_flat_field_dataset()
#     if moment0ff is None:
#         raise KeyError(
#             "Missing flat-field HoloDoppler moment0 dataset. Tried: "
#             "'moment0ff', 'M0FF'."
#         )
#     frame_count, height, width = (int(size) for size in moment0ff.shape)
#     if frame_count <= 0 or height <= 0 or width <= 0:
#         raise ValueError(f"moment0ff must have non-empty dimensions, got {moment0ff.shape}.")
#
#     # Preserve the diagnostic float32 cube for the AVI export. Segment profile
#     # analysis reads M0ff directly and applies filtering after resizing/rotation.
#     temporary_directory = tempfile.TemporaryDirectory(
#         prefix=".eyeflow-spatial-gradient-"
#     )
#     gradient_path = Path(temporary_directory.name) / "spatial_gradient_moment0.npy"
#     gradient_video = np.lib.format.open_memmap(
#         gradient_path,
#         mode="w+",
#         dtype=np.float32,
#         shape=(frame_count, height, width),
#     )
#     gradient_sum = np.zeros((height, width), dtype=np.float64)
#     observed_maximum = 0.0
#     filtered_frames = temporal_median_filter(
#         (moment0ff[frame_index] for frame_index in range(frame_count)),
#         window=TEMPORAL_MEDIAN_WINDOW,
#     )
#     filtered_frame_count = 0
#     for frame_index, filtered_frame in enumerate(filtered_frames):
#         gradient = spatial_gradient(filtered_frame)
#         gradient_video[frame_index] = gradient
#         gradient_sum += gradient
#         frame_maximum = float(np.max(gradient, initial=0.0))
#         observed_maximum = max(observed_maximum, frame_maximum)
#         filtered_frame_count += 1
#     if filtered_frame_count != frame_count:
#         raise RuntimeError(
#             "Temporal median filtering changed the spatial-gradient frame count: "
#             f"{filtered_frame_count} != {frame_count}."
#         )
#
#     mean_gradient = (gradient_sum / float(frame_count)).astype(np.float32)
#     gradient_video.flush()
#     display_maximum = _contrast_maximum(mean_gradient, observed_maximum)
#
#     fps = resolve_frame_rate(ctx)
#     avi_path = ctx.output.path_for(OutputType.AVI, AVI_FILENAME)
#     metadata = {
#         "title": "EyeFlow spatial gradient of moment0ff",
#         "artifact": "spatial_gradient_moment0",
#         "algorithm": (
#             "Centered temporal median followed by Sobel gradient magnitude "
#             "(ImageJ/Fiji Find Edges)"
#         ),
#         "source_dataset": "/moment0ff",
#         "temporal_median_window": TEMPORAL_MEDIAN_WINDOW,
#         "temporal_median": "centered pixel-wise median with replicated edges",
#         "display_range": [0.0, display_maximum],
#         "contrast_high_percentile": CONTRAST_HIGH_PERCENTILE,
#         "contrast_gamma": CONTRAST_GAMMA,
#         "fps": fps,
#         "frame_count": frame_count,
#     }
#     with MjpegAviWriter(
#         avi_path,
#         width=width,
#         height=height,
#         fps=fps,
#         metadata=metadata,
#     ) as video:
#         display_sum = np.zeros((height, width), dtype=np.float64)
#         for frame_index in range(frame_count):
#             display_frame = _display_frame(
#                 gradient_video[frame_index], display_maximum
#             )
#             display_sum += display_frame
#             video.write_frame(display_frame)
#
#     mean_display_frame = np.rint(display_sum / float(frame_count)).astype(np.uint8)
#     mean_png_path = ctx.output.write_png(
#         mean_display_frame,
#         PNG_FILENAME,
#     )
#     artifacts = SpatialGradientMoment0Artifacts(
#         avi_path=avi_path,
#         mean_png_path=mean_png_path,
#         gradient_path=gradient_path,
#         frame_count=frame_count,
#         display_maximum=display_maximum,
#         temporary_directory=temporary_directory,
#     )
#     del gradient_video
#     ctx.state.set(STATE_KEY, artifacts)
#     ctx.log(
#         "Exported moment0 spatial gradients: "
#         f"{avi_path} ({frame_count} frames) and {mean_png_path}."
#     )
#     return artifacts
#
#
