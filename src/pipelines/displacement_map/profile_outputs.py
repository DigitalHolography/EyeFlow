"""HDF5 packing for displacement profiles, metrics, and segment maps."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import find_peaks

from calculations.math import nanmean_float32
from calculations.vessel_segments.profiles.per_beat import (
    interpolate_profiles_per_beat,
    interpolate_segment_maps_per_beat,
)
from input_output.schema import EyeFlowOutputPaths
from input_output.writers.h5_storage import profile_h5_options, velocity_map_h5_options
from pipeline_engine.planning.base import DatasetValue
from pipelines.profile_outputs import profile_dataset, temporally_meaned_profile_dataset

from .constants import registration_method_output_name
from .segments import DisplacementSegmentResult

_PROFILE_FIELDS = (
    ("X", "x_sum_profile", "local_x"),
    ("Y", "y_sum_profile", "local_y"),
)
_VESSEL_OUTPUT_NAMES = {"artery": "Artery", "vein": "Vein"}


def pack_displacement_profile_outputs(
    vessels: Mapping[str, DisplacementSegmentResult],
    registration_method: str,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
    timing_attrs: Mapping[str, object] | None = None,
    output_paths: EyeFlowOutputPaths | None = None,
) -> dict[str, DatasetValue]:
    """Pack every restored profile and metric under the final shared roots."""

    schema = output_paths or EyeFlowOutputPaths.active()
    method = registration_method_output_name(registration_method)
    common_attrs = {
        "registration_method": str(registration_method),
        **dict(timing_attrs or {}),
    }
    outputs: dict[str, DatasetValue] = {}
    for vessel, result in sorted(vessels.items()):
        vessel_name = _vessel_output_name(vessel)
        profile_root = f"{schema.displacement_profiles_root}/{method}/{vessel_name}"
        metrics_root = f"{schema.displacement_metrics_root}/{method}/{vessel_name}"
        outputs.update(
            _pack_vessel_profiles(
                result,
                profile_root,
                metrics_root,
                cycle_boundary_indexes,
                index_base=index_base,
                common_attrs=common_attrs,
            )
        )
    return outputs


def pack_displacement_segment_map_outputs(
    vessels: Mapping[str, DisplacementSegmentResult],
    registration_method: str,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
    timing_attrs: Mapping[str, object] | None = None,
    output_paths: EyeFlowOutputPaths | None = None,
    array_allocator=None,
) -> dict[str, DatasetValue]:
    """Pack explicitly requested dense per-segment displacement maps."""

    schema = output_paths or EyeFlowOutputPaths.active()
    method = registration_method_output_name(registration_method)
    outputs: dict[str, DatasetValue] = {}
    for vessel, result in sorted(vessels.items()):
        if result.maps is None:
            raise RuntimeError(
                "Per-segment displacement maps were not retained; enable the "
                "segment_maps option before displacement analysis."
            )
        maps = interpolate_segment_maps_per_beat(
            result.maps,
            cycle_boundary_indexes,
            index_base=index_base,
            valid_segments=result.valid_segments,
            array_allocator=array_allocator,
        )
        path = (
            f"{schema.displacement_map_root}/{method}/"
            f"{_vessel_output_name(vessel)}/PerSegment"
        )
        outputs[path] = DatasetValue(
            data=maps,
            attrs={
                "unit": "pixels",
                "dimDesc": [
                    "x",
                    "y",
                    "time",
                    "beat",
                    "branch",
                    "radius",
                    "displacement_orientation",
                ],
                "coordinate_system": "rotated_segment_pixel",
                "components": ["local_x", "local_y"],
                "component_basis": "rotated_segment_local",
                "registration_method": str(registration_method),
                **dict(timing_attrs or {}),
            },
            h5_options={
                **velocity_map_h5_options(maps.shape),
                "sparse_nan_chunks": True,
            },
        )
    return outputs


def _pack_vessel_profiles(
    result: DisplacementSegmentResult,
    profile_root: str,
    metrics_root: str,
    cycle_boundary_indexes,
    *,
    index_base: int,
    common_attrs: Mapping[str, object],
) -> dict[str, DatasetValue]:
    outputs: dict[str, DatasetValue] = {}
    for axis_name, field_name, component in _PROFILE_FIELDS:
        outputs[f"{profile_root}/{axis_name}_sum_displacement_profile"] = (
            _summed_displacement_profile_dataset(
                np.asarray(getattr(result, field_name), dtype=np.float32),
                cycle_boundary_indexes,
                index_base=index_base,
                component=component,
                attrs=common_attrs,
            )
        )
    outputs[f"{profile_root}/DisplacementMagnitude"] = (
        _segment_displacement_magnitude_dataset(
            result.x_sum_profile,
            result.y_sum_profile,
            cycle_boundary_indexes,
            index_base=index_base,
            attrs=common_attrs,
        )
    )
    outputs[f"{profile_root}/Cross_sectional_radial_movement_amplitude_profile"] = (
        _segment_displacement_metric_dataset(
            result.radial_movement_amplitude,
            cycle_boundary_indexes,
            index_base=index_base,
            unit="pixels",
            attrs=common_attrs,
            metric_attrs={
                "component": "absolute_local_x",
                "centerline_source": "rotated_vessel_segment_mask",
                "spatial_region": "symmetric_wall_bands_extending_outside_mask",
                "spatial_reduction": "mean_per_side_then_mean",
            },
        )
    )
    outputs[f"{profile_root}/Cross_sectional_radial_asymmetry_index_profile"] = (
        _segment_displacement_metric_dataset(
            result.radial_asymmetry_index,
            cycle_boundary_indexes,
            index_base=index_base,
            unit="1",
            attrs=common_attrs,
            metric_attrs={
                "component": "absolute_local_x",
                "centerline_source": "rotated_vessel_segment_mask",
                "spatial_region": "symmetric_wall_bands_extending_outside_mask",
                "spatial_reduction": (
                    "(left_strength-right_strength)/"
                    "(left_strength+right_strength+1e-6)"
                ),
                "side_convention": "left_right_in_rotated_segment_coordinates",
            },
        )
    )

    derived: dict[tuple[str, str], tuple[DatasetValue, DatasetValue]] = {}
    for direction, spatial_axis in (("Longitudinal", "y"), ("Transverse", "x")):
        for mask_name in ("Masked", "Unmasked"):
            field_name = f"{direction.lower()}_profiles_{mask_name.lower()}"
            profile = profile_dataset(
                np.asarray(getattr(result, field_name), dtype=np.float32),
                cycle_boundary_indexes,
                index_base=index_base,
                spatial_axis=spatial_axis,
                unit="pixels",
                attrs=common_attrs,
            )
            meaned = temporally_meaned_profile_dataset(profile)
            squared_deviation = _temporally_centered_profile_power_dataset(
                profile,
                meaned,
            )
            meaned_squared_deviation = temporally_meaned_profile_dataset(
                squared_deviation
            )
            root = f"{profile_root}/{direction}/{mask_name}"
            outputs.update(
                {
                    f"{root}/tbkr/Profile": profile,
                    f"{root}/bkr/Profile": meaned,
                    f"{root}/b/Profile": _globally_meaned_profile_dataset(meaned),
                    f"{root}/tbkr/SquaredDeviation": squared_deviation,
                    f"{root}/bkr/SquaredDeviation": meaned_squared_deviation,
                }
            )
            derived[(direction, mask_name)] = (meaned, meaned_squared_deviation)

    meaned, temporal_variance = derived[("Transverse", "Masked")]
    x_max = _transverse_max_x_position_dataset(meaned)
    y_max = _mean_profile_values_at_peak_positions(meaned, x_max)
    y_variance = _mean_profile_values_at_peak_positions(
        temporal_variance,
        x_max,
        value_order=(
            "first_peak_temporal_variance",
            "second_peak_temporal_variance",
        ),
    )
    area_left, area_right = _profile_peak_area_datasets(meaned, x_max)
    variance_left, variance_right = _profile_peak_area_datasets(
        temporal_variance,
        x_max,
    )
    position_root = f"{metrics_root}/Transverse/Position"
    area_root = f"{metrics_root}/Transverse/Area"
    outputs.update(
        {
            f"{position_root}/XMax": x_max,
            f"{position_root}/YMax": y_max,
            f"{position_root}/YMaxTemporalVariance": y_variance,
            f"{position_root}/XMaxMeaned": _mean_peak_metric(x_max),
            f"{position_root}/YMaxMeaned": _mean_peak_metric(y_max),
            f"{position_root}/YMaxMeanedTemporalVariance": _mean_peak_metric(
                y_variance
            ),
            f"{area_root}/Left": area_left,
            f"{area_root}/Right": area_right,
            f"{area_root}/LeftTemporalVariance": variance_left,
            f"{area_root}/RightTemporalVariance": variance_right,
        }
    )
    return outputs


def _vessel_output_name(vessel: str) -> str:
    try:
        return _VESSEL_OUTPUT_NAMES[str(vessel).lower()]
    except KeyError as exc:
        raise ValueError(f"Unsupported displacement vessel: {vessel!r}.") from exc


def _options(data: np.ndarray) -> dict[str, object]:
    return {**profile_h5_options(data.shape), "sparse_nan_chunks": True}


def _provenance(attrs: Mapping[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in attrs.items()
        if key.startswith("cycle_") or key == "registration_method"
    }


def _globally_meaned_profile_dataset(profile: DatasetValue) -> DatasetValue:
    values = np.asarray(profile.data, dtype=np.float32)
    attrs = dict(profile.attrs or {})
    dimensions = list(attrs.get("dimDesc", ()))
    if values.ndim != 4 or dimensions[1:] != ["beat", "branch", "radius"]:
        raise ValueError(
            "Profile dataset must have spatial, beat, branch, and radius dimensions."
        )
    data = nanmean_float32(values, axis=(2, 3))
    attrs["dimDesc"] = dimensions[:2]
    attrs["segment_reduction"] = "mean_over_valid_branch_radius_segments"
    return DatasetValue(data=data, attrs=attrs, h5_options=_options(data))


def _transverse_max_x_position_dataset(profile: DatasetValue) -> DatasetValue:
    values = np.asarray(profile.data, dtype=np.float32)
    attrs = dict(profile.attrs or {})
    dimensions = list(attrs.get("dimDesc", ()))
    if values.ndim != 4 or dimensions[1:] != ["beat", "branch", "radius"]:
        raise ValueError(
            "Time-meaned profile must have spatial, beat, branch, and radius dimensions."
        )
    positions = np.full((2, *values.shape[1:]), np.nan, dtype=np.float32)
    for beat, branch, radius in np.ndindex(values.shape[1:]):
        peaks = _first_two_peak_positions(values[:, beat, branch, radius])
        positions[: peaks.size, beat, branch, radius] = peaks
    return DatasetValue(
        data=positions,
        attrs={
            "unit": "pixels",
            "dimDesc": ["peak", "beat", "branch", "radius"],
            "index_base": np.int32(0),
            "position_axis": dimensions[0],
            "position_representation": "zero_based_profile_sample_index",
            "value_order": ["first_peak_x", "second_peak_x"],
            "peak_selection": "first_two_local_maxima_left_to_right",
            "source_temporal_reduction": attrs.get("temporal_reduction", ""),
            **_provenance(attrs),
        },
        h5_options=_options(positions),
    )


def _mean_profile_values_at_peak_positions(
    profile: DatasetValue,
    peak_positions: DatasetValue,
    *,
    value_order: tuple[str, str] = ("first_peak_y", "second_peak_y"),
) -> DatasetValue:
    values, positions, spatial_axis = _validate_peak_sources(profile, peak_positions)
    extracted = np.full(positions.shape, np.nan, dtype=np.float32)
    for peak, beat, branch, radius in np.ndindex(positions.shape):
        position = positions[peak, beat, branch, radius]
        if np.isfinite(position) and 0 <= int(position) < values.shape[0]:
            extracted[peak, beat, branch, radius] = values[
                int(position), beat, branch, radius
            ]
    attrs = dict(profile.attrs or {})
    return DatasetValue(
        data=extracted,
        attrs={
            "unit": attrs.get("unit", "pixels"),
            "dimDesc": ["peak", "beat", "branch", "radius"],
            "peak_position_axis": spatial_axis,
            "peak_position_index_base": np.int32(0),
            "value_order": list(value_order),
            "source_temporal_reduction": attrs.get("temporal_reduction", ""),
            **_provenance(attrs),
        },
        h5_options=_options(extracted),
    )


def _profile_peak_area_datasets(
    profile: DatasetValue,
    peak_positions: DatasetValue,
) -> tuple[DatasetValue, DatasetValue]:
    values = np.asarray(profile.data, dtype=np.float32)
    attrs = dict(profile.attrs or {})
    dimensions = list(attrs.get("dimDesc", ()))
    positions = np.asarray(peak_positions.data, dtype=np.float32)
    if (
        values.ndim != 4
        or dimensions[1:] != ["beat", "branch", "radius"]
        or positions.shape != (2, *values.shape[1:])
    ):
        raise ValueError("Area sources have incompatible profile and peak dimensions.")
    left = np.full(values.shape[1:], np.nan, dtype=np.float32)
    right = np.full(values.shape[1:], np.nan, dtype=np.float32)
    for beat, branch, radius in np.ndindex(values.shape[1:]):
        first, second = positions[:, beat, branch, radius]
        if not np.isfinite(first) or not np.isfinite(second) or first >= second:
            continue
        curve = values[:, beat, branch, radius]
        finite_indexes = np.flatnonzero(np.isfinite(curve))
        if finite_indexes.size < 2:
            continue
        left_edge = float(finite_indexes[0])
        right_edge = float(finite_indexes[-1])
        midpoint = (float(first) + float(second)) / 2.0
        if left_edge < midpoint < right_edge:
            left[beat, branch, radius] = _finite_trapezoid_area(
                curve, left_edge, midpoint
            )
            right[beat, branch, radius] = _finite_trapezoid_area(
                curve, midpoint, right_edge
            )
    common = {
        "unit": _integrated_profile_unit(attrs.get("unit", "")),
        "dimDesc": ["beat", "branch", "radius"],
        "integration_axis": dimensions[0],
        "integration_method": "trapezoidal",
        "outer_boundaries": "first_and_last_finite_profile_samples",
        "shared_boundary": "midpoint_between_XMax_values",
        **_provenance(attrs),
    }
    return (
        DatasetValue(
            data=left,
            attrs={**common, "side": "left"},
            h5_options=_options(left),
        ),
        DatasetValue(
            data=right,
            attrs={**common, "side": "right"},
            h5_options=_options(right),
        ),
    )


def _finite_trapezoid_area(
    curve: np.ndarray,
    start: float,
    stop: float,
) -> np.float32:
    values = np.asarray(curve, dtype=np.float32)
    finite_indexes = np.flatnonzero(np.isfinite(values))
    if finite_indexes.size < 2 or stop <= start:
        return np.float32(np.nan)
    finite_runs = np.split(
        finite_indexes,
        np.flatnonzero(np.diff(finite_indexes) > 1) + 1,
    )
    total = 0.0
    integrated = False
    for run in finite_runs:
        run_start = max(start, float(run[0]))
        run_stop = min(stop, float(run[-1]))
        if run.size < 2 or run_stop <= run_start:
            continue
        interior = run[(run > run_start) & (run < run_stop)].astype(np.float64)
        x_values = np.concatenate(([run_start], interior, [run_stop]))
        y_values = np.interp(x_values, run, values[run])
        total += float(trapezoid(y_values, x_values))
        integrated = True
    return np.float32(total) if integrated else np.float32(np.nan)


def _integrated_profile_unit(unit: object) -> str:
    source = str(unit)
    if source == "pixels":
        return "pixels^2"
    if source == "pixels^2":
        return "pixels^3"
    return f"{source}*pixels" if source else "pixels"


def _validate_peak_sources(
    profile: DatasetValue,
    peak_positions: DatasetValue,
) -> tuple[np.ndarray, np.ndarray, str]:
    values = np.asarray(profile.data, dtype=np.float32)
    dimensions = list((profile.attrs or {}).get("dimDesc", ()))
    positions = np.asarray(peak_positions.data, dtype=np.float32)
    position_dimensions = list((peak_positions.attrs or {}).get("dimDesc", ()))
    if values.ndim != 4 or dimensions[1:] != ["beat", "branch", "radius"]:
        raise ValueError("Time-meaned profile dimensions are invalid.")
    if (
        positions.shape != (2, *values.shape[1:])
        or position_dimensions != ["peak", "beat", "branch", "radius"]
    ):
        raise ValueError("Peak position dimensions are invalid.")
    return values, positions, dimensions[0]


def _mean_peak_metric(metric: DatasetValue) -> DatasetValue:
    values = np.asarray(metric.data, dtype=np.float32)
    attrs = dict(metric.attrs or {})
    if (
        values.ndim != 4
        or list(attrs.get("dimDesc", ()))
        != ["peak", "beat", "branch", "radius"]
    ):
        raise ValueError("Peak metric dimensions are invalid.")
    data = nanmean_float32(values, axis=(2, 3))
    attrs["dimDesc"] = ["peak", "beat"]
    attrs["aggregation"] = "mean_over_valid_branch_radius_values"
    return DatasetValue(data=data, attrs=attrs, h5_options=_options(data))


def _first_two_peak_positions(curve: np.ndarray) -> np.ndarray:
    values = np.asarray(curve, dtype=np.float32)
    finite_indexes = np.flatnonzero(np.isfinite(values))
    if finite_indexes.size < 3:
        return np.empty(0, dtype=np.float32)
    finite_runs = np.split(
        finite_indexes,
        np.flatnonzero(np.diff(finite_indexes) > 1) + 1,
    )
    positions: list[int] = []
    for run in finite_runs:
        if run.size >= 3:
            local_positions, _ = find_peaks(values[run])
            positions.extend(int(run[index]) for index in local_positions)
        if len(positions) >= 2:
            break
    return np.asarray(positions[:2], dtype=np.float32)


def _temporally_centered_profile_power_dataset(
    profile: DatasetValue,
    temporal_mean: DatasetValue,
) -> DatasetValue:
    values = np.asarray(profile.data, dtype=np.float32)
    mean = np.asarray(temporal_mean.data, dtype=np.float32)
    if mean.shape != (values.shape[0], *values.shape[2:]):
        raise ValueError("Temporal mean must match the profile without its time axis.")
    data = np.square(values - mean[:, None, ...]).astype(np.float32, copy=False)
    attrs = dict(profile.attrs or {})
    attrs.update(
        {
            "unit": "pixels^2",
            "formula": "(D(t) - mean_t(D(t)))**2",
            "temporal_centering": "per_beat_mean",
        }
    )
    return DatasetValue(data=data, attrs=attrs, h5_options=_options(data))


def _summed_displacement_profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    component: str,
    attrs: Mapping[str, object],
) -> DatasetValue:
    values = np.asarray(profiles, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError(
            "Summed displacement profiles must have shape (radius, branch, frame)."
        )
    per_beat = interpolate_profiles_per_beat(
        values[..., None],
        cycle_boundary_indexes,
        index_base=index_base,
    )
    data = per_beat[0]
    return DatasetValue(
        data=data,
        attrs={
            "unit": "pixels",
            "dimDesc": ["time", "beat", "branch", "radius"],
            "coordinate_system": "rotated_segment_pixel",
            "component_basis": "rotated_segment_local",
            "component": component,
            "spatial_region": "full_unmasked_subimage",
            "spatial_reduction": "sum_over_valid_subimage_pixels",
            "displacement_reference": "temporal_mean_image",
            **dict(attrs),
        },
        h5_options=_options(data),
    )


def _segment_displacement_metric_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    unit: str,
    attrs: Mapping[str, object],
    metric_attrs: Mapping[str, object],
) -> DatasetValue:
    values = np.asarray(profiles, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError(
            "Segment displacement metrics must have shape (radius, branch, frame)."
        )
    per_beat = interpolate_profiles_per_beat(
        values[..., None],
        cycle_boundary_indexes,
        index_base=index_base,
    )
    data = per_beat[0]
    return DatasetValue(
        data=data,
        attrs={
            "unit": unit,
            "dimDesc": ["time", "beat", "branch", "radius"],
            "coordinate_system": "rotated_segment_local",
            "displacement_reference": "temporal_mean_image",
            **dict(metric_attrs),
            **dict(attrs),
        },
        h5_options=_options(data),
    )


def _segment_displacement_magnitude_dataset(
    x_profiles: np.ndarray,
    y_profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    attrs: Mapping[str, object],
) -> DatasetValue:
    x_values = np.asarray(x_profiles, dtype=np.float32)
    y_values = np.asarray(y_profiles, dtype=np.float32)
    if x_values.ndim != 3 or x_values.shape != y_values.shape:
        raise ValueError(
            "X and Y displacement profiles must have matching "
            "(radius, branch, frame) shapes."
        )
    magnitude = np.hypot(x_values, y_values).astype(np.float32, copy=False)
    per_beat = interpolate_profiles_per_beat(
        magnitude[..., None],
        cycle_boundary_indexes,
        index_base=index_base,
    )
    data = per_beat[0]
    return DatasetValue(
        data=data,
        attrs={
            "unit": "pixels",
            "dimDesc": ["time", "beat", "branch", "radius"],
            "coordinate_system": "rotation_invariant",
            "magnitude_formula": "sqrt(x**2 + y**2)",
            "spatial_region": "vessel_segment",
            "spatial_reduction": "magnitude_of_summed_displacement_components",
            "displacement_reference": "temporal_mean_image",
            **dict(attrs),
        },
        h5_options=_options(data),
    )


__all__ = [
    "pack_displacement_profile_outputs",
    "pack_displacement_segment_map_outputs",
]
