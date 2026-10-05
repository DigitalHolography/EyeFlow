"""HDF5 export for transverse and longitudinal velocity/displacement profiles."""

from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import find_peaks

from calculations.math import nanmean_float32
from calculations.topology import dilate_segment_masks, interpolate_profiles_per_beat
from input_output.writers.h5 import profile_h5_options
from input_output.schema import EyeFlowOutputPaths, VelocityProfileOutputPaths
from pipeline_engine.base import DatasetValue
from pipelines.displacement_map.constants import registration_method_output_name
from pipelines.shared.profile_datasets import (
    _profile_dataset,
    _temporally_meaned_profile_dataset,
)

_PROFILE_MASK_DILATION_ITERATIONS = 10
_DISPLACEMENT_PROFILE_ROOT = "Processing/Displacement/Profiles"
_DISPLACEMENT_METRICS_ROOT = "Processing/Displacement/Metrics"
_DISPLACEMENT_PROFILE_FIELDS = (
    ("X", "x_sum_profile", "local_x"),
    ("Y", "y_sum_profile", "local_y"),
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


def pack_displacement_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    """Pack per-segment and vessel-level displacement profiles."""

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


def pack_displacement_magnitude_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    """Pack one displacement-magnitude waveform per vessel segment."""

    outputs = _pack_vessel_displacement_magnitudes(
        artery_segments,
        "Artery",
        cycle_boundary_indexes,
        index_base=index_base,
    )
    outputs.update(
        _pack_vessel_displacement_magnitudes(
            vein_segments,
            "Vein",
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return outputs


def pack_cross_section_displacement_profile_outputs(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> dict[str, object]:
    """Pack displacement-magnitude profiles along each cross-section axis."""

    outputs = _pack_vessel_displacement_axis_profiles(
        artery_segments,
        "Artery",
        cycle_boundary_indexes,
        index_base=index_base,
    )
    outputs.update(
        _pack_vessel_displacement_axis_profiles(
            vein_segments,
            "Vein",
            cycle_boundary_indexes,
            index_base=index_base,
        )
    )
    return outputs


def _pack_vessel_displacement_axis_profiles(
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
        method = registration_method_output_name(raw_method)
        root = f"{_DISPLACEMENT_PROFILE_ROOT}/{method}/{vessel_name}"
        metrics_root = f"{_DISPLACEMENT_METRICS_ROOT}/{method}/{vessel_name}"
        outputs.update(
            _pack_displacement_axis_profiles_for_method(
                displacement,
                root,
                metrics_root=metrics_root,
                cycle_boundary_indexes=cycle_boundary_indexes,
                index_base=index_base,
            )
        )
    return outputs


def _pack_displacement_axis_profiles_for_method(
    displacement,
    root: str,
    *,
    metrics_root: str,
    cycle_boundary_indexes,
    index_base: int,
) -> dict[str, object]:
    outputs: dict[str, object] = {}
    derived_profiles: dict[tuple[str, str], tuple[DatasetValue, DatasetValue]] = {}
    for direction_name, spatial_axis in (
        ("Longitudinal", "y"),
        ("Transverse", "x"),
    ):
        for mask_name in ("Masked", "Unmasked"):
            profile_field = (
                f"{direction_name.lower()}_profiles_{mask_name.lower()}"
            )
            profile = _profile_dataset(
                np.asarray(getattr(displacement, profile_field), dtype=np.float32),
                cycle_boundary_indexes,
                index_base=index_base,
                spatial_axis=spatial_axis,
                unit="pixels",
            )
            meaned_profile = _temporally_meaned_profile_dataset(profile)
            displacement_power = _temporally_centered_profile_power_dataset(
                profile,
                meaned_profile,
            )
            meaned_displacement_power = _temporally_meaned_profile_dataset(
                displacement_power
            )
            profile_root = f"{root}/{direction_name}/{mask_name}"
            outputs.update(
                {
                    f"{profile_root}/tbkr/Profile": profile,
                    f"{profile_root}/bkr/Profile": meaned_profile,
                    f"{profile_root}/b/Profile": (
                        _globally_meaned_profile_dataset(meaned_profile)
                    ),
                    f"{profile_root}/tbkr/SquaredDeviation": (
                        displacement_power
                    ),
                    f"{profile_root}/bkr/SquaredDeviation": (
                        meaned_displacement_power
                    ),
                }
            )
            derived_profiles[(direction_name, mask_name)] = (
                meaned_profile,
                meaned_displacement_power,
            )

    transverse_meaned, transverse_temporal_variance = derived_profiles[
        ("Transverse", "Masked")
    ]
    max_x_position = _transverse_max_x_position_dataset(transverse_meaned)
    max_y_position = _mean_profile_values_at_peak_positions(
        transverse_meaned,
        max_x_position,
    )
    temporal_variance_at_peaks = _mean_profile_values_at_peak_positions(
        transverse_temporal_variance,
        max_x_position,
        value_order=(
            "first_peak_temporal_variance",
            "second_peak_temporal_variance",
        ),
    )
    area_l, area_r = _profile_peak_area_datasets(
        transverse_meaned,
        max_x_position,
    )
    temporal_variance_area_l, temporal_variance_area_r = (
        _profile_peak_area_datasets(
            transverse_temporal_variance,
            max_x_position,
        )
    )
    transverse_metrics_root = f"{metrics_root}/Transverse"
    position_metrics_root = f"{transverse_metrics_root}/Position"
    area_metrics_root = f"{transverse_metrics_root}/Area"
    outputs.update(
        {
            f"{position_metrics_root}/XMax": max_x_position,
            f"{position_metrics_root}/YMax": max_y_position,
            f"{position_metrics_root}/YMaxTemporalVariance": (
                temporal_variance_at_peaks
            ),
            f"{position_metrics_root}/XMaxMeaned": (
                _mean_peak_metric(max_x_position)
            ),
            f"{position_metrics_root}/YMaxMeaned": (
                _mean_peak_metric(max_y_position)
            ),
            f"{position_metrics_root}/YMaxMeanedTemporalVariance": (
                _mean_peak_metric(temporal_variance_at_peaks)
            ),
            f"{area_metrics_root}/Left": area_l,
            f"{area_metrics_root}/Right": area_r,
            f"{area_metrics_root}/LeftTemporalVariance": (
                temporal_variance_area_l
            ),
            f"{area_metrics_root}/RightTemporalVariance": (
                temporal_variance_area_r
            ),
        }
    )
    return outputs


def _pack_vessel_displacement_magnitudes(
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
        method = registration_method_output_name(raw_method)
        path = (
            f"{_DISPLACEMENT_PROFILE_ROOT}/{method}/{vessel_name}/"
            "DisplacementMagnitude"
        )
        outputs[path] = _segment_displacement_magnitude_dataset(
            np.asarray(
                displacement.x_sum_profile,
                dtype=np.float32,
            ),
            np.asarray(
                displacement.y_sum_profile,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
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
        method = registration_method_output_name(raw_method)
        root = f"{_DISPLACEMENT_PROFILE_ROOT}/{method}/{vessel_name}"
        for axis_name, profile_field, component in _DISPLACEMENT_PROFILE_FIELDS:
            outputs[f"{root}/{axis_name}_sum_displacement_profile"] = (
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
        outputs[
            f"{root}/Cross_sectional_radial_movement_amplitude_profile"
        ] = _segment_displacement_metric_dataset(
            np.asarray(
                displacement.radial_movement_amplitude,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            unit="pixels",
            metric_attrs={
                "component": "absolute_local_x",
                "centerline_source": "rotated_vessel_segment_mask",
                "spatial_region": "symmetric_wall_bands_extending_outside_mask",
                "spatial_reduction": "mean_per_side_then_mean",
            },
        )
        outputs[
            f"{root}/Cross_sectional_radial_asymmetry_index_profile"
        ] = _segment_displacement_metric_dataset(
            np.asarray(
                displacement.radial_asymmetry_index,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            unit="1",
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
    return outputs


def _pack_vessel_profiles(
    paths: VelocityProfileOutputPaths,
    segments,
    cycle_boundary_indexes,
    *,
    index_base: int,
    include_temporal_means: bool = False,
) -> dict[str, object]:
    valid_segments = np.asarray(segments.topology.valid_segments, dtype=bool)
    transverse_masked = np.asarray(
        segments.transverse_profiles_masked,
        dtype=np.float32,
    )
    outputs = {
        paths.transverse_velocity_profile_unmasked: _profile_dataset(
            np.asarray(segments.transverse_profiles_unmasked, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
            valid_segments=valid_segments,
        ),
        paths.transverse_velocity_profile_masked: _profile_dataset(
            transverse_masked,
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
            valid_segments=valid_segments,
        ),
        paths.longitudinal_velocity_profile_unmasked: _profile_dataset(
            np.asarray(
                segments.longitudinal_profiles_unmasked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
            valid_segments=valid_segments,
        ),
        paths.longitudinal_velocity_profile_masked: _profile_dataset(
            np.asarray(
                segments.longitudinal_profiles_masked,
                dtype=np.float32,
            ),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="y",
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


def _globally_meaned_profile_dataset(profile: DatasetValue) -> DatasetValue:
    """Merge time-meaned branch/radius profiles into one profile per beat."""

    values = np.asarray(profile.data)
    attrs = dict(profile.attrs or {})
    dim_desc = list(attrs.get("dimDesc", ()))
    if (
        values.ndim != 4
        or len(dim_desc) != 4
        or dim_desc[1:] != ["beat", "branch", "radius"]
    ):
        raise ValueError(
            "profile dataset must have dimensions "
            "(spatial_sample, beat, branch, radius)."
        )
    data = nanmean_float32(values, axis=(2, 3))
    attrs["dimDesc"] = dim_desc[:2]
    attrs["segment_reduction"] = "mean_over_valid_branch_radius_segments"
    return DatasetValue(
        data=data,
        attrs=attrs,
        h5_options=profile_h5_options(data.shape),
    )


def _transverse_max_x_position_dataset(profile: DatasetValue) -> DatasetValue:
    """Store the first two transverse-profile local maxima from left to right."""

    values = np.asarray(profile.data, dtype=np.float32)
    attrs = dict(profile.attrs or {})
    dim_desc = list(attrs.get("dimDesc", ()))
    if (
        values.ndim != 4
        or len(dim_desc) != 4
        or dim_desc[1:] != ["beat", "branch", "radius"]
    ):
        raise ValueError(
            "time-meaned profile dataset must have dimensions "
            "(spatial_sample, beat, branch, radius)."
        )

    positions = np.full((2, *values.shape[1:]), np.nan, dtype=np.float32)
    for beat_index, branch_index, radius_index in np.ndindex(values.shape[1:]):
        peaks = _first_two_peak_positions(
            values[:, beat_index, branch_index, radius_index]
        )
        positions[: peaks.size, beat_index, branch_index, radius_index] = peaks

    return DatasetValue(
        data=positions,
        attrs={
            "unit": "pixels",
            "dimDesc": ["peak", "beat", "branch", "radius"],
            "index_base": np.int32(0),
            "position_axis": dim_desc[0],
            "position_representation": "zero_based_profile_sample_index",
            "value_order": ["first_peak_x", "second_peak_x"],
            "peak_selection": "first_two_local_maxima_left_to_right",
            "source_temporal_reduction": attrs.get("temporal_reduction", ""),
        },
        h5_options=profile_h5_options(positions.shape),
    )


def _mean_profile_values_at_peak_positions(
    profile: DatasetValue,
    peak_positions: DatasetValue,
    *,
    value_order: tuple[str, str] = ("first_peak_y", "second_peak_y"),
) -> DatasetValue:
    """Extract both peak values from a time-meaned profile."""

    values, positions, spatial_axis = _validate_meaned_peak_metric_sources(
        profile,
        peak_positions,
    )
    extracted = np.full(positions.shape, np.nan, dtype=np.float32)
    for peak_index, beat_index, branch_index, radius_index in np.ndindex(
        positions.shape
    ):
        position = positions[peak_index, beat_index, branch_index, radius_index]
        if not np.isfinite(position):
            continue
        spatial_index = int(position)
        if spatial_index < 0 or spatial_index >= values.shape[0]:
            continue
        extracted[peak_index, beat_index, branch_index, radius_index] = (
            values[spatial_index, beat_index, branch_index, radius_index]
        )

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
        },
        h5_options=profile_h5_options(extracted.shape),
    )


def _profile_peak_area_datasets(
    profile: DatasetValue,
    peak_positions: DatasetValue,
) -> tuple[DatasetValue, DatasetValue]:
    """Integrate the left and right curve regions split between both peaks."""

    values = np.asarray(profile.data, dtype=np.float32)
    profile_attrs = dict(profile.attrs or {})
    profile_dims = list(profile_attrs.get("dimDesc", ()))
    positions = np.asarray(peak_positions.data, dtype=np.float32)
    position_dims = list((peak_positions.attrs or {}).get("dimDesc", ()))
    if (
        values.ndim != 4
        or len(profile_dims) != 4
        or profile_dims[1:] != ["beat", "branch", "radius"]
        or positions.shape != (2, *values.shape[1:])
        or position_dims != ["peak", "beat", "branch", "radius"]
    ):
        raise ValueError(
            "area sources must have compatible spatial/peak, beat, branch, "
            "and radius dimensions."
        )

    left_area = np.full(values.shape[1:], np.nan, dtype=np.float32)
    right_area = np.full(values.shape[1:], np.nan, dtype=np.float32)
    for beat_index, branch_index, radius_index in np.ndindex(values.shape[1:]):
        first_peak = positions[0, beat_index, branch_index, radius_index]
        second_peak = positions[1, beat_index, branch_index, radius_index]
        if not np.isfinite(first_peak) or not np.isfinite(second_peak):
            continue
        if first_peak >= second_peak:
            continue
        curve = values[:, beat_index, branch_index, radius_index]
        finite_indexes = np.flatnonzero(np.isfinite(curve))
        if finite_indexes.size < 2:
            continue
        left_edge = float(finite_indexes[0])
        right_edge = float(finite_indexes[-1])
        midpoint = (float(first_peak) + float(second_peak)) / 2.0
        if midpoint <= left_edge or midpoint >= right_edge:
            continue
        left_area[beat_index, branch_index, radius_index] = (
            _finite_trapezoid_area(curve, left_edge, midpoint)
        )
        right_area[beat_index, branch_index, radius_index] = (
            _finite_trapezoid_area(curve, midpoint, right_edge)
        )

    common_attrs = {
        "unit": _integrated_profile_unit(profile_attrs.get("unit", "")),
        "dimDesc": ["beat", "branch", "radius"],
        "integration_axis": profile_dims[0],
        "integration_method": "trapezoidal",
        "outer_boundaries": "first_and_last_finite_profile_samples",
        "shared_boundary": "midpoint_between_XMax_values",
    }
    return (
        DatasetValue(
            data=left_area,
            attrs={**common_attrs, "side": "left"},
            h5_options=profile_h5_options(left_area.shape),
        ),
        DatasetValue(
            data=right_area,
            attrs={**common_attrs, "side": "right"},
            h5_options=profile_h5_options(right_area.shape),
        ),
    )


def _finite_trapezoid_area(
    curve: np.ndarray,
    start: float,
    stop: float,
) -> np.float32:
    """Integrate finite curve runs without bridging across NaN gaps."""

    values = np.asarray(curve, dtype=np.float32)
    finite_indexes = np.flatnonzero(np.isfinite(values))
    if finite_indexes.size < 2 or stop <= start:
        return np.float32(np.nan)

    gap_indexes = np.flatnonzero(np.diff(finite_indexes) > 1) + 1
    finite_runs = np.split(finite_indexes, gap_indexes)
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
    source_unit = str(unit)
    if source_unit == "pixels":
        return "pixels^2"
    if source_unit == "pixels^2":
        return "pixels^3"
    return f"{source_unit}*pixels" if source_unit else "pixels"


def _validate_meaned_peak_metric_sources(
    profile: DatasetValue,
    peak_positions: DatasetValue,
) -> tuple[np.ndarray, np.ndarray, str]:
    values = np.asarray(profile.data, dtype=np.float32)
    profile_dims = list((profile.attrs or {}).get("dimDesc", ()))
    positions = np.asarray(peak_positions.data, dtype=np.float32)
    position_dims = list((peak_positions.attrs or {}).get("dimDesc", ()))
    if (
        values.ndim != 4
        or len(profile_dims) != 4
        or profile_dims[1:] != ["beat", "branch", "radius"]
    ):
        raise ValueError(
            "time-meaned profile dataset must have dimensions "
            "(spatial_sample, beat, branch, radius)."
        )
    if (
        positions.shape != (2, *values.shape[1:])
        or position_dims != ["peak", "beat", "branch", "radius"]
    ):
        raise ValueError(
            "peak positions must have dimensions "
            "(two_peaks, beat, branch, radius)."
        )
    return values, positions, profile_dims[0]


def _mean_peak_metric(metric: DatasetValue) -> DatasetValue:
    """Average one time-meaned peak metric over all branch/radius segments."""

    values = np.asarray(metric.data, dtype=np.float32)
    attrs = dict(metric.attrs or {})
    dim_desc = list(attrs.get("dimDesc", ()))
    if (
        dim_desc != ["peak", "beat", "branch", "radius"]
        or values.ndim != 4
    ):
        raise ValueError(
            "peak metric must have peak, beat, branch, and radius dimensions."
        )
    data = nanmean_float32(values, axis=(2, 3))
    attrs["dimDesc"] = ["peak", "beat"]
    attrs["aggregation"] = "mean_over_valid_branch_radius_values"
    return DatasetValue(
        data=data,
        attrs=attrs,
        h5_options=profile_h5_options(data.shape),
    )


def _first_two_peak_positions(curve: np.ndarray) -> np.ndarray:
    """Find local maxima without treating either side of a NaN gap as a peak."""

    values = np.asarray(curve, dtype=np.float32)
    finite_indexes = np.flatnonzero(np.isfinite(values))
    if finite_indexes.size < 3:
        return np.empty(0, dtype=np.float32)

    gap_indexes = np.flatnonzero(np.diff(finite_indexes) > 1) + 1
    finite_runs = np.split(finite_indexes, gap_indexes)
    positions: list[int] = []
    for run in finite_runs:
        if run.size < 3:
            continue
        local_positions, _ = find_peaks(values[run])
        positions.extend(int(run[index]) for index in local_positions)
        if len(positions) >= 2:
            break
    return np.asarray(positions[:2], dtype=np.float32)


def _temporally_centered_profile_power_dataset(
    profile: DatasetValue,
    temporal_mean: DatasetValue,
) -> DatasetValue:
    """Calculate squared displacement deviations from the temporal mean."""

    values = np.asarray(profile.data, dtype=np.float32)
    mean = np.asarray(temporal_mean.data, dtype=np.float32)
    expected_mean_shape = (values.shape[0], *values.shape[2:])
    if mean.shape != expected_mean_shape:
        raise ValueError(
            "temporal mean shape must match the profile without its time axis."
        )
    centered = values - mean[:, None, ...]
    data = np.square(centered).astype(np.float32, copy=False)
    attrs = dict(profile.attrs or {})
    attrs["unit"] = "pixels^2"
    attrs["formula"] = "(D(t) - mean_t(D(t)))**2"
    attrs["temporal_centering"] = "per_beat_mean"
    return DatasetValue(
        data=data,
        attrs=attrs,
        h5_options=profile_h5_options(data.shape),
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
    profiles_per_beat = interpolate_profiles_per_beat(
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
        h5_options=profile_h5_options(profiles_per_beat.shape[1:]),
    )


def _segment_displacement_metric_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    unit: str,
    metric_attrs: dict[str, object],
) -> DatasetValue:
    """Interpolate a scalar per-segment displacement metric per beat."""

    if profiles.ndim != 3:
        raise ValueError(
            "segment displacement metrics must have shape "
            "(radius, branch, frame)."
        )
    profiles_per_beat = interpolate_profiles_per_beat(
        profiles[..., None],
        cycle_boundary_indexes,
        index_base=index_base,
    )
    data = profiles_per_beat[0]
    return DatasetValue(
        data=data,
        attrs={
            "unit": unit,
            "dimDesc": ["time", "beat", "branch", "radius"],
            "coordinate_system": "rotated_segment_local",
            "displacement_reference": "temporal_mean_image",
            **metric_attrs,
        },
        h5_options=profile_h5_options(data.shape),
    )


def _segment_displacement_magnitude_dataset(
    x_profiles: np.ndarray,
    y_profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> DatasetValue:
    """Interpolate the vector magnitude trace of every vessel segment."""

    if x_profiles.ndim != 3 or x_profiles.shape != y_profiles.shape:
        raise ValueError(
            "X and Y displacement profiles must have matching "
            "(radius, branch, frame) shapes."
        )
    magnitude = np.hypot(x_profiles, y_profiles).astype(
        np.float32,
        copy=False,
    )
    profiles_per_beat = interpolate_profiles_per_beat(
        magnitude[..., None],
        cycle_boundary_indexes,
        index_base=index_base,
    )
    data = profiles_per_beat[0]
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
        },
        h5_options=profile_h5_options(data.shape),
    )


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
    unmasked_values = getattr(
        segments,
        "transverse_fft_profiles_unmasked",
        None,
    )
    masked_values = getattr(
        segments,
        "transverse_fft_profiles_masked",
        None,
    )
    if unmasked_values is None or masked_values is None:
        raise RuntimeError(
            "Velocity FFT profiles were not accumulated during streamed "
            "segment processing."
        )
    unmasked = np.asarray(unmasked_values, dtype=np.float32)
    masked = np.asarray(masked_values, dtype=np.float32)
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
            h5_options=profile_h5_options(unmasked.shape),
        ),
        masked_path: DatasetValue(
            data=masked,
            attrs={
                **shared_attrs,
                "mask_applied": True,
                "mask_dilation_iterations": int(mask_dilation_pixels),
                "mask_dilation_axis": "transverse_x",
            },
            h5_options=profile_h5_options(masked.shape),
        ),
    }
