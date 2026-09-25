"""Serialize blood-volume-rate calculations to the EyeFlow schema."""

from __future__ import annotations

import numpy as np

from calculations.blood_volume_rate import (
    TOTAL_MASKED_EDGES_WINDOW_SIZE,
    TOTAL_MASKED_EDGES_WINDOW_STRIDE,
    circular_lumen_profile_flow,
    mask_derived_lumen_geometry,
    masked_edges_flow,
    total_masked_edges_flow,
)
from calculations.math import nanmean_float32
from input_output.profile_datasets import _profile_dataset, _profile_h5_options
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue


def pack_gradient_edge_outputs(
    artery_velocity_segments,
    vein_velocity_segments,
    gradient_products,
    cycle_boundary_indexes,
    *,
    index_base: int,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, DatasetValue]:
    """Calculate dynamic- and static-edge flow for both vessel classes."""

    schema = _resolve_output_paths(output_paths)
    outputs: dict[str, DatasetValue] = {}
    vessels = (
        (
            "Artery",
            artery_velocity_segments,
            gradient_products.artery_segments,
            schema.blood_volume_rate.artery,
        ),
        (
            "Vein",
            vein_velocity_segments,
            gradient_products.vein_segments,
            schema.blood_volume_rate.vein,
        ),
    )
    for vessel_name, velocity, gradient, paths in vessels:
        _validate_profile_segment_alignment(vessel_name, velocity, gradient)
        profile = _profile_dataset(
            np.asarray(velocity.transverse_velocity_profiles_masked, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
            valid_segments=np.asarray(velocity.topology.valid_segments, dtype=bool),
        )
        metrics_root = (
            f"Processing/SpatialGradientMetrics/{vessel_name}/"
            "Transverse/Masked/tbkr"
        )
        left_path = f"{metrics_root}/left_edge_index"
        right_path = f"{metrics_root}/right_edge_index"
        left_value = gradient_products.outputs[left_path]
        right_value = gradient_products.outputs[right_path]
        pixel_size_mm = float(velocity.profile_pixel_size_mm)
        outputs[paths.dynamic_edges] = _gradient_edge_dataset(
            profile,
            left_value,
            right_value,
            profile_pixel_size_mm=pixel_size_mm,
            left_edge_path=left_path,
            right_edge_path=right_path,
            static_edges=False,
        )
        outputs[paths.static_edges] = _gradient_edge_dataset(
            profile,
            left_value,
            right_value,
            profile_pixel_size_mm=pixel_size_mm,
            left_edge_path=left_path,
            right_edge_path=right_path,
            static_edges=True,
        )
    return outputs


def _gradient_edge_dataset(
    profile: DatasetValue,
    left_edge: DatasetValue,
    right_edge: DatasetValue,
    *,
    profile_pixel_size_mm: float,
    left_edge_path: str,
    right_edge_path: str,
    static_edges: bool,
) -> DatasetValue:
    left = np.asarray(left_edge.data, dtype=np.float32)
    right = np.asarray(right_edge.data, dtype=np.float32)
    if static_edges:
        left = np.broadcast_to(nanmean_float32(left, axis=(0, 1)), left.shape)
        right = np.broadcast_to(nanmean_float32(right, axis=(0, 1)), right.shape)
    rate = circular_lumen_profile_flow(
        profile.data,
        left,
        right,
        profile_pixel_size_mm=profile_pixel_size_mm,
    )
    return DatasetValue(
        rate,
        {
            "unit": "mm^3/s",
            "dimDesc": ["time", "beat", "branch", "radius"],
            "definition": (
                "analytic piecewise-linear integration of transverse velocity "
                "weighted by the circular chord implied by the lumen edges"
            ),
            "source_velocity": "waveform_velocity_core.masked_transverse_profile",
            "source_left_edge_index": f"/{left_edge_path.lstrip('/')}",
            "source_right_edge_index": f"/{right_edge_path.lstrip('/')}",
            "cross_section_model": "circular_chord_from_gradient_edges",
            "velocity_interpolation": "piecewise_linear",
            "integration_method": "analytic_linear_velocity_times_circular_chord",
            "profile_pixel_size_mm": np.float32(profile_pixel_size_mm),
            "edge_temporal_reduction": (
                "mean_over_time_and_beats" if static_edges else "none"
            ),
        },
        h5_options=_profile_h5_options(rate.shape),
    )


def pack_mask_derived_outputs(
    prepared_topologies,
    velocity_per_beat_outputs: dict[str, object],
    *,
    pixel_size_mm: float,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, DatasetValue]:
    """Calculate masked-edge and total masked-edge flow for both vessels."""

    schema = _resolve_output_paths(output_paths)
    diameters, _, radial_widths = mask_derived_lumen_geometry(
        (prepared_topologies["artery"], prepared_topologies["vein"]),
        pixel_size_mm=pixel_size_mm,
    )
    vessel_sources = (
        (
            "Artery",
            schema.artery_per_beat_safe.velocity_signal,
            schema.blood_volume_rate.artery,
            diameters[0],
        ),
        (
            "Vein",
            schema.vein_per_beat_safe.velocity_signal,
            schema.blood_volume_rate.vein,
            diameters[1],
        ),
    )
    outputs: dict[str, DatasetValue] = {}
    for vessel_name, velocity_path, paths, diameter_mm in vessel_sources:
        if velocity_path is None or velocity_path not in velocity_per_beat_outputs:
            raise KeyError(
                f"Required safe per-beat velocity is unavailable for {vessel_name}."
            )
        velocity = _metric_data(velocity_per_beat_outputs[velocity_path])
        rate = masked_edges_flow(velocity, diameter_mm)
        masked = DatasetValue(
            rate,
            {
                "unit": "mm^3/s",
                "dimDesc": ["time", "beat", "branch", "radius"],
                "definition": (
                    "safe per-beat segment velocity multiplied by an equivalent "
                    "circular lumen area derived from native vessel-mask pixels"
                ),
                "source_velocity": "waveform_velocity_core.safe_per_beat_segment_velocity",
                "diameter_model": "masked_pixel_count_over_radial_width",
                "diameter_model_assumption": "locally_radial_vessel",
                "cross_section_model": "circular_pi_diameter_squared_over_4",
                "annulus_geometry": "native_pixel_center_section_mask",
                "annulus_edge_handling": "outer_radius_clipped_to_configured_limit",
                "native_pixel_size_mm": np.float32(pixel_size_mm),
                "radial_width_pixels": radial_widths,
            },
            h5_options=_profile_h5_options(rate.shape),
        )
        outputs[paths.masked_edges] = masked
        total = total_masked_edges_flow(rate)
        outputs[paths.total_masked_edges] = DatasetValue(
            total,
            {
                "unit": "mm^3/s",
                "dimDesc": ["time", "beat"],
                "definition": (
                    "median over radius of the sum over branches of masked-edge "
                    "blood-volume rate after a circular sliding average over tau"
                ),
                "source": f"/{paths.masked_edges.lstrip('/')}",
                "aggregation": "median_over_radius_of_sum_over_branches",
                "temporal_filter": "circular_sliding_average_over_tau",
                "temporal_window_size": np.int32(TOTAL_MASKED_EDGES_WINDOW_SIZE),
                "temporal_window_stride": np.int32(TOTAL_MASKED_EDGES_WINDOW_STRIDE),
                "temporal_boundary_mode": "circular",
                "temporal_window_alignment": "forward",
                "temporal_nan_policy": "propagate",
                "branch_reduction": "sum_over_finite_values",
                "radius_reduction": "median_over_finite_values",
            },
            h5_options=_profile_h5_options(total.shape),
        )
    return outputs


def _validate_profile_segment_alignment(vessel_name, velocity, gradient) -> None:
    for field in ("labels", "branch_ids"):
        if not np.array_equal(
            np.asarray(getattr(velocity, field)),
            np.asarray(getattr(gradient, field)),
        ):
            raise RuntimeError(f"{vessel_name} segment {field} do not match.")
    velocity_topology = velocity.topology
    gradient_topology = gradient.topology
    velocity_centers = np.transpose(
        np.asarray(velocity_topology.segment_center_xy),
        (1, 0, 2),
    )
    if not np.allclose(
        velocity_centers,
        np.asarray(gradient_topology.segment_centers_xy),
        equal_nan=True,
    ):
        raise RuntimeError(f"{vessel_name} segment centers do not match.")
    if not np.allclose(
        np.asarray(velocity_topology.profile_rotation_degrees),
        np.asarray(gradient_topology.profile_rotation_degrees),
        equal_nan=True,
    ):
        raise RuntimeError(f"{vessel_name} segment rotations do not match.")
    if velocity_topology.prepared_topology is not gradient_topology.prepared_topology:
        raise RuntimeError(f"{vessel_name} analyses did not share prepared topology.")


def _metric_data(value) -> np.ndarray:
    if isinstance(value, DatasetValue):
        value = value.data
    elif isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        value = value[0]
    return np.asarray(value)


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


__all__ = ["pack_gradient_edge_outputs", "pack_mask_derived_outputs"]
