"""Numerical contracts for blood-volume-rate calculations."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from calculations.blood_volume_rate import (
    circular_lumen_profile_flow,
    mask_derived_lumen_geometry,
    masked_edges_flow,
    total_masked_edges_flow,
)
from calculations.topology import AnnulusGeometry
from input_output.profile_datasets import _profile_dataset
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue
from pipelines.blood_volume_rate.outputs import (
    pack_gradient_edge_outputs,
    pack_mask_derived_outputs,
)


def _profile(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=np.float32).reshape((-1, 1, 1, 1, 1))


def _edge(value: float) -> np.ndarray:
    return np.full((1, 1, 1, 1), value, dtype=np.float32)


def test_constant_velocity_integrates_to_velocity_times_circle_area() -> None:
    velocity = 3.25
    left = 2.25
    right = 7.75
    pixel_size_mm = 0.02

    actual = circular_lumen_profile_flow(
        _profile(np.full(11, velocity)),
        _edge(left),
        _edge(right),
        profile_pixel_size_mm=pixel_size_mm,
    )

    radius_mm = (right - left) * 0.5 * pixel_size_mm
    expected = velocity * np.pi * radius_mm**2
    np.testing.assert_allclose(actual, expected, rtol=2e-6)


def test_linear_fractional_profile_matches_high_accuracy_integration() -> None:
    coordinates = np.arange(12, dtype=np.float64)
    samples = 1.75 - 0.32 * coordinates
    left = 1.4
    right = 9.65
    pixel_size_mm = 0.017

    actual = circular_lumen_profile_flow(
        _profile(samples),
        _edge(left),
        _edge(right),
        profile_pixel_size_mm=pixel_size_mm,
    )[0, 0, 0, 0]

    x = np.linspace(left, right, 1_000_001)
    center = 0.5 * (left + right)
    radius = 0.5 * (right - left)
    chord = 2.0 * np.sqrt(np.maximum(radius**2 - (x - center) ** 2, 0.0))
    expected = np.trapezoid((1.75 - 0.32 * x) * chord, x) * pixel_size_mm**2
    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=1e-9)


def test_invalid_edges_and_profiles_without_finite_intervals_are_nan() -> None:
    profiles = np.broadcast_to(
        _profile(np.arange(6)),
        (6, 1, 1, 1, 4),
    ).copy()
    profiles[:, :, :, :, 3] = np.nan
    left = np.asarray([[[[4.0, 2.0, 2.5, 1.0]]]], dtype=np.float32)
    right = np.asarray([[[[2.0, 2.0, np.nan, 4.0]]]], dtype=np.float32)

    actual = circular_lumen_profile_flow(
        profiles,
        left,
        right,
        profile_pixel_size_mm=0.01,
    )

    assert np.all(np.isnan(actual))


def test_mask_geometry_and_signed_flow_keep_established_model() -> None:
    labels = np.asarray(
        [
            [0, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 1, 1, 0],
        ],
        dtype=np.int32,
    )
    topology = SimpleNamespace(
        labels=labels,
        branch_ids=np.asarray([1], dtype=np.int32),
        annulus_masks=np.ones((1, 4, 4), dtype=bool),
        ring_settings=AnnulusGeometry(0.0, 0.5, 0.5, 1, 0.5),
    )
    diameters, areas, widths = mask_derived_lumen_geometry(
        (topology,),
        pixel_size_mm=0.1,
    )
    expected_diameter = 8 * 0.1 / widths[0]
    np.testing.assert_array_equal(areas[0], [[8]])
    np.testing.assert_allclose(diameters[0], [[expected_diameter]])

    velocity = np.full((8, 1, 1, 1), -2.0, dtype=np.float32)
    rate = masked_edges_flow(velocity, diameters[0])
    expected_rate = -2.0 * np.pi / 4.0 * expected_diameter**2
    np.testing.assert_allclose(rate, expected_rate)
    np.testing.assert_allclose(total_masked_edges_flow(rate), expected_rate)


def test_total_masked_edges_flow_uses_centered_periodic_nine_point_window() -> None:
    rate = np.zeros((10, 1, 1, 1), dtype=np.float32)
    rate[0] = 9.0

    expected = np.ones((10, 1), dtype=np.float32)
    expected[5, 0] = 0.0
    np.testing.assert_array_equal(total_masked_edges_flow(rate), expected)


def test_output_packers_keep_paths_units_and_valid_provenance() -> None:
    prepared = object()
    velocity_topology = SimpleNamespace(
        valid_segments=np.ones((1, 1), dtype=bool),
        segment_center_xy=np.asarray([[[2.0, 3.0]]], dtype=np.float32),
        profile_rotation_degrees=np.asarray([[17.0]], dtype=np.float32),
        prepared_topology=prepared,
    )
    profiles = np.full((1, 1, 3, 6), -2.0, dtype=np.float32)
    velocity_segments = SimpleNamespace(
        labels=np.asarray([[1]], dtype=np.int32),
        branch_ids=np.asarray([1], dtype=np.int32),
        topology=velocity_topology,
        transverse_velocity_profiles_masked=profiles,
        profile_pixel_size_mm=0.02,
    )
    gradient_segments = SimpleNamespace(
        labels=velocity_segments.labels,
        branch_ids=velocity_segments.branch_ids,
        topology=SimpleNamespace(
            segment_centers_xy=np.asarray([[[2.0, 3.0]]], dtype=np.float32),
            profile_rotation_degrees=np.asarray([[17.0]], dtype=np.float32),
            prepared_topology=prepared,
        ),
    )
    cycle_boundaries = np.asarray([0, 2], dtype=np.int32)
    profile_dataset = _profile_dataset(
        profiles,
        cycle_boundaries,
        index_base=0,
        valid_segments=np.ones((1, 1), dtype=bool),
    )
    edge_shape = profile_dataset.data.shape[1:]
    edges = {}
    for vessel in ("Artery", "Vein"):
        root = (
            f"Processing/SpatialGradientMetrics/{vessel}/"
            "Transverse/Masked/tbkr"
        )
        edges[f"{root}/left_edge_index"] = DatasetValue(
            np.full(edge_shape, 0.5, dtype=np.float32)
        )
        edges[f"{root}/right_edge_index"] = DatasetValue(
            np.full(edge_shape, 4.5, dtype=np.float32)
        )
    gradient_outputs = pack_gradient_edge_outputs(
        velocity_segments,
        velocity_segments,
        SimpleNamespace(
            artery_segments=gradient_segments,
            vein_segments=gradient_segments,
            outputs=edges,
        ),
        cycle_boundaries,
        index_base=0,
    )

    topology = SimpleNamespace(
        labels=np.ones((4, 4), dtype=np.int32),
        branch_ids=np.asarray([1], dtype=np.int32),
        annulus_masks=np.ones((1, 4, 4), dtype=bool),
        ring_settings=AnnulusGeometry(0.0, 0.5, 0.5, 1, 0.5),
    )
    schema = EyeFlowOutputPaths.active()
    safe_velocity = np.full((8, 1, 1, 1), -2.0, dtype=np.float32)
    mask_outputs = pack_mask_derived_outputs(
        {"artery": topology, "vein": topology},
        {
            schema.artery_per_beat_safe.velocity_signal: safe_velocity,
            schema.vein_per_beat_safe.velocity_signal: safe_velocity,
        },
        pixel_size_mm=0.1,
    )
    outputs = {**gradient_outputs, **mask_outputs}

    for vessel_paths in (
        schema.blood_volume_rate.artery,
        schema.blood_volume_rate.vein,
    ):
        for path in (
            vessel_paths.dynamic_edges,
            vessel_paths.static_edges,
            vessel_paths.masked_edges,
            vessel_paths.total_masked_edges,
        ):
            assert path in outputs
            assert outputs[path].attrs["unit"] == "mm^3/s"
        assert not gradient_outputs[vessel_paths.dynamic_edges].attrs[
            "source_velocity"
        ].startswith("/")
        assert (
            mask_outputs[vessel_paths.total_masked_edges].attrs["source"]
            == f"/{vessel_paths.masked_edges}"
        )

