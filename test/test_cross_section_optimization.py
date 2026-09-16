"""Bounded topology-based cross-section regression checks."""
import importlib
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pytest

from calculations.topology import (
    SegmentRingSettings, dilate_segment_masks, prepare_segments, prepare_topology,
)
from calculations.blood_flow_velocity.cross_section.reusable_cross_section_signals import (
    fit_cross_section_plan, project_cross_section_cube,
)

cs = importlib.import_module(
    "calculations.blood_flow_velocity.cross_section.generate_cross_section_signals"
)


@pytest.fixture(autouse=True)
def cpu(monkeypatch):
    from calculations.compute_backend import optional_cupy_backend
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cpu")
    optional_cupy_backend.cache_clear()
    yield
    optional_cupy_backend.cache_clear()


@pytest.fixture
def geometry():
    vessel = np.zeros((61, 61), bool)
    vessel[27:34, 5:56] = True
    disc = np.zeros_like(vessel)
    disc[27:34, 27:34] = True
    return vessel, disc, SegmentRingSettings(.1, .7, .25, 2)


@pytest.mark.parametrize("field,value", [
    ("pixel_size_mm", 0), ("pixel_size_mm", np.nan),
    ("working_memory_mb", -1), ("submask_size_percentile_kept", np.nan),
])
def test_invalid_settings_are_rejected(field, value):
    with pytest.raises(ValueError):
        replace(cs.CrossSectionSignalSettings(.01), **{field: value})


def test_missing_frames_remain_missing_and_zero_flow_remains_zero(geometry):
    vessel, disc, rings = geometry
    cube = np.zeros((2, *vessel.shape), np.float32)
    cube[0] = np.nan
    result = cs.generate_cross_section_signals(
        cube, vessel, (30, 30), rings, cs.CrossSectionSignalSettings(.01),
        retain_velocity_maps=False,
    )
    valid = result.topology.valid_segments
    assert valid.any()
    assert np.isnan(result.velocity[valid, 0]).all()
    np.testing.assert_array_equal(result.velocity[valid, 1], 0)
    assert result.velocity_maps_per_segment is None


def test_missing_intersections_allocate_no_movie_payload():
    buffers = cs._CrossSectionBuffers.allocate(
        frame_count=1000, ring_count=10, branch_count=20,
        velocity_map_segment_indexes=np.empty((0, 2), dtype=np.int32),
        retain_velocity_maps=True,
    )
    assert buffers.velocity_maps_per_segment.nbytes == 0
    assert buffers.velocity_maps_per_segment.shape == (0, 1000, 181, 181)


def test_spatial_gradient_mask_expands_five_pixels_horizontally():
    mask = np.zeros((31, 31), bool)
    mask[15, 15] = True
    dilated = dilate_segment_masks(mask, iterations=5, horizontal_only=True)
    expected = np.zeros_like(mask)
    expected[15, 10:21] = True
    np.testing.assert_array_equal(dilated, expected)


def test_windows_are_extracted_only_when_iteration_starts(geometry):
    vessel, disc, rings = geometry
    prepared = prepare_topology(vessel, disc, rings)
    workflow = importlib.import_module("calculations.topology.workflow")
    cube = np.ones((2, *vessel.shape), np.float32)
    with patch.object(workflow, "extract_segment", wraps=workflow.extract_segment) as extract:
        segments = prepare_segments(cube, prepared)
        extract.assert_not_called()
        next(segments)
        extract.assert_called_once()


def test_worker_count_obeys_memory_budget():
    with patch.object(cs, "cap_parallel_jobs", side_effect=lambda n, **kw: n):
        assert cs._cross_section_worker_count(
            100, frame_count=1000, working_memory_mb=1
        ) == 1


def test_cuda_profile_measurement_matches_cpu(geometry, monkeypatch):
    from calculations.compute_backend import optional_cupy_backend
    vessel, disc, rings = geometry
    cube = np.ones((3, *vessel.shape), np.float32)
    cube[0] = np.nan
    cpu_result = cs.generate_cross_section_signals(
        cube, vessel, (30, 30), rings, cs.CrossSectionSignalSettings(.01),
        retain_velocity_maps=False,
    )
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
    optional_cupy_backend.cache_clear()
    if optional_cupy_backend() is None:
        pytest.skip("CuPy/CUDA unavailable")
    gpu_result = cs.generate_cross_section_signals(
        cube, vessel, (30, 30), rings, cs.CrossSectionSignalSettings(.01),
        retain_velocity_maps=False,
    )
    for field in (
        "velocity", "safe_velocity", "velocity_profiles",
        "transverse_velocity_profiles_masked",
        "longitudinal_velocity_profiles_masked",
    ):
        np.testing.assert_allclose(
            getattr(gpu_result, field), getattr(cpu_result, field),
            rtol=1e-5, atol=1e-5, equal_nan=True,
        )
    assert gpu_result.velocity_maps_per_segment is None


def test_velocity_and_gradient_share_topology_and_segment_axes(geometry):
    from pipelines.waveform_velocity_core.segments import analyze_velocity_segments
    from pipelines.waveform_velocity.runner import _validate_profile_segment_alignment
    vessel, disc, rings = geometry
    cube = np.ones((2, *vessel.shape), np.float32)
    masks = {"artery": vessel, "vein": np.zeros_like(vessel)}
    cache = {}
    common = dict(optic_disc_mask=disc, source_id="registered", topology_cache=cache)
    velocity = analyze_velocity_segments(
        cube, masks, (30, 30), rings, cs.CrossSectionSignalSettings(.01), **common,
    )
    from calculations.topology import prepare_topologies
    from pipelines.spatial_gradient_moment0.profiles import _project_spatial_gradient_segments
    topologies = prepare_topologies(
        masks, disc, rings, source_id="registered", cache=cache, optic_disc_center=(30, 30),
    )
    gradient = {
        name: _project_spatial_gradient_segments(cube * 3, prepared)
        for name, prepared in topologies.items()
    }
    for name in masks:
        _validate_profile_segment_alignment(name, velocity[name], gradient[name])
        if velocity[name].branch_ids.size:
            assert velocity[name].topology.prepared_topology is gradient[name].prepared_topology
        assert not hasattr(gradient[name], "velocity_maps_per_segment")


@pytest.mark.parametrize("mode", ["reference", "per_cube"])
def test_registered_cubes_reuse_one_prepared_topology(geometry, mode):
    vessel, disc, rings = geometry
    cube = np.ones((2, *vessel.shape), np.float32)
    plan, first = fit_cross_section_plan(
        cube, vessel, (30, 30), rings, cs.CrossSectionSignalSettings(.01)
    )
    with patch(
        "calculations.blood_flow_velocity.cross_section.reusable_cross_section_signals.prepare_topologies",
        side_effect=AssertionError("topology must not be refitted"),
    ):
        second = project_cross_section_cube(cube * 2, plan, limits_mode=mode)
    assert second.topology.prepared_topology is plan.prepared_topology
    np.testing.assert_allclose(second.velocity, first.velocity * 2, equal_nan=True)
    np.testing.assert_array_equal(second.branch_ids, first.branch_ids)
    with pytest.raises(ValueError, match="spatial"):
        project_cross_section_cube(np.ones((2, 5, 5)), plan)
