"""Bounded topology-based cross-section regression checks."""
import importlib
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pytest

from calculations.topology import (
    AnnulusGeometry, OpticDisc, dilate_segment_masks, prepare_segments,
    prepare_topology,
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
    return vessel, OpticDisc(disc, (30.0, 30.0), None, None), AnnulusGeometry(.1, .7, .25, 2)


@pytest.mark.parametrize("field,value", [
    ("pixel_size_mm", 0), ("pixel_size_mm", np.nan),
    ("working_memory_mb", -1), ("submask_size_percentile_kept", np.nan),
])
def test_invalid_settings_are_rejected(field, value):
    with pytest.raises(ValueError):
        replace(cs.CrossSectionSignalSettings(.01), **{field: value})


def test_missing_frames_remain_missing_and_zero_flow_remains_zero(geometry):
    vessel, optic_disc, rings = geometry
    cube = np.zeros((2, *vessel.shape), np.float32)
    cube[0] = np.nan
    result = cs.generate_cross_section_signals(
        cube, vessel, optic_disc, rings, cs.CrossSectionSignalSettings(.01),
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
    vessel, optic_disc, rings = geometry
    prepared = prepare_topology(vessel, optic_disc, rings)
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
    vessel, optic_disc, rings = geometry
    cube = np.ones((3, *vessel.shape), np.float32)
    cube[0] = np.nan
    cpu_result = cs.generate_cross_section_signals(
        cube, vessel, optic_disc, rings, cs.CrossSectionSignalSettings(.01),
        retain_velocity_maps=False,
    )
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
    optional_cupy_backend.cache_clear()
    if optional_cupy_backend() is None:
        pytest.skip("CuPy/CUDA unavailable")
    gpu_result = cs.generate_cross_section_signals(
        cube, vessel, optic_disc, rings, cs.CrossSectionSignalSettings(.01),
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
    from calculations.segment_profiles import analyze_segment_profiles
    from pipelines.spatial_gradient_moment0.runner import (
        _validate_profile_segment_alignment,
    )
    vessel, optic_disc, rings = geometry
    cube = np.ones((2, *vessel.shape), np.float32)
    masks = {"artery": vessel, "vein": np.zeros_like(vessel)}
    cache = {}
    common = dict(source_id="registered", topology_cache=cache)
    velocity = analyze_segment_profiles(
        cube, masks, optic_disc, rings, cs.CrossSectionSignalSettings(.01), **common,
    )
    from pipelines.spatial_gradient_moment0.profiles import _spatial_gradient_chain
    topologies = {
        name: result.topology.prepared_topology
        for name, result in velocity.items()
    }
    gradient = analyze_segment_profiles(
        cube * 3,
        masks,
        optic_disc,
        rings,
        cs.CrossSectionSignalSettings(.01),
        prepared_topologies=topologies,
        transform_mode="staged",
        post_interpolation=_spatial_gradient_chain,
        temporal_halo=6,
        scratch_array_count=7,
    )
    for name in masks:
        assert not hasattr(velocity[name], "velocity")
        assert not hasattr(velocity[name], "displacements")
        _validate_profile_segment_alignment(name, velocity[name], gradient[name])
        if velocity[name].branch_ids.size:
            assert (
                velocity[name].topology.prepared_topology
                is gradient[name].topology.prepared_topology
            )
        assert gradient[name].segment_maps is None


def test_velocity_adapter_preserves_legacy_profile_values(geometry):
    from pipelines.waveform_velocity_core.segments import (
        analyze_velocity_segment_profiles,
    )

    vessel, optic_disc, rings = geometry
    y, x = np.indices(vessel.shape, dtype=np.float32)
    cube = np.stack((x + y, 2 * x - y), axis=0)
    settings = cs.CrossSectionSignalSettings(.01)
    legacy = cs.generate_cross_section_signals(
        cube,
        vessel,
        optic_disc,
        rings,
        settings,
        retain_velocity_maps=False,
    )
    adapted = analyze_velocity_segment_profiles(
        cube,
        {"vessel": vessel},
        optic_disc,
        rings,
        settings,
        retain_velocity_maps=False,
        transverse_mask_dilation_pixels=0,
    )["vessel"]

    for field in (
        "velocity",
        "safe_velocity",
        "velocity_profiles",
        "transverse_velocity_profiles_masked",
        "longitudinal_velocity_profiles_unmasked",
        "longitudinal_velocity_profiles_masked",
        "rotated_mean_images",
        "rotated_mean_images_masked",
    ):
        np.testing.assert_allclose(
            getattr(adapted, field),
            getattr(legacy, field),
            rtol=1e-6,
            atol=1e-6,
            equal_nan=True,
        )
    np.testing.assert_array_equal(adapted.segment_masks, legacy.segment_masks)
    np.testing.assert_array_equal(
        adapted.topology.valid_segments,
        legacy.topology.valid_segments,
    )


def test_legacy_profile_dilation_does_not_change_segment_velocity(geometry):
    from pipelines.waveform_velocity_core.segments import (
        analyze_velocity_segment_profiles,
    )

    vessel, optic_disc, rings = geometry
    cube = np.broadcast_to(
        np.where(vessel, np.float32(10.0), np.float32(1.0)),
        (3, *vessel.shape),
    ).copy()
    results = analyze_velocity_segment_profiles(
        cube,
        {"artery": vessel, "vein": vessel},
        optic_disc,
        rings,
        cs.CrossSectionSignalSettings(.01),
    )
    artery = results["artery"]
    vein = results["vein"]
    np.testing.assert_allclose(artery.velocity, vein.velocity, equal_nan=True)
    np.testing.assert_allclose(
        artery.rotated_mean_images_masked,
        vein.rotated_mean_images_masked,
        equal_nan=True,
    )
    assert np.count_nonzero(
        np.isfinite(artery.transverse_velocity_profiles_masked)
    ) > np.count_nonzero(
        np.isfinite(vein.transverse_velocity_profiles_masked)
    )


@pytest.mark.parametrize("mode", ["reference", "per_cube"])
def test_registered_cubes_reuse_one_prepared_topology(geometry, mode):
    vessel, optic_disc, rings = geometry
    cube = np.ones((2, *vessel.shape), np.float32)
    plan, first = fit_cross_section_plan(
        cube, vessel, optic_disc, rings, cs.CrossSectionSignalSettings(.01)
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
