"""Contracts for the optional displacement profile and map families."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import h5py

from calculations.math import interpft_real
from calculations.vessel_segments.profiles.per_beat import (
    interpolate_segment_maps_per_beat,
)
from input_output.schema import EyeFlowOutputPaths
from input_output.writers.h5_output import PipelineH5Output
from pipeline_engine import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.displacement_map.runner import DisplacementMapPipelineConfig
from pipelines.displacement_map.graph import (
    displacement_graph,
    displacement_profiles_graph,
)
from pipelines.displacement_map import timing as displacement_timing
from pipelines.displacement_map.profile_outputs import (
    pack_displacement_profile_outputs,
    pack_displacement_segment_map_outputs,
)
from pipelines.displacement_map import segments as displacement_segments
from pipelines.displacement_map.segments import DisplacementSegmentResult


def _result(*, retain_maps: bool = True) -> DisplacementSegmentResult:
    signal_shape = (1, 1, 5)
    curve = np.asarray([0.0, 2.0, 0.0, 1.0, 0.0], dtype=np.float32)
    profiles = np.broadcast_to(curve, (*signal_shape, curve.size)).copy()
    scalar = np.ones(signal_shape, dtype=np.float32)
    maps = (
        np.ones((*signal_shape, curve.size, curve.size, 2), dtype=np.float32)
        if retain_maps
        else None
    )
    return DisplacementSegmentResult(
        maps=maps,
        transverse_profiles_unmasked=profiles,
        transverse_profiles_masked=profiles,
        longitudinal_profiles_unmasked=profiles,
        longitudinal_profiles_masked=profiles,
        x_sum_profile=scalar,
        y_sum_profile=scalar,
        radial_movement_amplitude=scalar,
        radial_asymmetry_index=scalar,
    )


def test_all_restored_profile_and_metric_paths_use_the_final_names() -> None:
    schema = EyeFlowOutputPaths.active()
    timing_attrs = {"cycle_boundary_workflow": "doppler_moments"}
    outputs = pack_displacement_profile_outputs(
        {"artery": _result()},
        "level_set_motion",
        np.asarray([0, 4]),
        timing_attrs=timing_attrs,
    )
    profile_root = f"{schema.displacement_profiles_root}/LevelSetMotion/Artery"
    metric_root = f"{schema.displacement_metrics_root}/LevelSetMotion/Artery"
    expected = {
        f"{profile_root}/{name}"
        for name in (
            "X_sum_displacement_profile",
            "Y_sum_displacement_profile",
            "DisplacementMagnitude",
            "Cross_sectional_radial_movement_amplitude_profile",
            "Cross_sectional_radial_asymmetry_index_profile",
        )
    }
    for direction in ("Longitudinal", "Transverse"):
        for mask in ("Masked", "Unmasked"):
            root = f"{profile_root}/{direction}/{mask}"
            expected.update(
                {
                    f"{root}/tbkr/Profile",
                    f"{root}/bkr/Profile",
                    f"{root}/b/Profile",
                    f"{root}/tbkr/SquaredDeviation",
                    f"{root}/bkr/SquaredDeviation",
                }
            )
    expected.update(
        {
            f"{metric_root}/Transverse/Position/{name}"
            for name in (
                "XMax",
                "YMax",
                "YMaxTemporalVariance",
                "XMaxMeaned",
                "YMaxMeaned",
                "YMaxMeanedTemporalVariance",
            )
        }
    )
    expected.update(
        {
            f"{metric_root}/Transverse/Area/{name}"
            for name in (
                "Left",
                "Right",
                "LeftTemporalVariance",
                "RightTemporalVariance",
            )
        }
    )

    assert set(outputs) == expected
    assert all(not path.endswith("/value") for path in outputs)
    profile = outputs[f"{profile_root}/Transverse/Masked/tbkr/Profile"]
    assert profile.data.shape == (5, 4, 1, 1, 1)
    assert profile.attrs["dimDesc"] == ["x", "time", "beat", "branch", "radius"]
    assert profile.attrs["cycle_boundary_workflow"] == "doppler_moments"


def test_vector_segment_map_interpolation_preserves_components_and_axis_order() -> None:
    maps = np.arange(1 * 1 * 5 * 2 * 3 * 2, dtype=np.float32).reshape(
        1, 1, 5, 2, 3, 2
    )

    actual = interpolate_segment_maps_per_beat(maps, np.asarray([0, 4]))
    expected = interpft_real(maps[0, 0, :, 1, 2, 1], 5)[:-1]

    assert actual.shape == (3, 2, 4, 1, 1, 1, 2)
    np.testing.assert_allclose(actual[2, 1, :, 0, 0, 0, 1], expected, rtol=1e-6)


def test_segment_maps_are_packed_only_by_the_separate_map_packer() -> None:
    schema = EyeFlowOutputPaths.active()
    assert schema.displacement_map_root == "ProcessingWidth/Displacement/Map"
    assert schema.displacement_profiles_root == "ProcessingWidth/Displacement/Profiles"
    assert schema.displacement_metrics_root == "ProcessingWidth/Displacement/Metrics"
    result = _result()
    profile_outputs = pack_displacement_profile_outputs(
        {"vein": result},
        "level_set_motion",
        np.asarray([0, 4]),
    )
    map_outputs = pack_displacement_segment_map_outputs(
        {"vein": result},
        "level_set_motion",
        np.asarray([0, 4]),
    )
    path = f"{schema.displacement_map_root}/LevelSetMotion/Vein/PerSegment"

    assert path not in profile_outputs
    assert set(map_outputs) == {path}
    dataset = map_outputs[path]
    assert dataset.data.shape == (5, 5, 4, 1, 1, 1, 2)
    assert dataset.attrs["components"] == ["local_x", "local_y"]
    assert dataset.h5_options["sparse_nan_chunks"] is True


def test_displacement_registration_and_analysis_are_separate_pipelines() -> None:
    load_pipeline_catalog()
    options = {
        option.name: option
        for option in PIPELINE_REGISTRY["displacement_profiles"].options
    }

    base_descriptor = PIPELINE_REGISTRY["displacement_map"]
    profile_descriptor = PIPELINE_REGISTRY["displacement_profiles"]
    assert options["segment_maps"].default_enabled is False
    assert base_descriptor.dag_requires == ()
    assert base_descriptor.execution_mode == "process"
    assert base_descriptor.native_parallelism is True
    assert profile_descriptor.execution_mode == "process"
    assert profile_descriptor.run_once is True
    assert profile_descriptor.native_parallelism is False
    assert profile_descriptor.dag_requires == (
        "prepared_topology",
        "velocity_workflows",
    )

    assert profile_descriptor.product_requires[0].product == "registered_fields"
    assert profile_descriptor.product_requires[0].producer_pipeline == "displacement_map"

    base_graph = displacement_graph(DisplacementMapPipelineConfig())
    assert [node.id for node in base_graph.stages] == [
        "prepare",
        "artery",
        "vein",
        "finalize",
    ]

    profile_graph = displacement_profiles_graph()
    assert [node.id for node in profile_graph.stages] == ["analysis", "profiles"]
    assert [product.id for product in profile_graph.stages[1].products] == [
        "profiles"
    ]

    map_graph = displacement_profiles_graph(retain_segment_maps=True)
    assert [product.id for stage in map_graph.stages[1:] for product in stage.products] == [
        "profiles",
        "segment_maps",
    ]


def test_segment_results_use_the_run_scratch_allocator(monkeypatch) -> None:
    topology = SimpleNamespace(
        segment_shape=(1, 1),
        rotated_masks=np.ones((1, 1, 3, 3), dtype=bool),
    )
    allocated = []

    def allocate(shape):
        allocated.append(shape)
        return np.empty(shape, dtype=np.float32)

    monkeypatch.setattr(
        displacement_segments,
        "prepare_segment_chunks",
        lambda *_args, **_kwargs: iter(()),
    )
    monkeypatch.setattr(
        displacement_segments,
        "dilate_segment_masks",
        lambda masks, **_kwargs: masks,
    )

    result = displacement_segments.analyze_displacement_segments(
        np.empty((2, 4, 4, 2), dtype=np.float32),
        topology,
        retain_maps=True,
        working_memory_mb=1.0,
        array_allocator=allocate,
    )

    assert len(allocated) == 9
    assert result.maps is not None
    assert all(np.isnan(value).all() for value in result.__dict__.values())


def test_displacement_timing_prefers_moments_then_falls_back_to_bands(
    monkeypatch,
) -> None:
    moments = SimpleNamespace(
        cardiac_cycle_source="none",
        cycle_boundary_indexes=np.asarray([], dtype=np.int32),
    )
    bands = SimpleNamespace(
        cardiac_cycle_source="artery",
        cycle_boundary_indexes=np.asarray([0, 4], dtype=np.int32),
        cycle_durations_seconds=np.asarray([0.4], dtype=np.float32),
    )
    monkeypatch.setattr(
        displacement_timing,
        "velocity_workflows",
        lambda _ctx: SimpleNamespace(moments=moments, bandratio=bands),
    )

    selected = displacement_timing.select_displacement_timing(object())

    assert selected is not None
    assert selected.workflow == "frequency_bands"
    np.testing.assert_array_equal(selected.boundaries, [0, 4])
    assert selected.attrs["cycle_boundary_vessel_source"] == "artery"


def test_restored_payloads_write_through_the_current_h5_adapter() -> None:
    boundaries = np.asarray([0, 4])
    result = _result()
    outputs = pack_displacement_profile_outputs(
        {"artery": result},
        "level_set_motion",
        boundaries,
    )
    outputs.update(
        pack_displacement_segment_map_outputs(
            {"artery": result},
            "level_set_motion",
            boundaries,
        )
    )
    schema = EyeFlowOutputPaths.active()
    profile_path = (
        f"{schema.displacement_profiles_root}/LevelSetMotion/Artery/"
        "Transverse/Masked/tbkr/Profile"
    )
    map_path = f"{schema.displacement_map_root}/LevelSetMotion/Artery/PerSegment"

    with h5py.File("displacement-output.h5", "w", driver="core", backing_store=False) as h5:
        PipelineH5Output(h5).write_many(outputs)

        assert h5[profile_path].shape == (5, 4, 1, 1, 1)
        assert h5[profile_path].attrs["unit"] == "pixels"
        assert h5[map_path].shape == (5, 5, 4, 1, 1, 1, 2)
        assert h5[map_path].compression == "gzip"


def test_real_process_registration_product_feeds_the_production_profile_graph(tmp_path):
    """Exercise both actual pipelines, calculations, handoff and final schemas."""
    from test_topology_chunks import _prepared
    from test_velocity_pipeline import _retinal_velocity
    from input_output.output_manager import OutputManager
    from input_output.runs.layout import HoloRunLayout
    from input_output.schema.doppler_view import DV_ARTERY_MASK_PATH, DV_VEIN_MASK_PATH
    from pipeline_engine.base import PipelineDescriptor
    from pipeline_engine.resources import ResourceBudget
    from pipeline_engine.runtime import run_pipelines_to_output
    from pipelines.topology_core.runner import TOPOLOGY_CORE_STATE
    from pipelines.velocity.models import VelocityWorkflows
    from pipelines.velocity.runner import VELOCITY_WORKFLOWS_STATE

    catalog, _ = load_pipeline_catalog()
    registration = next(item for item in catalog if item.name == "displacement_map")
    profile = next(item for item in catalog if item.name == "displacement_profiles")
    if not registration.available:
        import pytest
        pytest.skip("Optional displacement registration dependencies unavailable")
    hd_path, dv_path = tmp_path / "HD.h5", tmp_path / "DV.h5"
    y, x = np.mgrid[:9, :9]
    image = np.exp(-((x - 4)**2 + (y - 4)**2) / 4).astype(np.float32)
    with h5py.File(hd_path, "w") as hd:
        hd["moment0"] = np.stack([np.roll(image, frame % 2, axis=1) for frame in range(6)])
    with h5py.File(dv_path, "w") as dv:
        mask = np.zeros((9, 9), np.uint8)
        mask[2:7, 2:7] = 1
        dv[DV_ARTERY_MASK_PATH] = mask
        dv[DV_VEIN_MASK_PATH] = mask

    def inputs(ctx):
        ctx.state.set(TOPOLOGY_CORE_STATE, {name: _prepared(angle=0)
                                          for name in ("artery", "vein")})
        ctx.state.set(VELOCITY_WORKFLOWS_STATE, VelocityWorkflows(_retinal_velocity(), None))

    source = PipelineDescriptor("shared_fixture", "shared_fixture", True,
        dag_produces=("prepared_topology", "velocity_workflows"),
        pipeline_factory=lambda: SimpleNamespace(name="shared_fixture", run=inputs))
    budget = ResourceBudget(cpu_threads=8, memory_bytes=2 * 1024**3)
    with budget.activate():
        final = run_pipelines_to_output(output_manager=OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo")),
            pipelines=[source, registration, profile], pipeline_options={"displacement_profiles": ("segment_maps",)},
            holodoppler_h5=hd_path, doppler_vision_h5=dv_path, copy_hd_passthrough=False)
    schema = EyeFlowOutputPaths.active()
    with h5py.File(final) as result:
        for vessel in ("Artery", "Vein"):
            group = result[f"{schema.displacement_profiles_root}/LevelSetMotion/{vessel}"]
            assert group["Transverse/Masked/tbkr/Profile"].attrs["unit"] == "pixels"
            assert result[f"{schema.displacement_map_root}/LevelSetMotion/{vessel}/PerSegment"].shape[-1] == 2
    assert budget.cpu_used == 0 and not any(budget.used.values())
