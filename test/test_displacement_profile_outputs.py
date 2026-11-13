"""Contracts for the optional displacement profile and map families."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from calculations.math import interpft_real
from calculations.vessel_segments.profiles.per_beat import (
    interpolate_segment_maps_per_beat,
)
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.displacement_map.runner import DisplacementMapPipelineConfig
from pipelines.displacement_map.stages import (
    displacement_graph,
    displacement_profiles_graph,
)
from pipelines.displacement_map import timing as displacement_timing
from pipelines.displacement_map.profile_outputs import (
    pack_displacement_profile_outputs,
)
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
        "displacement_map",
        "prepared_topology",
        "velocity_workflows",
    )

    base_graph = displacement_graph(DisplacementMapPipelineConfig())
    assert [node.name for node in base_graph.nodes] == [
        "prepare",
        "artery",
        "vein",
        "finalize",
    ]

    profile_graph = displacement_profiles_graph()
    assert [node.name for node in profile_graph.nodes] == ["analysis"]
    assert [product.name for product in profile_graph.nodes[0].outputs] == [
        "profiles"
    ]

    map_graph = displacement_profiles_graph(retain_segment_maps=True)
    assert [product.name for product in map_graph.nodes[0].outputs] == [
        "profiles",
        "segment_maps",
    ]


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
