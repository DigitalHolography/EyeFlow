"""DAG and metadata contracts for the selectable BVR pipeline."""

from __future__ import annotations

from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines import load_pipeline_catalog
from pipelines.velocity_analysis.builder import _segments_required


def _dag() -> PipelineDAG:
    load_pipeline_catalog()
    return PipelineDAG(PIPELINE_REGISTRY.values())


def test_bvr_option_dependencies_are_resolved_independently() -> None:
    dag = _dag()

    both = dag.resolve_targets(
        ["blood_volume_rate"],
        pipeline_options={
            "blood_volume_rate": ("gradient_edges", "masked_edges")
        },
    ).names
    gradient = dag.resolve_targets(
        ["blood_volume_rate"],
        pipeline_options={"blood_volume_rate": ("gradient_edges",)},
    ).names
    masked = dag.resolve_targets(
        ["blood_volume_rate"],
        pipeline_options={"blood_volume_rate": ("masked_edges",)},
    ).names
    neither = dag.resolve_targets(
        ["blood_volume_rate"],
        pipeline_options={"blood_volume_rate": ()},
    ).names

    assert "spatial_gradient_moment0" in both
    assert "velocity_analysis" in both
    assert "spatial_gradient_moment0" in gradient
    assert "velocity_analysis" in gradient
    assert "spatial_gradient_moment0" not in masked
    assert "velocity_analysis" in masked
    assert neither == (
        "topology_core",
        "velocity",
        "velocity_analysis",
        "blood_volume_rate",
    )


def test_spatial_gradient_alone_does_not_schedule_velocity_analysis() -> None:
    names = _dag().resolve_targets(["spatial_gradient_moment0"]).names

    assert "velocity" in names
    assert "topology_core" in names
    assert "velocity_analysis" not in names


def test_bvr_defaults_to_mask_derived_outputs_only() -> None:
    load_pipeline_catalog()
    descriptor = PIPELINE_REGISTRY["blood_volume_rate"]

    assert descriptor.visibility == "visible"
    assert descriptor.default_selected
    assert {option.name for option in descriptor.options if option.default_enabled} == {
        "masked_edges",
    }


def test_each_bvr_family_requests_canonical_segment_analysis() -> None:
    class Context:
        def __init__(self, bvr_options):
            self.bvr_options = frozenset(bvr_options)

        def pipeline_scheduled(self, name):
            return name in {"blood_volume_rate", "velocity_analysis"}

        def options_for(self, name):
            return self.bvr_options if name == "blood_volume_rate" else frozenset()

    assert _segments_required(Context(("gradient_edges",)))
    assert _segments_required(Context(("masked_edges",)))

