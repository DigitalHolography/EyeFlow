"""DAG and metadata contracts for the selectable BVR pipeline."""

from __future__ import annotations

from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines import load_pipeline_catalog
from pipelines.waveform_velocity_core.runner import _per_beat_required


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
    assert "waveform_velocity_core" in both
    assert "spatial_gradient_moment0" in gradient
    assert "waveform_velocity_core" in gradient
    assert "spatial_gradient_moment0" not in masked
    assert "waveform_velocity_core" in masked
    assert neither == ("blood_volume_rate",)


def test_spatial_gradient_alone_does_not_schedule_velocity_core() -> None:
    names = _dag().resolve_targets(["spatial_gradient_moment0"]).names

    assert "heartbeat_core" in names
    assert "topology_core" in names
    assert "waveform_velocity_core" not in names


def test_bvr_is_visible_default_selected_with_default_enabled_families() -> None:
    load_pipeline_catalog()
    descriptor = PIPELINE_REGISTRY["blood_volume_rate"]

    assert descriptor.visibility == "visible"
    assert descriptor.default_selected
    assert {option.name for option in descriptor.options if option.default_enabled} == {
        "gradient_edges",
        "masked_edges",
    }


def test_only_mask_family_requests_per_beat_segment_velocity() -> None:
    class Context:
        def __init__(self, bvr_options):
            self.bvr_options = frozenset(bvr_options)

        def pipeline_scheduled(self, name):
            return name in {"blood_volume_rate", "waveform_velocity_core"}

        def option_enabled(self, name, *, pipeline):
            return pipeline == "blood_volume_rate" and name in self.bvr_options

        def options_for(self, name):
            return self.bvr_options if name == "blood_volume_rate" else frozenset()

    assert not _per_beat_required(Context(("gradient_edges",)))
    assert _per_beat_required(Context(("masked_edges",)))

