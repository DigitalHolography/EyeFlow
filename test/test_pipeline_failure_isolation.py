"""Failure recovery follows active DAG dependencies and retains earlier products."""

import json
from types import SimpleNamespace

import h5py
import pytest
from test_run_service import _write_input

from input_output.output_manager import OutputType
from pipeline_engine.base import ExecutionVariant, PipelineDescriptor, PipelineOption
from pipeline_engine.run_service import execute_run, resolve_run_spec


def _descriptor(name, run, dependencies=(), *, variants=False, options=()):
    return PipelineDescriptor(
        name=name, description=name, available=True, dag_requires=dependencies,
        produces_execution_variants=variants, options=options,
        pipeline_factory=lambda: SimpleNamespace(name=name, run=run, source_path=__file__),
    )


@pytest.mark.parametrize("variants,fail_both", [(False, False), (True, False), (True, True)])
@pytest.mark.parametrize("failure_stage", ["calculate", "write"])
def test_failure_rolls_back_only_current_pipeline_and_continues_safe_nodes(
    tmp_path, variants, fail_both, failure_stage,
):
    holo = _write_input(tmp_path)
    calls = []

    def source(ctx):
        ctx.state.set("stable", 7)
        ctx.output.h5.set_attr("stable_attribute", "original")
        ctx.output.h5.write("ProcessingMoments/Upstream/value", 7)
        ctx.output.h5.set_attr("stable_attribute", "original")
        artifact = ctx.output.path_for(OutputType.PNG, "earlier.png")
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(b"earlier artifact")
        if variants:
            ctx.output.h5.write("ProcessingBandRatio/Upstream/value", 9)
            ctx.publish_execution_variants([
                ExecutionVariant(
                    name=name, state={"stable": value}, output_namespace=root,
                    artifact_namespace=name, provenance={},
                )
                for name, value, root in (
                    ("primary", 7, "ProcessingMoments"), ("alternate", 9, "ProcessingBandRatio"),
                )
            ])

    def failing(ctx):
        if variants and ctx.execution_variant.name == "alternate" and not fail_both:
            return {"ProcessingMoments/Product/value": 10}
        ctx.state.set("stable", -1)
        ctx.state.set("incomplete", True)
        ctx.output.h5.write("ProcessingMoments/Upstream/value", -1)
        ctx.output.h5.write("ProcessingMoments/Partial/value", -1)
        ctx.output.h5.set_attr("stable_attribute", "changed")
        artifact = ctx.output.path_for(OutputType.PNG, "partial.png")
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(b"incomplete artifact")
        if failure_stage == "write":
            return {"ProcessingMoments/Product/value": 1, "ProcessingMoments/Invalid/value": object()}
        raise ValueError("injected failure")

    def independent(ctx):
        branch = ctx.execution_variant.name if variants else "shared"
        calls.append(branch)
        assert "incomplete" not in ctx.state
        assert ctx.state["stable"] == (9 if branch == "alternate" else 7)
        assert ctx.output.h5.read("ProcessingMoments/Upstream/value") == ctx.state["stable"]
        return {"ProcessingMoments/Sibling/value": ctx.state["stable"]}

    def dependent(ctx):
        assert variants and not fail_both and ctx.execution_variant.name == "alternate"
        return {"ProcessingMoments/Dependent/value": 1}

    descriptors = [
        _descriptor("source", source, variants=variants),
        _descriptor("failing", failing, ("source",)),
        _descriptor("dependent", dependent, options=(
            PipelineOption("needs_product", "Needs product", dag_requires=("failing",)),
        )),
        _descriptor("sibling", independent, ("source",)),
    ]
    progress = []
    spec = resolve_run_spec(
        input_paths=[holo], target_names=["dependent", "sibling"], pipelines=descriptors,
        pipeline_options={"dependent": ("needs_product",)},
    )
    result = execute_run(spec, on_progress=lambda: progress.append(True))
    assert result.succeeded == (variants and not fail_both)
    assert len(progress) == len(spec.plan.names)
    assert calls == (["primary", "alternate"] if variants else ["shared"])
    assert len(result.outputs) == 1
    with h5py.File(result.outputs[0], "r") as output:
        assert output["ProcessingMoments/Upstream/value"][()] == 7
        assert output["ProcessingMoments/Sibling/value"][()] == 7
        assert "ProcessingMoments/Partial" not in output
        assert "ProcessingMoments/Product" not in output
        assert "ProcessingMoments/Dependent" not in output
        assert output.attrs["stable_attribute"] == "original"
        statuses = json.loads(output.attrs["pipeline_status"])
        status = statuses["variants"]["primary"] if variants else statuses["shared"]
        assert status["failing"] == "failed"
        assert status["dependent"] == "skipped"
        assert status["sibling"] == "succeeded"
        assert json.loads(output.attrs["pipeline_failures"])
        assert json.loads(output.attrs["pipeline_skips"])
        if variants and not fail_both:
            assert output["ProcessingBandRatio/Dependent/value"][()] == 1
            assert list(output.attrs["execution_variants_completed"]) == ["alternate"]
    root = spec.requests[0].output_manager.layout.ef_dir
    assert list(root.rglob("*earlier.png"))
    assert not list(root.rglob("*partial.png"))


def test_failed_variant_producer_does_not_prevent_unrelated_pipeline(tmp_path):
    holo = _write_input(tmp_path)

    def fail(ctx):
        raise ValueError("no variants available")

    def dependent(ctx):
        pytest.fail("A consumer must not run without its producer")

    spec = resolve_run_spec(
        input_paths=[holo], target_names=["dependent", "unrelated"], pipelines=[
            _descriptor("producer", fail, variants=True),
            _descriptor("dependent", dependent, ("producer",)),
            _descriptor("unrelated", lambda ctx: {"Unrelated/value": 42}),
        ],
    )
    result = execute_run(spec)
    assert not result.succeeded
    with h5py.File(result.outputs[0], "r") as output:
        assert output["Unrelated/value"][()] == 42
        assert json.loads(output.attrs["pipeline_status"])["shared"]["dependent"] == "skipped"
