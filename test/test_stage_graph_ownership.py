"""Graph validation and physical ownership contracts without scientific mocks."""

from dataclasses import dataclass
from threading import Thread

import numpy as np
import pytest

from data_manager import DataManager
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec
from pipeline_engine.storage import ProductStorage


def noop(ctx):
    pass


def test_plan_is_deterministic_and_keeps_dependency_kinds_separate():
    graph = StageGraph((
        StageSpec("a", noop, products=(ProductSpec("x"),)),
        StageSpec("b", noop, requires=("x",), products=(ProductSpec("y"),)),
        StageSpec("c", noop, optional_requires=("x",), after=("b",)),
        StageSpec("independent", noop),
    ))
    plan = graph.plan()
    assert plan.order == ("a", "b", "c", "independent")
    assert plan.required["b"] == ("x",)
    assert plan.optional["c"] == ("x",)
    assert plan.successful_stages["c"] == ("b",)
    assert plan.dependencies("c") == ("a", "b")
    assert plan.dependents("a") == ("b", "c")
    assert plan.ancestors("c") == ("a", "b")
    assert plan.producer_by_product["x"] == "a"
    with pytest.raises(TypeError):
        plan.stage_by_id["other"] = graph.stages[0]


def test_local_planning_does_not_construct_global_descriptors(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("A local stage graph used the global pipeline planner")
    monkeypatch.setattr("pipeline_engine.dag.PipelineDAG", forbidden)
    monkeypatch.setattr("pipeline_engine.base.PipelineDescriptor", forbidden)
    assert StageGraph((StageSpec("a", noop),)).plan().order == ("a",)


@pytest.mark.parametrize("stages", [
    (StageSpec("", noop),),
    (StageSpec("a", noop), StageSpec("a", noop)),
    (StageSpec("a", None),),
    (StageSpec("a", noop, products=(ProductSpec(""),)),),
    (StageSpec("a", noop, products=(ProductSpec("x"), ProductSpec("x"))),),
    (StageSpec("a", noop, requires=("missing",)),),
    (StageSpec("a", noop, after=("missing",)),),
    (StageSpec("a", noop, after=("a",)),),
    (StageSpec("a", noop, requires=("x",), products=(ProductSpec("x"),)),),
    (StageSpec("a", noop, after=("b",)), StageSpec("b", noop, after=("a",))),
    (StageSpec("a", noop, requires=("x",), optional_requires=("x",)),
     StageSpec("b", noop, products=(ProductSpec("x"),))),
    (StageSpec("a", noop, products=(ProductSpec("x", "Data/x", pack=noop),)),),
    (StageSpec("a", noop, products=(ProductSpec("x", "Data/x"), ProductSpec("y", "Data"))),),
])
def test_invalid_definitions_are_rejected(stages):
    with pytest.raises(ValueError):
        StageGraph(stages).plan()


def test_definition_metadata_is_recursively_immutable():
    metadata = {"nested": {"axes": ["t", "x"]}}
    spec = ProductSpec("x", attrs=metadata)
    metadata["nested"]["axes"].append("y")
    assert spec.attrs["nested"]["axes"] == ("t", "x")
    with pytest.raises(TypeError):
        spec.attrs["nested"]["axes"] = ()
    with pytest.raises((TypeError, AttributeError)):
        spec.attrs["nested"]["axes"].append("y")


@pytest.mark.parametrize("large", [False, True])
def test_snapshot_protects_producer_aliases_and_reuses_consumers(large):
    size = 1024 * 1024 if large else 8
    producer = np.arange(size, dtype=np.float32)
    alias = producer.view()
    alias.flags.writeable = False
    with ProductStorage() as storage:
        manager = DataManager(storage)
        handle = manager.accept(alias)
        producer[:] = -1
        first, second = manager.read(handle), manager.read(handle)
        np.testing.assert_array_equal(first, np.arange(size, dtype=np.float32))
        assert first is second
        assert not first.flags.writeable
        with pytest.raises(ValueError):
            first[:] = 1
        with pytest.raises(ValueError):
            first.flags.writeable = True
        writable = manager.writable(handle)
        writable[:] = 5
        assert manager.read(handle)[0] == 0
        if large:
            assert storage.metrics.snapshot()["npy_writes"] == 1
            adopted = manager.accept(first)
            assert adopted.storage == handle.storage
            assert storage.metrics.snapshot()["npy_writes"] == 1


def test_mutable_scientific_envelopes_are_private_without_array_copies():
    @dataclass
    class Result:
        array: np.ndarray
        label: str

    with ProductStorage() as storage:
        manager = DataManager(storage)
        original = Result(np.ones(1024 * 1024, dtype=np.float32), "original")
        handle = manager.accept(original)
        original.label = "producer mutation"
        first = manager.read(handle)
        first.label = "consumer mutation"
        second = manager.read(handle)
        assert second.label == "original"
        assert first.array is second.array
        assert storage.metrics.snapshot()["npy_writes"] == 1


def test_exposed_local_scientific_types_travel_by_small_stable_reference():
    import pickle
    from dataclasses import dataclass
    from pipeline_engine.storage import ValueReference
    @dataclass
    class Result:
        values: np.ndarray
    with ProductStorage() as storage:
        data = DataManager(storage)
        handle = data.accept(Result(np.ones(1024 * 1024, dtype=np.float32)), transportable=True)
        assert isinstance(handle.storage, ValueReference)
        assert len(pickle.dumps(handle, protocol=5)) < 1024
        received = pickle.loads(pickle.dumps(handle, protocol=5))
        consumer = DataManager(storage)
        consumer.acquire(received)
        np.testing.assert_array_equal(consumer.read(received).values, 1)
        assert storage.metrics.snapshot()["npy_writes"] == 1


def test_scientific_result_reads_do_not_rerun_constructors_or_copy_arrays():
    from dataclasses import dataclass, field
    constructed = []
    @dataclass
    class Result:
        array: np.ndarray
        computed: int = field(init=False)
        def __post_init__(self):
            constructed.append(True)
            self.array = self.array.copy()
            self.computed = 17
    with ProductStorage() as storage:
        data = DataManager(storage)
        handle = data.accept(Result(np.ones(32768, dtype=np.float32)))
        first, second = data.read(handle), data.read(handle)
        assert constructed == [True]
        assert first.array is second.array
        assert second.computed == 17
        assert data.writable(handle).computed == 17
        assert constructed == [True]


def test_storage_owner_and_shape_are_validated(tmp_path):
    from dataclasses import replace
    with ProductStorage() as storage:
        manager = DataManager(storage)
        handle = manager.accept(np.ones(1024 * 1024, dtype=np.float32))
        assert manager.acquire(handle) is handle
        with pytest.raises(ValueError, match="owner"):
            manager.acquire(replace(handle, owner_id="other-run"))
        with pytest.raises(ValueError, match="shape/dtype"):
            manager.acquire(replace(handle, storage=replace(handle.storage, shape=(1,))))


def test_repeated_aliases_in_one_product_use_one_snapshot_and_one_mapping():
    with ProductStorage() as storage:
        manager = DataManager(storage)
        value = np.ones(1024 * 1024, dtype=np.float32)
        handle = manager.accept({"a": value, "b": value})
        accepted = manager.read(handle)
        assert accepted["a"] is accepted["b"]
        assert storage.metrics.snapshot()["npy_writes"] == 1


def test_owned_resources_all_close_even_if_one_close_fails():
    from data_manager import InvocationResource
    closed = []
    class Resource:
        def __init__(self, name, fail=False):
            self.name, self.fail = name, fail
        def close(self):
            closed.append(self.name)
            if self.fail:
                raise OSError("close failed")
    with ProductStorage() as storage:
        manager = DataManager(storage)
        manager.accept(InvocationResource(Resource("first")))
        manager.accept(InvocationResource(Resource("second", True)))
        with pytest.raises(OSError, match="close failed"):
            manager.close_resources()
        assert closed == ["second", "first"]
        manager.close_resources()
        assert closed == ["second", "first"]
