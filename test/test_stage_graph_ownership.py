"""Graph validation and physical ownership contracts without scientific mocks."""


import numpy as np
import pytest

from data_manager import DataManager
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec
from pipeline_engine.data.storage import ProductStorage


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
        first = manager.read(handle)
        counters = storage.metrics.snapshot()
        second = manager.read(handle)
        np.testing.assert_array_equal(first, np.arange(size, dtype=np.float32))
        np.testing.assert_array_equal(second, first)
        assert storage.metrics.snapshot()["array_copies"] == counters["array_copies"]
        assert storage.metrics.snapshot().get("npy_writes", 0) == counters.get("npy_writes", 0)
        assert not first.flags.writeable
        with pytest.raises(ValueError):
            first[:] = 1
        with pytest.raises(ValueError):
            first.flags.writeable = True
        writable = manager.writable(handle)
        writable[:] = 5
        assert manager.read(handle)[0] == 0
        copies = storage.metrics.snapshot()["array_copies"]
        manager.accept(first)
        assert storage.metrics.snapshot()["array_copies"] == copies
        aliases = manager.read(manager.accept({"first": first, "second": first}))
        assert np.shares_memory(aliases["first"], aliases["second"])
        assert storage.metrics.snapshot()["array_copies"] == copies
        if large:
            assert storage.metrics.snapshot()["npy_writes"] == 1


def test_exposed_local_scientific_types_travel_by_small_stable_reference():
    import pickle
    from dataclasses import dataclass
    from pipeline_engine.data.storage import ValueReference
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


@pytest.mark.parametrize("values", [np.array([3, 7], dtype=np.int32), np.array([-0., np.nan], dtype=np.float32)])
def test_small_exposed_arrays_keep_read_only_ownership_across_pickle_protocols(values):
    import pickle
    from pipeline_engine.data.storage import ArrayReference

    with ProductStorage() as storage:
        handle = DataManager(storage).accept(values, transportable=True)
        assert isinstance(handle.storage, ArrayReference)
        event_handle = pickle.loads(pickle.dumps(handle, protocol=4))
        manifest_handle = pickle.loads(pickle.dumps(handle, protocol=5))
        assert pickle.dumps(event_handle, protocol=5) == pickle.dumps(manifest_handle, protocol=5)
        consumer = DataManager(storage)
        consumer.acquire(event_handle)
        received = consumer.read(event_handle)
        assert received.dtype == values.dtype and received.tobytes() == values.tobytes()
        with pytest.raises(ValueError):
            received[0] = 1
        with pytest.raises(ValueError):
            received.flags.writeable = True


def test_scientific_result_reads_do_not_rerun_constructors_or_copy_arrays():
    from dataclasses import dataclass, field
    constructed = []
    @dataclass
    class Result:
        array: np.ndarray
        label: str = "original"
        computed: int = field(init=False)
        def __post_init__(self):
            constructed.append(True)
            self.array = self.array.copy()
            self.computed = 17
    with ProductStorage() as storage:
        data = DataManager(storage)
        original = Result(np.ones(32768, dtype=np.float32))
        handle = data.accept(original)
        original.label = "producer mutation"
        first = data.read(handle)
        first.label = "consumer mutation"
        counters = storage.metrics.snapshot()
        second = data.read(handle)
        assert constructed == [True]
        assert second.label == "original"
        np.testing.assert_array_equal(first.array, second.array)
        assert storage.metrics.snapshot()["array_copies"] == counters["array_copies"]
        assert storage.metrics.snapshot()["npy_writes"] == counters["npy_writes"]
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
        resource = Resource("first")
        manager.accept(InvocationResource(resource))
        handle = manager.accept((InvocationResource(resource), InvocationResource(resource)))
        with pytest.raises(TypeError, match="writable scientific copies"):
            manager.writable(handle)
        with pytest.raises(ValueError, match="cannot cross"):
            manager.accept({"nested": InvocationResource(resource)}, transportable=True)
        manager.accept(InvocationResource(Resource("second", True)))
        with pytest.raises(OSError, match="close failed"):
            manager.close_resources()
        assert closed == ["second", "first"]
        manager.close_resources()
        assert closed == ["second", "first"]
