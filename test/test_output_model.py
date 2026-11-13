"""Output serialization must preserve metadata and publish atomically."""

from pathlib import Path

import h5py
import numpy as np
import pytest

from input_output.writers.h5_initialization import _copy_source_object
from pipeline_engine.output.output_model import OutputModel, OutputImmutableError, OutputOwnershipError, OutputPathConflictError
from pipeline_engine.data.storage import ProductStorage


def test_incremental_dense_shape_is_sparse_atomic_immutable_and_reused(tmp_path):
    """A 1 TiB logical output needs only its two written chunks, never dense RAM."""
    import psutil
    from data_manager import DataManager
    from input_output.writers.h5 import write_value_dataset
    from input_output.schema.output_payload import DatasetValue

    process = psutil.Process()
    before = process.memory_info().rss
    shape = (1024, 1024, 1024, 256)
    options = {"chunks": (2, 2, 2, 1), "compression": "gzip", "compression_opts": 4, "shuffle": True}
    values = np.array([0., -0., np.nan, np.inf, -np.inf, 1., 2., 3.], np.float32).reshape(2, 2, 2)
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "sparse.h5")
        data = DataManager(storage)
        with model.checkpoint(owner="producer"):
            with model.build_arrays() as allocate:
                array = allocate(shape, h5_options=options)
                with pytest.raises(RuntimeError, match="Unfinished"):
                    data.accept(array)
                array[:2, :2, :2, 0] = values
                array[-2:, -2:, -2:, -1] = values
            accepted = data.read(data.accept(array))
            write_value_dataset(model, "Shared/Maps", DatasetValue(accepted, {"unit": "pixels"}, options))
        assert storage.store(accepted) == array.reference
        assert not storage.metrics.snapshot().get("npy_writes", 0)
        with pytest.raises(MemoryError, match="bounded slices"):
            np.asarray(accepted)
        with pytest.raises(OutputImmutableError), model.writing_as("producer"):
            model["Shared/Maps"][0, 0, 0, 0] = 42
        with pytest.raises(RuntimeError, match="construction failed"), model.checkpoint(owner="producer"):
            with model.build_arrays() as allocate:
                broken = allocate(shape, h5_options=options)
                broken[:2, :2, :2, 0] = values
                broken_path = broken._reference.path
                raise RuntimeError("construction failed")
        assert not Path(broken_path).exists()
        model.publish()
        with h5py.File(model.filename) as output:
            dataset = output["Shared/Maps"]
            assert dataset.shape == shape and dataset.dtype == np.dtype(np.float32)
            assert dataset.id.get_num_chunks() == 2
            assert dataset[:2, :2, :2, 0].tobytes() == values.tobytes()
            assert dataset[100:102, 100:102, 100:102, 99].tobytes() == np.full((2, 2, 2), np.nan, np.float32).tobytes()
            assert dataset.attrs["unit"] == "pixels"
            assert dataset.compression_opts == 4 and dataset.shuffle
        assert Path(model.filename).stat().st_size < 128 * 1024
        assert process.memory_info().rss - before < 64 * 2**20


def test_replacing_attribute_infers_dtype_again(tmp_path):
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "output.h5")
        model.attrs.create("value", 1, dtype=np.uint8)
        model.attrs["value"] = 300
        model.publish()
        with h5py.File(model.filename) as output:
            assert output.attrs["value"] == 300


def test_managed_source_copy_does_not_filter_string_metadata(tmp_path):
    with h5py.File(tmp_path / "source.h5", "w") as source, ProductStorage() as storage:
        dataset = source.create_dataset(
            "labels", data=["artery", "vein"], dtype=h5py.string_dtype("utf-8")
        )
        dataset.attrs["description"] = "source labels"
        model = OutputModel(storage, tmp_path / "output.h5")
        _copy_source_object(dataset, model, "copied")
        assert not model.records["copied"].options
        model.publish()
        with h5py.File(model.filename) as output:
            assert output["copied"].asstr()[:].tolist() == ["artery", "vein"]
            assert output["copied"].attrs["description"] == "source labels"


def test_failed_destination_copy_keeps_previous_output_and_cleans_staging(tmp_path, monkeypatch):
    final = tmp_path / "output.h5"
    final.write_bytes(b"previous output")

    def fail_copy(source, destination):
        Path(destination).write_bytes(b"partial")
        raise OSError("destination unavailable")

    monkeypatch.setattr("input_output.writers.h5_publication.shutil.copyfile", fail_copy)
    with ProductStorage() as storage:
        model = OutputModel(storage, final)
        model.create_dataset("value", data=1)
        with pytest.raises(OSError, match="destination unavailable"):
            model.publish()
        assert final.read_bytes() == b"previous output"
        assert not list(tmp_path.glob(".ef-*.h5"))
        assert not list(storage.directory.glob("*.h5"))


def test_dataset_commits_are_immutable_while_groups_remain_shared_and_authorized(tmp_path):
    """Two producers share a namespace; rejected writes roll back only new work."""
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "result.h5")
        with model.checkpoint(owner="A"):
            model.create_dataset("Temporary/value", data=1)
            del model["Temporary/value"]
            del model["Temporary"]
            data = model.create_dataset("Shared/Mean", data=[1, 2])
            data[0] = 3
            data.attrs["unit"] = "pixels"
            model.create_dataset("Shared/Mean", data=[3, 4])
            model["Shared/Mean"].attrs["unit"] = "pixels"
        records, cursor = model.committed_delta({})
        assert "Shared/Mean" in records
        assert not any(path.startswith("Temporary") for path in records)
        with model.checkpoint(owner="B"):
            model.create_dataset("Shared/Maximum", data=[9])
            model.require_group("Shared/Subgroup")
        records, cursor = model.committed_delta(cursor)
        assert "Shared/Maximum" in records and "Shared/Mean" not in records
        for owner, expected in (("A", OutputImmutableError), ("B", OutputOwnershipError)):
            for mutate in (
                lambda: model.create_dataset("Shared/Mean", data=[99]),
                lambda: model["Shared/Mean"].__setitem__(0, 99),
                lambda: model["Shared/Mean"].attrs.__setitem__("unit", "wrong"),
                lambda: model.__delitem__("Shared/Mean"),
            ):
                with pytest.raises(expected), model.checkpoint(owner=owner):
                    model.create_dataset("Rejected/value", data=99)
                    records, _ = model.committed_delta(cursor)
                    assert "Rejected/value" not in records
                    mutate()
                assert "Rejected" not in model
        with model.writing_as("B"):
            for mutate, expected in (
                (lambda: model.__delitem__("Shared"), OutputImmutableError),
                (lambda: model.create_dataset("Shared", data=1), OutputPathConflictError),
                (lambda: model.create_dataset("Shared/Mean/child", data=1), OutputPathConflictError),
                (lambda: model.require_group("Shared/Mean"), OutputPathConflictError),
                (lambda: model["Shared"].attrs.__setitem__("title", "unauthorized"), OutputOwnershipError),
                (lambda: model.attrs.__setitem__("pipeline_status", "unauthorized"), OutputOwnershipError),
            ):
                with pytest.raises(expected):
                    mutate()
        with model.checkpoint(owner="A"):
            model["Shared"].attrs["title"] = "authorized"
        np.testing.assert_array_equal(model["Shared/Mean"][:], [3, 4])
        np.testing.assert_array_equal(model["Shared/Maximum"][:], [9])
        assert model["Shared/Mean"].attrs["unit"] == "pixels"
        model.publish()
        with h5py.File(model.filename) as result:
            assert result["Shared"].attrs["title"] == "authorized"
            assert result["Shared/Mean"].attrs["unit"] == "pixels"


def _observe_committed_datasets(model):
    batches, cursor = [], {}

    def observed(current):
        records, next_cursor = current.committed_delta(cursor)
        datasets = {path: record for path, record in records.items() if record.kind == "dataset"}
        if datasets:
            batches.append(datasets)
        cursor.clear()
        cursor.update(next_cursor)

    model.bind_commit_listener(observed)
    return batches, observed


def test_commit_listener_hides_child_commits_until_outer_acceptance(tmp_path):
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "output.h5")
        batches, _ = _observe_committed_datasets(model)
        with pytest.raises(RuntimeError, match="outer rejected"), model.checkpoint(owner="producer"):
            with model.checkpoint() as child:
                model.create_dataset("Rejected/value", data=[1, 2])
                child.commit(("Rejected/value",))
                assert model.records["Rejected/value"].committed
                assert not batches
            assert not batches
            raise RuntimeError("outer rejected")
        assert "Rejected" not in model and not batches
        with model.checkpoint(owner="producer"):
            model.create_dataset("Accepted/value", data=[3, 4])
            assert not batches
        assert [set(batch) for batch in batches] == [{"Accepted/value"}]


def test_nonexposed_handoff_notifies_before_producer_returns_and_survives_later_failure(tmp_path):
    from input_output.output_manager import OutputManager
    from input_output.runs.layout import HoloRunLayout
    from pipeline_engine.context import PipelineContext, managed_invocation
    from pipeline_engine.execution.invocation import _pipeline_checkpoint
    from pipeline_engine.stage_graph import StageGraph, StageSpec, ProductSpec
    from pipeline_engine.stage_graph.executor import StageGraphExecutor

    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "output.h5")
        manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
        batches, _ = _observe_committed_datasets(model)
        exposed = []
        context = PipelineContext(work_h5=model, holodoppler_h5=None, doppler_vision_h5=None,
                                  output_manager=manager, pipeline_name="science", variables={},
                                  on_product_commit=lambda *args: exposed.append(args))

        def produce(ctx):
            ctx.dag.handoff("value", np.array([7., -0., np.nan]))
            assert [set(batch) for batch in batches] == [{"Scientific/value"}]
            assert not exposed  # A non-exposed handoff has no acceptance callback.
            raise ValueError("later calculation failed")

        graph = StageGraph((StageSpec("produce", produce, products=(
            ProductSpec("value", "Scientific/value", {"unit": "Hz"}),
        )),))
        with pytest.raises(RuntimeError, match="later calculation failed"), managed_invocation(context):
            with _pipeline_checkpoint(model, context.state.raw, manager, owner="science") as outer:
                StageGraphExecutor(graph, context, outer, _pipeline_checkpoint).run()
        assert [set(batch) for batch in batches] == [{"Scientific/value"}]
        assert model.records["Scientific/value"].committed
        assert batches[0]["Scientific/value"].attrs["unit"] == "Hz"
        np.testing.assert_array_equal(model["Scientific/value"][:], [7., -0., np.nan])


def test_commit_listener_late_binding_and_single_owner_do_not_repeat_dataset_deltas(tmp_path):
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "output.h5")
        with model.checkpoint(owner="producer"):
            model.create_dataset("value", data=3)
        batches, callback = _observe_committed_datasets(model)
        assert [set(batch) for batch in batches] == [{"value"}]
        model.bind_commit_listener(callback)
        with model.checkpoint(owner="producer"):
            model.require_group("Empty")
        assert len(batches) == 1
        with pytest.raises(RuntimeError, match="already bound"):
            model.bind_commit_listener(lambda current: None)
        model.bind_commit_listener(None)
        with model.checkpoint(owner="producer"):
            model.create_dataset("later", data=4)
        assert len(batches) == 1


def test_commit_listener_failure_preserves_accepted_records_and_state(tmp_path):
    from input_output.output_manager import OutputManager
    from input_output.runs.layout import HoloRunLayout
    from pipeline_engine.execution.invocation import _pipeline_checkpoint

    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "output.h5")
        failure, attempts = OSError("writer is unavailable"), []

        def failed(current):
            if "value" in current.committed_records():
                attempts.append(1)
                raise failure

        model.bind_commit_listener(failed)
        manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
        state = {"before": 1}
        with _pipeline_checkpoint(model, state, manager, owner="producer"):
            model.create_dataset("value", data=7)
            state["accepted"] = 2
        assert model.commit_listener_error is failure
        assert model.records["value"].committed and state == {"before": 1, "accepted": 2}
        with model.checkpoint(owner="producer"):
            model.create_dataset("later", data=8)
        assert attempts == [1] and model.records["later"].committed
