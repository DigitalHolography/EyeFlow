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
