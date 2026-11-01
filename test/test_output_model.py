"""Output serialization must preserve metadata and publish atomically."""

from pathlib import Path

import h5py
import numpy as np
import pytest

from input_output.writers.h5_initialization import _copy_source_object
from pipeline_engine.output_model import OutputModel
from pipeline_engine.storage import ProductStorage


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
