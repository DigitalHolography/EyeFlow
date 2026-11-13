"""Local-space admission and the completed-file publication boundary."""

from hashlib import sha256
from dataclasses import replace
from types import SimpleNamespace
import pickle

import h5py
import numpy as np
import pytest

from input_output.writers.h5_publication import OutputPublisher
from pipeline_engine.data.metrics import ExecutionMetrics
from pipeline_engine.data.storage import ArrayReference, H5ArrayReference
from pipeline_engine.io.output.protocol import WriteCommand
from pipeline_engine.io.output.writer import _required_free_bytes, writer_main
from pipeline_engine.output.output_model import Record
from scratch_location import is_local_path


def test_publication_copies_once_to_destination_staging_and_reports_each_phase(tmp_path, monkeypatch):
    local, destination = tmp_path / "local.h5", tmp_path / "nas" / "final.h5"
    with h5py.File(local, "w") as file:
        file["value"] = np.arange(16)
    destination.parent.mkdir()
    destination.write_bytes(b"previous result")
    import shutil
    import os
    copyfile, replace = shutil.copyfile, os.replace
    events = []

    def copy(source, staging):
        assert source == local and staging.parent == destination.parent
        assert staging.name.startswith(".") and staging != destination
        assert destination.read_bytes() == b"previous result"
        events.append("copy")
        return copyfile(source, staging)

    def publish(staging, final):
        assert events == ["copy"] and final == destination
        with h5py.File(staging) as file:
            np.testing.assert_array_equal(file["value"][:], np.arange(16))
        events.append("replace")
        return replace(staging, final)

    monkeypatch.setattr("input_output.writers.h5_publication.shutil.copyfile", copy)
    monkeypatch.setattr("input_output.writers.h5_publication.os.replace", publish)
    metrics = ExecutionMetrics()
    assert OutputPublisher(destination, metrics).publish(local) == destination
    assert events == ["copy", "replace"] and local.exists()
    assert not list(destination.parent.glob(".ef-*"))
    values = metrics.snapshot()
    assert values["h5_publication_file_bytes"] == values["h5_publication_copy_bytes"] == local.stat().st_size
    assert values["h5_publication_copy_seconds"] >= 0
    assert values["h5_publication_replace_seconds"] >= 0
    assert values["h5_publication_copy_bytes_per_second"] > 0


def test_cancelled_publication_does_not_begin_network_copy(tmp_path, monkeypatch):
    destination = tmp_path / "final.h5"
    destination.write_bytes(b"previous result")
    monkeypatch.setattr("input_output.writers.h5_publication.shutil.copyfile",
                        lambda *_: pytest.fail("cancelled publication started copying"))
    metrics = ExecutionMetrics()
    with pytest.raises(OSError, match="cancelled"):
        OutputPublisher(destination, metrics, cancelled=lambda: True).publish(tmp_path / "absent.h5")
    assert destination.read_bytes() == b"previous result"
    assert metrics.snapshot()["h5_publication_failures"] == 1


def test_copy_failure_reports_elapsed_time_and_cleans_partial_staging(tmp_path, monkeypatch):
    local, destination = tmp_path / "local.h5", tmp_path / "final.h5"
    with h5py.File(local, "w") as file:
        file["value"] = 1
    destination.write_bytes(b"previous result")

    def failed_copy(source, staging):
        staging.write_bytes(b"incomplete copy")
        raise OSError("copy interrupted")

    monkeypatch.setattr("input_output.writers.h5_publication.shutil.copyfile", failed_copy)
    metrics = ExecutionMetrics()
    with pytest.raises(OSError, match="copy interrupted"):
        OutputPublisher(destination, metrics).publish(local)
    assert destination.read_bytes() == b"previous result" and local.exists()
    assert not list(tmp_path.glob(".ef-*"))
    values = metrics.snapshot()
    assert values["h5_publication_copy_seconds"] >= 0
    assert values["h5_publication_failures"] == values["h5_publication_staging_cleanups"] == 1


def test_dense_write_preflight_requires_full_payload_space_and_chunk_overhead(tmp_path):
    backing = tmp_path / "array.npy"
    backing.write_bytes(b"placeholder; do not read array contents")
    record = Record("dataset", ArrayReference(str(backing), (100_000_000,), np.dtype("f4")), committed=True)
    assert _required_free_bytes(record) >= 400_000_000 + 16 * 2**20
    record = replace(record, options={"dtype": "f8"})
    assert _required_free_bytes(record) >= 800_000_000 + 16 * 2**20
    record = replace(record, options={"dtype": "f8", "chunks": (1,)})
    assert _required_free_bytes(record) > 100_000_000 * 128


def test_native_sparse_hdf5_preflight_uses_allocated_storage(tmp_path):
    backing = tmp_path / "sparse.h5"
    shape = (1024, 1024, 256)
    with h5py.File(backing, "w") as file:
        data = file.create_dataset("data", shape=shape, dtype="f4", chunks=(2, 2, 2), fillvalue=np.nan)
        data[:2, :2, :2] = 1
    record = Record("dataset", H5ArrayReference(str(backing), shape, np.dtype("f4")), committed=True)
    assert backing.stat().st_size < _required_free_bytes(record) < 20 * 2**20


def test_writer_rejects_large_write_before_opening_backing_when_space_is_low(tmp_path, monkeypatch):
    backing = tmp_path / "array.npy"
    backing.write_bytes(b"not a numpy file")
    record = Record("dataset", ArrayReference(str(backing), (100_000_000,), np.dtype("f4")), committed=True)
    manifest = tmp_path / "manifest.pickle"
    manifest.write_bytes(pickle.dumps([("large", record)]))
    command = WriteCommand(1, str(manifest), sha256(manifest.read_bytes()).hexdigest())
    acknowledgements = []
    monkeypatch.setattr("pipeline_engine.io.output.writer.configure_native_threads", lambda *_: None)
    monkeypatch.setattr("pipeline_engine.io.output.writer.shutil.disk_usage",
                        lambda *_: SimpleNamespace(free=100 * 2**20))
    writer_main(SimpleNamespace(recv=lambda: command), SimpleNamespace(send=acknowledgements.append),
                str(tmp_path / "incomplete.h5"), 2**20)
    assert acknowledgements[0].sequence == 0
    assert "Insufficient local scratch" in acknowledgements[1].error
    with h5py.File(tmp_path / "incomplete.h5") as file:
        assert "large" not in file


def test_writer_rejects_network_scratch_before_creating_hdf5(tmp_path, monkeypatch):
    acknowledgements = []
    local = tmp_path / "remote" / "incomplete.h5"
    monkeypatch.setattr("pipeline_engine.io.output.writer.configure_native_threads", lambda *_: None)
    monkeypatch.setattr("pipeline_engine.io.output.writer.is_local_path", lambda *_: False)
    writer_main(None, SimpleNamespace(send=acknowledgements.append), str(local), 2**20)
    assert "scratch must be local" in acknowledgements[0].error
    assert not local.exists()


def test_writer_reports_separate_process_counters_for_each_operation(tmp_path, monkeypatch):
    contents = [
        [("value", Record("dataset", np.arange(16), committed=True))],
        ({"": Record(committed=True)}, ("value",)),
    ]
    commands = []
    for sequence, payload in enumerate(contents, start=1):
        manifest = tmp_path / f"{sequence}.pickle"
        manifest.write_bytes(pickle.dumps(payload))
        commands.append(WriteCommand(sequence, str(manifest), sha256(manifest.read_bytes()).hexdigest(), final=sequence == 2))
    iterator = iter(commands)
    acknowledgements = []
    monkeypatch.setattr("pipeline_engine.io.output.writer.configure_native_threads", lambda *_: None)
    writer_main(SimpleNamespace(recv=lambda: next(iterator)), SimpleNamespace(send=acknowledgements.append),
                str(tmp_path / "complete.h5"), 2**20)
    assert len(acknowledgements) == 3 and not any(ack.error for ack in acknowledgements)
    for ack in acknowledgements[1:]:
        assert ack.metrics["output_writer_cpu_seconds"] >= 0
        assert ack.metrics["output_writer_sampled_rss_bytes"] > 0
        assert all(value >= 0 for name, value in ack.metrics.items() if name.startswith("output_writer_os_"))


@pytest.mark.parametrize("filesystem", ["ext4", "xfs", "btrfs", "zfs", "tmpfs", "overlay", "apfs", "exfat"])
def test_posix_scratch_accepts_recognized_local_mounts(tmp_path, monkeypatch, filesystem):
    monkeypatch.setattr("scratch_location.os", SimpleNamespace(name="posix"))
    monkeypatch.setattr("psutil.disk_partitions", lambda **_: [SimpleNamespace(mountpoint=str(tmp_path), fstype=filesystem)])
    assert is_local_path(tmp_path / "scratch")


@pytest.mark.parametrize("filesystem", ["nfs", "nfs4", "cifs", "smb3", "ceph", "afs", "9p", "afpfs", "webdav", "fuse.rclone", "unrecognized"])
def test_posix_scratch_rejects_remote_and_unrecognized_mounts(tmp_path, monkeypatch, filesystem):
    monkeypatch.setattr("scratch_location.os", SimpleNamespace(name="posix"))
    monkeypatch.setattr("psutil.disk_partitions", lambda **_: [SimpleNamespace(mountpoint=str(tmp_path), fstype=filesystem)])
    assert not is_local_path(tmp_path / "scratch")


def test_posix_scratch_uses_longest_mount_and_fails_closed_without_mount_data(tmp_path, monkeypatch):
    monkeypatch.setattr("scratch_location.os", SimpleNamespace(name="posix"))
    mounts = [SimpleNamespace(mountpoint=str(tmp_path), fstype="ext4"),
              SimpleNamespace(mountpoint=str(tmp_path / "nas"), fstype="smb3")]
    monkeypatch.setattr("psutil.disk_partitions", lambda **_: mounts)
    assert not is_local_path(tmp_path / "nas" / "scratch")
    assert is_local_path(tmp_path / "local" / "scratch")
    mounts.clear()
    assert not is_local_path(tmp_path / "scratch")

    def unavailable(**_):
        raise OSError("mount inspection unavailable")

    monkeypatch.setattr("psutil.disk_partitions", unavailable)
    assert not is_local_path(tmp_path / "scratch")
