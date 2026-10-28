"""Rollback follows explicit ownership, never directory creation time."""

from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
import pytest

from input_output.artifact_ownership import ArtifactOwnership, artifact_transaction
from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.png import write_png_file
from pipeline_engine.runtime import _pipeline_checkpoint


def test_checkpoint_never_deletes_unregistered_files_created_during_a_failed_pipeline(tmp_path):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    manager.prepare()
    unrelated = manager.layout.ef_dir / "unrelated.txt"
    with h5py.File(tmp_path / "result.h5", "w") as h5:
        with pytest.raises(ValueError, match="fail"):
            with _pipeline_checkpoint(h5, {}, manager):
                owned = manager.path_for(OutputType.PNG, "owned.png")
                owned.parent.mkdir(parents=True)
                owned.write_bytes(b"owned by this invocation")
                # Another producer is not inside this invocation's scope.
                with ThreadPoolExecutor(max_workers=1) as other:
                    other.submit(unrelated.write_bytes, b"unrelated producer").result()
                raise ValueError("fail")
    assert not owned.exists()
    assert unrelated.read_bytes() == b"unrelated producer"


def test_reserving_another_stages_pending_path_does_not_transfer_ownership(tmp_path):
    registry = ArtifactOwnership()
    producer = registry.transaction(tmp_path)
    failing = registry.transaction(tmp_path)
    path = tmp_path / "committed_stage.png"
    producer.reserve(path, writing=True)
    failing.reserve(path)  # Still absent: the earlier export has not finished.
    path.write_bytes(b"export published during later stage")
    failing.rollback()
    assert path.read_bytes() == b"export published during later stage"


def test_shared_writer_cannot_overwrite_another_invocations_artifact(tmp_path):
    registry = ArtifactOwnership()
    path = tmp_path / "shared.png"
    with artifact_transaction(registry.transaction(tmp_path)):
        write_png_file(path, np.ones((2, 2), dtype=np.uint8))
    original = path.read_bytes()
    with artifact_transaction(registry.transaction(tmp_path)):
        with pytest.raises(FileExistsError, match="another pipeline"):
            write_png_file(path, np.zeros((2, 2), dtype=np.uint8))
    assert path.read_bytes() == original


def test_preexisting_and_outside_files_are_never_claimed_for_rollback(tmp_path):
    root = tmp_path / "result"
    root.mkdir()
    earlier = root / "earlier.png"
    outside = tmp_path / "external.png"
    earlier.write_bytes(b"earlier")
    outside.write_bytes(b"outside")
    transaction = ArtifactOwnership().transaction(root)
    transaction.reserve(earlier)
    transaction.reserve(outside)
    transaction.rollback()
    assert earlier.read_bytes() == b"earlier"
    assert outside.read_bytes() == b"outside"
