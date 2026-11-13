"""Artifact naming, archive contents, and migration of existing results."""

import zipfile

import pytest

from input_output.archives.zip_archive import write_zip_artifact
from scripts.migrate_artifacts import migrate_artifacts
from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.artifact_names import prefixed_artifact_path


def test_manager_prefixes_png_leaf_filename_once(tmp_path):
    kind = OutputType.PNG
    manager = OutputManager(
        HoloRunLayout.from_holo(tmp_path / "scan.v2.holo")
    ).for_artifact_namespace("moments")
    expected = manager.dir_for(kind) / "nested" / f"scan.v2_plot.{kind.value}"
    assert manager.path_for(kind, f"nested/plot.{kind.value}") == expected
    assert manager.path_for(kind, f"nested/scan.v2_plot.{kind.value}") == expected
    assert manager.path_for(kind, f"nested\\plot.{kind.value}") == expected
    assert manager.path_for(kind).name.startswith("scan.v2_")
    assert manager.path_for(OutputType.H5).name == "scan.v2_EF.h5"


def test_zip_artifact_writer_names_published_archive_and_preserves_members(tmp_path):
    filename = "outputs.zip"
    tree = tmp_path / "tree"
    tree.mkdir()
    (tree / "scan_plot.png").write_bytes(b"plot")
    archive = write_zip_artifact(tree, tmp_path / filename, stem="scan")
    assert archive == prefixed_artifact_path(tmp_path / filename, "scan")
    with zipfile.ZipFile(archive) as contents:
        assert contents.read("scan_plot.png") == b"plot"
    assert not list(tmp_path.glob(".*eyeflow-staging-*"))


def test_migration_detects_collision_before_renaming_anything(tmp_path):
    root = tmp_path / "scan_EF"
    root.mkdir()
    (root / "first.png").write_bytes(b"first")
    (root / "plot.png").write_bytes(b"original")
    (root / "scan_plot.png").write_bytes(b"collision")
    with pytest.raises(FileExistsError, match="overwrite"):
        migrate_artifacts(root)
    assert (root / "first.png").read_bytes() == b"first"
    assert (root / "plot.png").read_bytes() == b"original"
    assert (root / "scan_plot.png").read_bytes() == b"collision"
