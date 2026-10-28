"""Artifact naming, archive contents, and migration of existing results."""

import zipfile

import numpy as np
import pytest
from matplotlib.figure import Figure

from input_output.archives.zip_archive import write_zip_artifact
from scripts.migrate_artifacts import migrate_artifacts
from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.eps import write_eps_file
from input_output.writers.png import write_png_file


@pytest.mark.parametrize("kind, namespace", [
    (OutputType.H5, None),
    (OutputType.PNG, "moments"),
])
def test_manager_prefixes_leaf_filename_once_for_every_format(tmp_path, kind, namespace):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.v2.holo"))
    if namespace is not None:
        manager = manager.for_artifact_namespace(namespace)
    expected = manager.dir_for(kind) / "nested" / f"scan.v2_plot.{kind.value}"
    assert manager.path_for(kind, f"nested/plot.{kind.value}") == expected
    assert manager.path_for(kind, f"nested/scan.v2_plot.{kind.value}") == expected
    assert manager.path_for(kind, f"nested\\plot.{kind.value}") == expected
    assert manager.path_for(kind).name.startswith("scan.v2_")
    assert manager.path_for(OutputType.H5).name == "scan.v2_EF.h5"


def test_raw_png_eps_writers_accept_explicit_acquisition_stem(tmp_path):
    png = write_png_file(tmp_path / "image.png", np.zeros((2, 2), dtype=np.uint8), stem="scan")
    fig = Figure()
    fig.subplots().plot([0, 1], [1, 0])
    eps = write_eps_file(tmp_path / "figure.eps", fig, stem="scan")
    assert png.name == "scan_image.png" and png.is_file()
    assert eps.name == "scan_figure.eps" and eps.is_file()


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
