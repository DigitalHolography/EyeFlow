"""Artifact naming at writer boundaries and migration of existing results."""

import zipfile

import numpy as np
import pytest
from matplotlib.figure import Figure

from input_output.archives.zip_archive import create_zip_from_tree, write_zip_artifact
from input_output.artifact_migration import find_output_dirs, migrate_artifacts
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.avi import AviArtifactWriter, MjpegAviWriter
from input_output.writers.eps import EpsArtifactWriter, write_eps_file
from input_output.writers.png import PngArtifactWriter, write_png_file


@pytest.mark.parametrize("kind", list(OutputType))
@pytest.mark.parametrize("method", [None, "doppler_moments", "frequency_bands"])
def test_manager_prefixes_leaf_filename_once_for_every_format(tmp_path, kind, method):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.v2.holo"))
    if method is not None:
        manager = manager.for_workflow(method)
    expected = manager.dir_for(kind) / "nested" / f"scan.v2_plot.{kind.value}"
    assert manager.path_for(kind, f"nested/plot.{kind.value}") == expected
    assert manager.path_for(kind, f"nested/scan.v2_plot.{kind.value}") == expected
    assert manager.path_for(kind, f"nested\\plot.{kind.value}") == expected
    assert manager.path_for(kind).name.startswith("scan.v2_")
    assert manager.path_for(OutputType.H5).name == "scan.v2_EF.h5"


@pytest.mark.parametrize("writer_class", [PngArtifactWriter, EpsArtifactWriter, AviArtifactWriter])
def test_format_writer_keeps_subfolders_and_avoids_duplicate_prefix(tmp_path, writer_class):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "sample.holo"))
    writer = writer_class(manager)
    kind = {PngArtifactWriter: "png", EpsArtifactWriter: "eps", AviArtifactWriter: "avi"}[
        writer_class
    ]
    expected = manager.layout.ef_dir / kind / "nested" / f"sample_figure.{kind}"
    assert writer.path(f"nested/figure.{kind}") == expected
    assert writer.path(f"nested/sample_figure.{kind}") == expected
    labeled = writer_class(manager, "nested/artery")
    assert labeled.path(f"figure.{kind}") == expected.with_name(f"sample_artery_figure.{kind}")
    assert labeled.path(f"sample_figure.{kind}") == expected.with_name(
        f"sample_artery_figure.{kind}"
    )


def test_raw_png_eps_writers_accept_explicit_acquisition_stem(tmp_path):
    png = write_png_file(tmp_path / "image.png", np.zeros((2, 2), dtype=np.uint8), stem="scan")
    fig = Figure()
    fig.subplots().plot([0, 1], [1, 0])
    eps = write_eps_file(tmp_path / "figure.eps", fig, stem="scan")
    assert png.name == "scan_image.png" and png.is_file()
    assert eps.name == "scan_figure.eps" and eps.is_file()


def test_raw_avi_zip_writers_accept_explicit_acquisition_stem(tmp_path):
    tree = tmp_path / "tree"
    tree.mkdir()
    with MjpegAviWriter(tree / "video.avi", width=2, height=2, fps=10, stem="scan") as avi:
        avi.write_frame(np.zeros((2, 2, 3), dtype=np.uint8))
    assert avi.path.name == "scan_video.avi" and avi.path.is_file()
    archive = create_zip_from_tree(tree, tmp_path / "outputs.zip", stem="scan")
    assert archive.name == "scan_outputs.zip"
    with zipfile.ZipFile(archive) as contents:
        assert contents.namelist() == ["scan_video.avi"]


@pytest.mark.parametrize("filename", ["outputs.zip", "scan_outputs.zip", "custom.zip"])
def test_zip_artifact_writer_names_published_archive_and_preserves_members(tmp_path, filename):
    tree = tmp_path / "tree"
    tree.mkdir()
    (tree / "scan_plot.png").write_bytes(b"plot")
    archive = write_zip_artifact(tree, tmp_path / filename, stem="scan")
    assert archive == prefixed_artifact_path(tmp_path / filename, "scan")
    with zipfile.ZipFile(archive) as contents:
        assert contents.read("scan_plot.png") == b"plot"
    assert not list(tmp_path.glob(".*eyeflow-staging-*"))


def test_migration_renames_all_artifacts_preserving_content_and_directories(tmp_path):
    root = tmp_path / "scan_EF"
    originals = {
        "h5/scan_EF.h5": b"main HDF5",
        "png/moments/lumen/artery.png": b"PNG",
        "eps/bandratio/vein.eps": b"EPS",
        "avi/bandratio/movie.avi": b"AVI",
        "mp4/movie.mp4": b"MP4",
        "pdf/report.pdf": b"PDF",
        "h5/scratch.h5": b"auxiliary HDF5",
        "outputs.zip": b"ZIP",
        "diagnostic.json": b"JSON",
        "png/moments/scan_existing.png": b"already named",
    }
    for relative, contents in originals.items():
        file = root / relative
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes(contents)
    assert len(migrate_artifacts(root, dry_run=True)) == 8
    assert (root / "outputs.zip").exists()
    changes = migrate_artifacts(root)
    assert len(changes) == 8
    for relative, contents in originals.items():
        destination = prefixed_artifact_path(root / relative, "scan")
        assert destination.read_bytes() == contents
        assert destination.parent == (root / relative).parent
    assert migrate_artifacts(root) == ()
    assert find_output_dirs(tmp_path) == (root,)


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
