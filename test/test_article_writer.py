"""Article pairs share dimensions/layout without changing generic PNG policy."""

import re

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from PIL import Image
import pytest

from input_output.article_formatting import create_article_figure, prepare_article_fill_layer
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager
from input_output.writers.article import save_article_figure_pair
from input_output.writers.png import PngArtifactWriter


@pytest.fixture(autouse=True)
def headless_figures():
    with mpl.rc_context({"backend": "Agg"}):
        yield


@pytest.mark.parametrize("name, namespace, stem, expected", [
    ("nested/plot", None, None, "nested/scan_plot"),
    ("nested/scan_plot.png", "moments", None, "nested/scan_plot"),
    ("nested/scan_plot.eps", "spectrograms", None, "nested/scan_plot"),
    ("plot", "bandratio", "velocity/artery", "velocity/scan_artery_plot"),
])
def test_article_pair_preserves_names_namespaces_and_prefixes(tmp_path, name, namespace, stem, expected):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    if namespace is not None:
        manager = manager.for_artifact_namespace(namespace)
    fig, ax = create_article_figure(width=1, square=True)
    ax.set_axis_off()
    paths = save_article_figure_pair(manager, fig, name, stem=stem, dpi=100)
    for path, kind in zip(paths, ("png", "eps")):
        base = manager.layout.ef_dir / kind
        if namespace:
            base /= namespace
        assert path == base / f"{expected}.{kind}"
        assert path.is_file()


def test_article_pair_freezes_layout_preserves_full_canvas_and_closes_once(tmp_path, monkeypatch, caplog):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    fig, ax = plt.subplots(figsize=(4, 2), constrained_layout=True)
    ax.plot([0, 1], [1, 2])
    prepare_article_fill_layer(ax)
    ax.fill_between([0, 1], [0.5, 1.5], [1.5, 2.5], alpha=0.3, zorder=0)
    ax.set(xlabel="Time (s)", ylabel="Power (a.u.)")
    original_save = fig.savefig
    saves = []

    def record_save(path, **options):
        assert plt.fignum_exists(fig.number)  # Neither save closes the shared figure early.
        saves.append((options, ax.get_position().bounds))
        return original_save(path, **options)

    monkeypatch.setattr(fig, "savefig", record_save)
    with mpl.rc_context({"savefig.bbox": "tight", "savefig.pad_inches": 0.6,
                         "savefig.transparent": True}):
        with caplog.at_level("WARNING", logger="matplotlib.backends.backend_ps"):
            png, eps = save_article_figure_pair(manager, fig, "article", dpi=100)
    assert len(saves) == 2 and not plt.fignum_exists(fig.number)
    assert [options["format"] for options, _ in saves] == ["png", "eps"]
    assert all(options["dpi"] == 100 and options["pad_inches"] == 0 for options, _ in saves)
    np.testing.assert_allclose(saves[0][1], saves[1][1])
    with Image.open(png) as image:
        assert image.size == (400, 200)
        assert image.convert("RGBA").getpixel((0, 0))[3] == 255
    contents = eps.read_text(encoding="latin-1")
    match = re.search(r"^%%HiResBoundingBox: ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)", contents, re.M)
    left, bottom, right, top = map(float, match.groups())
    np.testing.assert_allclose([right - left, top - bottom], [4 * 72, 2 * 72])
    assert "colorimage" in contents
    assert "does not support transparency" not in caplog.text


@pytest.mark.parametrize("failed_format", ["png", "eps"])
def test_article_pair_closes_figure_on_export_failure(tmp_path, monkeypatch, failed_format):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    def failing_save(path, **options):
        if options["format"] == failed_format:
            raise OSError("export failed")

    monkeypatch.setattr(fig, "savefig", failing_save)
    with pytest.raises(OSError, match="export failed"):
        save_article_figure_pair(manager, fig, "failure")
    assert not plt.fignum_exists(fig.number)


def test_generic_png_export_still_supports_tight_cropping_without_eps(tmp_path):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    fig = plt.figure(figsize=(4, 2))
    ax = fig.add_axes((0.4, 0.4, 0.2, 0.2))
    ax.plot([0, 1], [1, 2])
    ax.set_axis_off()
    png = PngArtifactWriter(manager).save_figure(fig, "diagnostic.png", dpi=100)
    with Image.open(png) as image:
        assert image.width < 400 and image.height < 200
    assert not (manager.layout.ef_dir / "eps").exists()
