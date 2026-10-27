"""Tests for PDF report image lookup paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from PIL import Image

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.reports.pdf_report import (  # noqa: E402
    _try_load_ri_image,
    _try_load_systole_image,
    _try_load_vessel_image,
)
from input_output.writers.png import FigureArtifactWriter  # noqa: E402
from pipelines.velocity_analysis.artifacts.figures.waveforms import _export_ri_pi_plots  # noqa: E402


class PdfReportImageLookupTests(unittest.TestCase):
    def test_loads_stem_prefixed_waveform_png_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            png_dir = Path(temp_dir) / "scan_EF" / "png"
            stem = "scan"

            segmentation = png_dir / f"{stem}_artery_seg_map_bkg.png"
            ri = png_dir / f"{stem}_RI_v_artery.png"
            systoles = png_dir / f"{stem}_find_systoles_indices_artery.png"
            _write_rgb(segmentation, (10, 20, 30))
            _write_rgb(ri, (40, 50, 60))
            _write_rgb(systoles, (70, 80, 90))

            np.testing.assert_array_equal(
                _try_load_vessel_image(png_dir, None, None, stem, "artery")[0, 0],
                [10, 20, 30],
            )
            np.testing.assert_array_equal(
                _try_load_ri_image(png_dir, stem, "artery")[0, 0],
                [40, 50, 60],
            )
            np.testing.assert_array_equal(
                _try_load_systole_image(png_dir, stem, "artery")[0, 0],
                [70, 80, 90],
            )

    def test_registered_paths_follow_actual_variant_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir) / "png" / "moments"
            images = {
                ("vessel_map", "artery"): folder / "custom-map.png",
                ("ri", "artery"): folder / "custom-ri.png",
                ("systole", "artery"): folder / "custom-systole.png",
            }
            for key, color in zip(images, ((10, 20, 30), (40, 50, 60), (70, 80, 90))):
                _write_rgb(images[key], color)

            np.testing.assert_array_equal(
                _try_load_vessel_image(None, None, None, "scan", "artery", report_images=images)[0, 0],
                [10, 20, 30],
            )
            np.testing.assert_array_equal(
                _try_load_ri_image(None, "scan", "artery", report_images=images)[0, 0],
                [40, 50, 60],
            )
            np.testing.assert_array_equal(
                _try_load_systole_image(None, "scan", "artery", report_images=images)[0, 0],
                [70, 80, 90],
            )

    def test_ri_export_records_actual_writer_paths(self) -> None:
        writer = FigureArtifactWriter(SimpleNamespace(), "scan")
        context = SimpleNamespace(
            cycle_boundary_indexes=np.asarray([1, 3]),
            time=np.arange(5, dtype=np.float32),
            velocity=SimpleNamespace(continuous=lambda _vessel: np.arange(5)),
            velocity_semantics=None,
        )
        cycles = SimpleNamespace(artery=np.arange(3), vein=np.arange(3))

        def fake_plot(_writer, suffix, *_args):
            return Path(f"actual-{suffix}")

        with (
            patch(
                "pipelines.velocity_analysis.artifacts.figures.waveforms.paired_vessel_cycles",
                return_value=cycles,
            ),
            patch(
                "pipelines.velocity_analysis.artifacts.figures.waveforms.pulse_metric_from_signal",
                return_value=object(),
            ),
            patch(
                "pipelines.velocity_analysis.artifacts.figures.waveforms._ri_pi_plot",
                side_effect=fake_plot,
            ),
        ):
            _export_ri_pi_plots(writer, context)

        self.assertEqual(
            {
                ("ri", "artery"): Path("actual-RI_v_artery.png"),
                ("ri", "vein"): Path("actual-RI_v_vein.png"),
            },
            writer.artifacts,
        )


def _write_rgb(path: Path, color: tuple[int, int, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image = np.full((2, 2, 3), color, dtype=np.uint8)
    Image.fromarray(image, mode="RGB").save(path)


if __name__ == "__main__":
    unittest.main()
