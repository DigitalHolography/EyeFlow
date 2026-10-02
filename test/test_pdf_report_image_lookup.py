"""Tests for PDF report image lookup paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from PIL import Image

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.reports.pdf_report import (  # noqa: E402
    _try_load_ri_image,
    _try_load_systole_image,
    _try_load_vessel_image,
    generate_a4_report,
)


class PdfReportImageLookupTests(unittest.TestCase):
    def test_report_accepts_registered_images_with_arbitrary_names(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            output_h5 = root / "output.h5"
            with h5py.File(output_h5, "w"):
                pass
            vessel_map = root / "unrelated-actual-name.png"
            _write_rgb(vessel_map, (10, 20, 30))

            report = generate_a4_report(
                output_h5_path=output_h5,
                output_dir=root / "pdf",
                folder_name="scan",
                report_images={("vessel_map", "artery"): vessel_map},
            )

            self.assertTrue(report.is_file())
            self.assertGreater(report.stat().st_size, 1000)

    def test_loads_registered_png_paths_without_rebuilding_names(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            png_dir = Path(temp_dir) / "scan_EF" / "png"
            segmentation = png_dir / "actual-map.png"
            ri = png_dir / "actual-ri.png"
            systoles = png_dir / "actual-systoles.png"
            _write_rgb(segmentation, (10, 20, 30))
            _write_rgb(ri, (40, 50, 60))
            _write_rgb(systoles, (70, 80, 90))
            report_images = {
                ("vessel_map", "artery"): segmentation,
                ("ri", "artery"): ri,
                ("systole", "artery"): systoles,
            }

            np.testing.assert_array_equal(
                _try_load_vessel_image(report_images, None, None, "scan", "artery")[0, 0],
                [10, 20, 30],
            )
            np.testing.assert_array_equal(
                _try_load_ri_image(report_images, "artery")[0, 0],
                [40, 50, 60],
            )
            np.testing.assert_array_equal(
                _try_load_systole_image(report_images, "artery")[0, 0],
                [70, 80, 90],
            )


def _write_rgb(path: Path, color: tuple[int, int, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image = np.full((2, 2, 3), color, dtype=np.uint8)
    Image.fromarray(image, mode="RGB").save(path)


if __name__ == "__main__":
    unittest.main()
