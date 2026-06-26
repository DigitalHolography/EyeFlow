"""Tests for PDF report image lookup paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

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


def _write_rgb(path: Path, color: tuple[int, int, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image = np.full((2, 2, 3), color, dtype=np.uint8)
    Image.fromarray(image, mode="RGB").save(path)


if __name__ == "__main__":
    unittest.main()
