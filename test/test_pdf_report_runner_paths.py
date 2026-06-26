"""Tests for PDF report runner output paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.holo_run_layout import HoloRunLayout  # noqa: E402
from input_output.output_manager import OutputManager, OutputType  # noqa: E402
from pipelines.pdf_report.runner import run_pdf_report  # noqa: E402


class PdfReportRunnerPathTests(unittest.TestCase):
    def test_uses_layout_stem_and_pdf_report_dir(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_root = root / "input"
            output_root = root / "output"
            holo_path = input_root / "scan.holo"
            hd_h5_path = input_root / "scan" / "scan_HD" / "h5" / "scan_HD_output.h5"
            manager = OutputManager(
                HoloRunLayout.from_holo(holo_path, output_root=output_root)
            )
            output_h5_path = manager.path_for(OutputType.H5)
            (manager.layout.ef_dir / "png").mkdir(parents=True)
            ctx = _fake_context(manager, output_h5_path, hd_h5_path)

            def fake_generate_a4_report(**kwargs):
                return kwargs["output_dir"] / f"{kwargs['folder_name']}_report.pdf"

            with patch(
                "pipelines.pdf_report.runner.generate_a4_report",
                side_effect=fake_generate_a4_report,
            ) as generate:
                result = run_pdf_report(ctx)

            kwargs = generate.call_args.kwargs
            self.assertEqual("scan", kwargs["folder_name"])
            self.assertEqual(manager.layout.ef_dir / "pdf", kwargs["output_dir"])
            self.assertEqual(manager.layout.ef_dir / "png", kwargs["png_dir"])
            self.assertNotIn("_HD_output_EF", str(kwargs["output_dir"]))
            self.assertEqual(1, int(result.metrics["pdf_report_generated"]))


def _fake_context(manager: OutputManager, output_h5_path: Path, hd_h5_path: Path):
    logs: list[str] = []
    return SimpleNamespace(
        require_inputs=lambda *names: None,
        log=logs.append,
        output=SimpleNamespace(
            manager=manager,
            h5=SimpleNamespace(filename=str(output_h5_path)),
        ),
        inputs=SimpleNamespace(
            hd=SimpleNamespace(filename=str(hd_h5_path)),
            dv=SimpleNamespace(filename="dv.h5"),
        ),
    )


if __name__ == "__main__":
    unittest.main()
