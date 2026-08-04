"""Tests for PDF report runner output paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.holo_run_layout import HoloRunLayout  # noqa: E402
from input_output.output_manager import OutputManager, OutputType  # noqa: E402
from input_output.reports.pdf_report import _extract_parameters_from_h5  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
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
            (hd_h5_path.parent.parent / "png").mkdir(parents=True)
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
            self.assertEqual(hd_h5_path.parent.parent / "png", kwargs["hd_png_dir"])
            self.assertNotIn("_HD_output_EF", str(kwargs["output_dir"]))
            self.assertIsNone(result)

    def test_extracts_active_waveform_output_schema(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.artery_per_beat.velocity_signal,
                    data=np.full((2, 3), 4.0),
                )
                output_h5.create_dataset(
                    schema.beat_period_seconds,
                    data=np.asarray([[0.5, 0.5]]),
                )
                output_h5.create_dataset(
                    schema.heartbeat.spectral_heart_rate_bpm,
                    data=np.asarray(96.795),
                )
                output_h5.create_dataset(
                    f"{schema.waveform_shape_metrics_root}/artery/global/raw/RI",
                    data=np.asarray([0.2, 0.4]),
                )

            parameters = _extract_parameters_from_h5([path])

        self.assertEqual(4.0, parameters["Average_Arterial_Velocity"]["value"])
        self.assertAlmostEqual(96.795, parameters["heart_beat"]["value"])
        self.assertAlmostEqual(0.3, parameters["ARI"]["value"])


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
