"""Tests for PDF report runner output paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.reports.pdf_data import extract_parameters_from_h5  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5_output import PipelineH5Output  # noqa: E402


class PdfReportRunnerPathTests(unittest.TestCase):
    def test_report_inherits_velocity_unit_with_dataset_override(self) -> None:
        from input_output.reports.pdf_data import _velocity_dataset_unit

        with h5py.File("report_units.h5", "w", driver="core", backing_store=False) as f:
            f.attrs["velocity_unit"] = "wrong root fallback"
            output = PipelineH5Output(f, processing_root="ProcessingBandRatio",
                                      provenance={"velocity_unit": "mm/s"})
            output.write("ProcessingMoments/signal", [1])
            dataset = output.get("ProcessingMoments/signal")
            self.assertEqual("mm/s", _velocity_dataset_unit(f, dataset))
            dataset.attrs["unit"] = "local override"
            self.assertEqual("local override", _velocity_dataset_unit(f, dataset))


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
                    schema.cardiac_cycle.systolic_cycle_duration_seconds,
                    data=np.asarray([0.5, 0.5]),
                )
                output_h5.create_dataset(
                    schema.cardiac_cycle.spectral_heart_rate_bpm,
                    data=np.asarray(96.795),
                )
                output_h5.create_dataset(
                    f"{schema.waveform_shape_metrics_root}/artery/global/raw/RI",
                    data=np.asarray([0.2, 0.4]),
                )

            parameters = extract_parameters_from_h5(path)

        self.assertEqual(4.0, parameters["Average_Arterial_Velocity"]["value"])
        self.assertAlmostEqual(96.795, parameters["heart_beat"]["value"])
        self.assertAlmostEqual(0.3, parameters["ARI"]["value"])


if __name__ == "__main__":
    unittest.main()
