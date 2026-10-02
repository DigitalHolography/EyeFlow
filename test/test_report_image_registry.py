"""Check producer-side image registration for downstream PDF reports."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.writers.png import FigureArtifactWriter  # noqa: E402
from pipelines.waveform_velocity_core.figures.maps import _export_maps  # noqa: E402
from pipelines.waveform_velocity_core.figures.systole import _export_systole_plots  # noqa: E402
from pipelines.waveform_velocity_core.figures.waveforms import _export_ri_pi_plots  # noqa: E402


class ReportImageRegistryTests(unittest.TestCase):
    def test_vessel_maps_register_actual_paths(self) -> None:
        writer = FigureArtifactWriter(SimpleNamespace(), "scan")
        paths = [Path("custom-artery.png"), Path("custom-vein.png")]
        context = SimpleNamespace(
            velocity_analysis={"fRMS_avg": np.ones((2, 2), dtype=np.float32)},
            artery_mask=np.ones((2, 2), dtype=bool),
            vein_mask=np.ones((2, 2), dtype=bool),
        )
        with (
            patch(
                "pipelines.waveform_velocity_core.figures.maps._heatmap_with_colorbar",
                return_value=[],
            ),
            patch(
                "pipelines.waveform_velocity_core.figures.maps._mask_background_map",
                side_effect=paths,
            ),
        ):
            written = _export_maps(writer, context)

        self.assertEqual(paths, written)
        self.assertEqual(
            {("vessel_map", "artery"): paths[0], ("vessel_map", "vein"): paths[1]},
            writer.artifacts,
        )

    def test_systole_plots_register_actual_paths(self) -> None:
        writer = FigureArtifactWriter(SimpleNamespace(), "scan")
        paths = [Path("custom-artery.png"), Path("custom-vein.png")]
        context = SimpleNamespace(
            cycle_boundary_indexes=np.asarray([1, 3]),
            time=np.arange(5, dtype=np.float32),
            velocity_analysis={
                "retinal_artery_velocity_signal_filtered": np.arange(5),
                "retinal_artery_velocity_signal_derivative": np.arange(5),
                "retinal_vein_velocity_signal_filtered": np.arange(5),
                "retinal_vein_velocity_signal_derivative": np.arange(5),
            },
        )
        with (
            patch(
                "pipelines.waveform_velocity_core.figures.systole._cycle_extrema",
                return_value=(np.asarray([1]), np.asarray([3])),
            ),
            patch(
                "pipelines.waveform_velocity_core.figures.systole._systole_plot",
                side_effect=paths,
            ),
        ):
            written = _export_systole_plots(writer, context)

        self.assertEqual(paths, written)
        self.assertEqual(
            {("systole", "artery"): paths[0], ("systole", "vein"): paths[1]},
            writer.artifacts,
        )

    def test_ri_plots_register_actual_paths(self) -> None:
        writer = FigureArtifactWriter(SimpleNamespace(), "scan")
        context = SimpleNamespace(
            cycle_boundary_indexes=np.asarray([1, 3]),
            time=np.arange(5, dtype=np.float32),
            velocity_analysis={
                "retinal_artery_velocity_signal_filtered": np.arange(5),
                "retinal_vein_velocity_signal_filtered": np.arange(5),
            },
        )
        cycles = SimpleNamespace(artery=np.arange(3), vein=np.arange(3))

        def fake_plot(_writer, suffix, *_args):
            return Path(f"actual-{suffix}")

        with (
            patch(
                "pipelines.waveform_velocity_core.figures.waveforms.paired_vessel_cycles",
                return_value=cycles,
            ),
            patch(
                "pipelines.waveform_velocity_core.figures.waveforms.pulse_metric_from_signal",
                return_value=object(),
            ),
            patch(
                "pipelines.waveform_velocity_core.figures.waveforms._ri_pi_plot",
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


if __name__ == "__main__":
    unittest.main()
