"""Tests for waveform-shape pulse PNG exporters."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity.signal_analysis.waveform import (  # noqa: E402
    average_cycle,
    pulse_metric_from_cycles,
    pulse_metric_from_signal,
)
from pipelines.velocity_analysis.artifacts.figures.systole import (  # noqa: E402
    _systole_plot,
)
from pipelines.velocity_analysis.artifacts.figures.waveforms import (  # noqa: E402
    _export_ri_pi_plots,
)


def _has_matplotlib() -> bool:
    try:
        import matplotlib  # noqa: F401
    except ModuleNotFoundError:
        return False
    return True


class PulsePngExporterTests(unittest.TestCase):
    def test_pulse_metrics_use_extrema_of_interpolated_mean_cycle(self) -> None:
        values = np.asarray(
            [1.0, 4.0, 1.0, 2.0, 1.0, 2.0, 1.0, 4.0, 1.0],
            dtype=np.float32,
        )
        peaks = np.asarray([0, 4, 9], dtype=np.int32)
        cycle = average_cycle(values, peaks, samples=5)
        metric = pulse_metric_from_signal(values, peaks, "PI", samples=5)

        np.testing.assert_allclose(
            cycle,
            [1.0, 2.625, 1.75, 2.625, 1.5],
            atol=1e-6,
        )
        self.assertIsNotNone(metric)
        self.assertAlmostEqual(2.625, metric.maximum)
        self.assertAlmostEqual(1.0, metric.minimum)
        self.assertAlmostEqual(1.9, metric.mean)
        self.assertAlmostEqual((2.625 - 1.0) / 1.9, metric.value)
        self.assertNotEqual(4.0, metric.maximum)

    @unittest.skipUnless(_has_matplotlib(), "matplotlib is not installed")
    def test_ri_pi_pngs_use_pdf_filtered_per_beat_sampling(self) -> None:
        beats = np.asarray(
            [[-2.0, 4.0, 0.0, 2.0], [2.0, 2.0, 2.0, 2.0]],
            dtype=np.float32,
        )
        expected = pulse_metric_from_cycles(beats, "RI", beat_axis=0)
        vessel = SimpleNamespace(
            signal=SimpleNamespace(velocity_signal_per_beat=beats),
        )
        context = SimpleNamespace(
            per_beat_result=SimpleNamespace(
                artery_filtered=vessel,
                vein_filtered=vessel,
            ),
            velocity_semantics=SimpleNamespace(
                axis_label="Velocity (mm/s)",
                value_suffix=" mm/s",
            ),
        )
        figures = {}
        artifacts = {}

        def save_figure(fig, suffix):
            figures[suffix] = fig
            return Path(suffix)

        writer = SimpleNamespace(
            save_figure=save_figure,
            register_artifact=lambda kind, vessel_name, path: artifacts.update(
                {(kind, vessel_name): path}
            ),
        )

        paths = _export_ri_pi_plots(writer, context)
        fig = figures["RI_v_artery.png"]
        try:
            ax = fig.axes[0]
            np.testing.assert_allclose(ax.lines[0].get_xdata(), [0.0, 0.25, 0.5, 0.75])
            np.testing.assert_allclose(ax.lines[0].get_ydata(), [0.0, 3.0, 1.0, 2.0])
            self.assertEqual("Cardiac phase", ax.get_xlabel())
            self.assertIn(
                f"RI = {expected.value:.2f}",
                [text.get_text() for text in ax.texts],
            )
            self.assertEqual(4, len(paths))
            self.assertEqual(Path("RI_v_artery.png"), artifacts[("ri", "artery")])
        finally:
            import matplotlib.pyplot as plt

            for created in figures.values():
                plt.close(created)

    @unittest.skipUnless(_has_matplotlib(), "matplotlib is not installed")
    def test_systole_plot_has_title_and_no_acceleration_axis(self) -> None:
        captured = {}

        def save_figure(fig, suffix):
            captured["figure"] = fig
            return Path(suffix)

        writer = SimpleNamespace(save_figure=save_figure)
        path = _systole_plot(
            writer,
            "find_systoles_indices_artery.png",
            np.arange(4, dtype=float),
            np.array([1.0, 3.0, 2.0, 1.0]),
            np.array([1]),
            np.array([1]),
            np.array([3]),
            title="Artery systole indices",
        )
        fig = captured["figure"]
        try:
            self.assertEqual(Path("find_systoles_indices_artery.png"), path)
            self.assertEqual(1, len(fig.axes))
            self.assertEqual("Artery systole indices", fig.axes[0].get_title())
            self.assertEqual("black", fig.axes[0].lines[0].get_color())
            self.assertEqual(10, fig.axes[0].title.get_fontsize())
            self.assertEqual("normal", fig.axes[0].title.get_fontweight())
            self.assertEqual(["sans-serif"], fig.axes[0].title.get_fontfamily())
            self.assertNotIn("Acceleration", fig.axes[0].get_ylabel())
        finally:
            import matplotlib.pyplot as plt

            plt.close(fig)

if __name__ == "__main__":
    unittest.main()
