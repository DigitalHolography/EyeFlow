"""Tests for waveform-shape pulse PNG exporters."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pipelines.waveform_shape_metrics.velocity.figures.signal_inputs import (  # noqa: E402
    display_frequency,
    display_velocity,
    histogram_matrix,
    masked_video_signal,
)
from pipelines.waveform_shape_metrics.velocity.figures.spectrum import (  # noqa: E402
    correlation_data,
    paired_spectrum_analysis,
    spectrum_signal_analysis,
    synthetic_spectrum_analysis,
)
from calculations.blood_flow_velocity.signal_analysis.waveform import (  # noqa: E402
    arterial_waveform_analysis,
    average_cycle,
    cycle_extrema,
    pulse_metric,
)
from input_output.output_manager import OutputType  # noqa: E402
from input_output.writers.png import write_png_file  # noqa: E402
from calculations.dopplerview_analysis.vessel_velocity_estimator import (  # noqa: E402
    _masked_signal as _velocity_masked_signal,
)
from pipelines.waveform_shape_metrics.velocity.figures import (  # noqa: E402
    PULSE_PNG_SUFFIXES,
    export_pulse_pngs,
)
from pipelines.waveform_shape_metrics.velocity.figures.plotting import (  # noqa: E402
    _velocity_gradient_values,
)
from pipelines.waveform_shape_metrics.velocity.figures.velocity_maps import (  # noqa: E402
    _velocity_colorbar_vmax,
    _vessel_histogram_colormap,
)


class FakeOutput:
    available = True

    def __init__(self, root: Path, stem: str = "sample") -> None:
        self.root = root
        self.manager = SimpleNamespace(layout=SimpleNamespace(stem=stem))

    def path_for(self, output_type: OutputType, filename: str | None = None) -> Path:
        assert output_type is OutputType.PNG
        return self.root / "png" / (filename or "sample")

    def write_png(self, output, filename: str | None = None) -> Path:
        return write_png_file(self.path_for(OutputType.PNG, filename), output)


def _has_matplotlib() -> bool:
    try:
        import matplotlib  # noqa: F401
    except ModuleNotFoundError:
        return False
    return True


class PulsePngExporterTests(unittest.TestCase):
    def test_expected_suffixes_include_core_pulse_outputs(self) -> None:
        self.assertIn("f_artery_graph.png", PULSE_PNG_SUFFIXES)
        self.assertIn("find_systoles_indices_artery.png", PULSE_PNG_SUFFIXES)
        self.assertIn("ArterialSpectralAnalysis_v_artery.png", PULSE_PNG_SUFFIXES)
        self.assertIn("artery_seg_map_bkg.png", PULSE_PNG_SUFFIXES)
        self.assertIn("vein_seg_map_bkg.png", PULSE_PNG_SUFFIXES)
        self.assertIn("AVGflowVideoCombined.png", PULSE_PNG_SUFFIXES)

    @unittest.skipUnless(_has_matplotlib(), "matplotlib is not installed")
    def test_histogram_colormap_matches_matlab_vessel_ramp(self) -> None:
        artery_cmap = _vessel_histogram_colormap("artery")(np.linspace(0.0, 1.0, 4))
        vein_cmap = _vessel_histogram_colormap("vein")(np.linspace(0.0, 1.0, 4))

        np.testing.assert_allclose(artery_cmap[0, :3], [0.0, 0.0, 0.0], atol=1e-6)
        self.assertGreater(artery_cmap[1, 0], artery_cmap[1, 1])
        self.assertGreater(artery_cmap[2, 0], 0.8)
        self.assertGreater(artery_cmap[2, 1], 0.8)
        np.testing.assert_allclose(vein_cmap[0, :3], [0.0, 0.0, 0.0], atol=1e-6)
        self.assertGreater(vein_cmap[1, 2], vein_cmap[1, 0])

    def test_histogram_matrix_counts_masked_pixels_per_frame(self) -> None:
        velocity = np.asarray(
            [
                [[0.0, 10.0], [20.0, 30.0]],
                [[0.0, 10.0], [10.0, 30.0]],
            ],
            dtype=np.float32,
        )
        mask = np.asarray([[False, True], [True, False]])

        histo = histogram_matrix(velocity, mask, bins=2)

        self.assertEqual((2, 2), histo.counts.shape)
        np.testing.assert_array_equal(histo.counts[:, 0], [1, 1])
        np.testing.assert_array_equal(histo.counts[:, 1], [2, 0])
        self.assertEqual(2.0, histo.count_max)

    def test_correlation_data_uses_matlab_style_msc_settings(self) -> None:
        time = np.arange(96, dtype=np.float32) * 0.05
        artery = np.sin(2 * np.pi * 1.0 * time).astype(np.float32)
        vein = np.sin(2 * np.pi * 1.0 * time + 0.35).astype(np.float32)

        corr = correlation_data(artery, vein, 0.05)

        self.assertEqual(129, corr.coherence_freq.size)
        self.assertTrue(np.isfinite(corr.gamma_0))
        self.assertGreaterEqual(corr.gamma_0, 0.0)
        self.assertLessEqual(corr.gamma_0, 1.0)

    def test_display_unit_scaling_matches_png_labels(self) -> None:
        raw = np.asarray([0.0, 1000.0, 6000.0], dtype=np.float32)

        np.testing.assert_allclose(display_frequency(raw), [0.0, 1.0, 6.0])
        np.testing.assert_allclose(display_velocity(raw), raw)

    def test_average_cycle_and_extrema_are_generic_signal_helpers(self) -> None:
        values = np.asarray([0.0, 2.0, 1.0, -1.0, 0.0, 3.0, 1.0, -2.0], dtype=np.float32)
        peaks = np.asarray([1, 5], dtype=np.int32)

        cycle = average_cycle(values, peaks, samples=4)
        maxima, minima = cycle_extrema(values, peaks)

        np.testing.assert_allclose(cycle, values[1:5], atol=1e-6)
        np.testing.assert_array_equal(maxima, [1, 5])
        np.testing.assert_array_equal(minima, [0, 3])

    def test_spectrum_reports_frequency_bins_and_ordered_peaks(self) -> None:
        time = np.arange(128, dtype=np.float32) * 0.05
        values = np.sin(2 * np.pi * 1.25 * time).astype(np.float32)

        data = spectrum_signal_analysis(values, 0.05)

        self.assertEqual(data.frequencies.shape, data.magnitude.shape)
        self.assertEqual(data.frequencies.shape, data.phase.shape)
        self.assertTrue(np.all(np.diff(data.frequencies[data.peak_indexes]) >= 0))
        self.assertTrue(np.isfinite(data.heart_rate_bpm))

    def test_synthetic_and_paired_spectrum_analysis_outputs(self) -> None:
        time = np.arange(96, dtype=np.float32) * 0.05
        first = np.sin(2 * np.pi * 1.0 * time).astype(np.float32)
        second = np.sin(2 * np.pi * 1.0 * time + 0.35).astype(np.float32)
        beat_indexes = np.asarray([0, 24, 48, 72], dtype=np.int32)
        cycle = average_cycle(first, beat_indexes, 32)

        synthetic = synthetic_spectrum_analysis(cycle, 0.05, beat_indexes)
        paired = paired_spectrum_analysis(first, second, 0.05, beat_indexes)

        self.assertEqual(synthetic.frequencies.shape, synthetic.magnitude.shape)
        self.assertEqual(synthetic.frequencies.shape, synthetic.phase.shape)
        self.assertTrue(synthetic.peak_indexes.size > 0)
        self.assertEqual(paired.transfer.frequencies.shape, paired.transfer.transfer.shape)
        self.assertIsNotNone(paired.delay)
        self.assertGreaterEqual(paired.correlation.gamma_0, 0.0)
        self.assertLessEqual(paired.correlation.gamma_0, 1.0)

    def test_waveform_metric_and_marker_analysis_outputs(self) -> None:
        cycle = np.asarray([0.0, 3.0, 1.0, 2.0, 0.5, -1.0, 0.0, 1.0], dtype=np.float32)
        beat_indexes = np.asarray([0, 8, 16], dtype=np.int32)

        metric = pulse_metric(cycle, "RI")
        markers = arterial_waveform_analysis(cycle, beat_indexes, 0.1)

        self.assertAlmostEqual((metric.maximum - metric.minimum) / metric.maximum, metric.value)
        self.assertEqual(cycle.size, markers.gradient.size)
        self.assertGreaterEqual(markers.peak_indexes.size, 1)

    def test_masked_video_signal_averages_masked_pixels_per_frame(self) -> None:
        video = np.asarray(
            [
                [[1.0, 3.0], [5.0, 7.0]],
                [[2.0, 4.0], [6.0, 8.0]],
            ],
            dtype=np.float32,
        )
        mask = np.asarray([[True, False], [False, True]])

        np.testing.assert_allclose(masked_video_signal(video, mask), [4.0, 5.0])

    def test_velocity_masked_signal_ignores_nans_outside_and_inside_mask(self) -> None:
        velocity = np.asarray(
            [
                [[np.nan, 1.0], [3.0, np.nan]],
                [[np.nan, np.nan], [5.0, np.nan]],
            ],
            dtype=np.float32,
        )
        mask = np.asarray([[False, True], [True, False]])

        np.testing.assert_allclose(_velocity_masked_signal(velocity, mask), [2.0, 5.0])

    def test_velocity_gradient_uses_averaged_map_not_3d_peak(self) -> None:
        velocity_map = np.asarray(
            [
                [[0.0, 50.0, 70.0]],
                [[0.0, 50.0, 70.0]],
                [[90.0, 50.0, 70.0]],
            ],
            dtype=np.float32,
        )
        velocity_avg = np.mean(velocity_map, axis=0)
        section_mask = np.ones_like(velocity_avg, dtype=bool)

        values = _velocity_gradient_values(velocity_avg, section_mask, velocity_map)

        np.testing.assert_allclose(values, [[30.0 / 70.0, 50.0 / 70.0, 1.0]])
        self.assertEqual(70.0, _velocity_colorbar_vmax(velocity_avg, section_mask, velocity_map))

    @unittest.skipUnless(_has_matplotlib(), "matplotlib is not installed")
    def test_export_pulse_pngs_writes_non_empty_pngs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output = FakeOutput(Path(temp_dir))
            context = _synthetic_context()
            written = export_pulse_pngs(output, context, SimpleNamespace(), log=None)

            self.assertGreaterEqual(len(written), 35)
            for required in (
                "sample_f_artery_graph.png",
                "sample_find_systoles_indices_artery.png",
                "sample_ArterialSpectralAnalysis_v_artery.png",
                "sample_v_map.png",
                "sample_AVGflowVideoCombined.png",
            ):
                path = Path(temp_dir) / "png" / required
                self.assertTrue(path.is_file(), required)
                self.assertGreater(path.stat().st_size, 0, required)


def _synthetic_context():
    frames = 32
    height = 8
    width = 8
    t = np.linspace(0, 4 * np.pi, frames, dtype=np.float32)
    artery_signal = 18.0 + 3.0 * np.sin(t)
    vein_signal = 12.0 + 2.0 * np.sin(t + 0.8)
    artery_mask = np.zeros((height, width), dtype=bool)
    vein_mask = np.zeros((height, width), dtype=bool)
    artery_mask[2:5, 2:4] = True
    vein_mask[3:6, 5:7] = True
    section_mask = np.ones((height, width), dtype=bool)
    base_map = np.linspace(0.0, 1.0, height * width, dtype=np.float32).reshape(height, width)
    f_video = np.stack([base_map + 0.05 * i for i in range(frames)]).astype(np.float32)
    f_bkg = (f_video * 0.5).astype(np.float32)
    delta = (f_video - f_bkg).astype(np.float32)
    velocity = np.stack(
        [
            np.where(artery_mask, artery_signal[i], 0.0)
            + np.where(vein_mask, vein_signal[i], 0.0)
            for i in range(frames)
        ]
    ).astype(np.float32)
    beat_indices = np.asarray([2, 10, 18, 26], dtype=np.int32)
    analysis = {
        "fRMS": f_video,
        "fRMS_bkg": f_bkg,
        "fRMS_avg": np.mean(f_video, axis=0),
        "fRMS_bkg_avg": np.mean(f_bkg, axis=0),
        "deltafRMS": delta,
        "velocity_section_mask": section_mask,
        "retinal_vessel_velocity": velocity,
        "velocity_map_avg": np.mean(velocity, axis=0),
        "retinal_artery_velocity_signal": artery_signal,
        "retinal_vein_velocity_signal": vein_signal,
        "retinal_artery_velocity_signal_filtered": artery_signal,
        "retinal_vein_velocity_signal_filtered": vein_signal,
        "retinal_artery_velocity_signal_derivative": np.gradient(artery_signal).astype(np.float32),
        "retinal_vein_velocity_signal_derivative": np.gradient(vein_signal).astype(np.float32),
        "beat_indices": beat_indices,
    }
    source_data = SimpleNamespace(
        timing=SimpleNamespace(dt_seconds=0.1),
        moment0=np.stack([base_map + 1.0 for _ in range(frames)]).astype(np.float32),
        retinal_artery_mask=artery_mask,
        retinal_vein_mask=vein_mask,
    )
    return SimpleNamespace(source_data=source_data, dopplerview_analysis=analysis)


if __name__ == "__main__":
    unittest.main()
