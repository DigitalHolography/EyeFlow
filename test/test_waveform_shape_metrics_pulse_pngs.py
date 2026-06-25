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

from input_output.output_manager import OutputType  # noqa: E402
from pipelines.waveform_shape_metrics.velocity.figures import (  # noqa: E402
    PULSE_PNG_SUFFIXES,
    export_pulse_pngs,
)


class FakeOutput:
    available = True

    def __init__(self, root: Path, stem: str = "sample") -> None:
        self.root = root
        self.manager = SimpleNamespace(layout=SimpleNamespace(stem=stem))

    def path_for(self, output_type: OutputType, filename: str | None = None) -> Path:
        assert output_type is OutputType.PNG
        return self.root / "png" / (filename or "sample")


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
        self.assertIn("AVGflowVideoCombined.png", PULSE_PNG_SUFFIXES)

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
