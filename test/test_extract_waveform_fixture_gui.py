"""Tests for the executable's non-visual extraction job."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


TOOLS_DIR = Path(__file__).resolve().parents[1] / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))
SCRIPT_PATH = TOOLS_DIR / "extract_waveform_fixture_gui.py"
SPEC = importlib.util.spec_from_file_location("extract_waveform_fixture_gui", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load {SCRIPT_PATH}")
gui = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = gui
SPEC.loader.exec_module(gui)


class WaveformFixtureGuiJobTests(unittest.TestCase):
    def test_validation_does_not_write_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            source_path = root / "source.h5"
            output_dir = root / "fixtures"
            _write_minimal_source(source_path)
            messages: list[str] = []

            summary = gui.run_extraction_job(
                source_path,
                output_dir,
                recursive=False,
                dry_run=True,
                log=messages.append,
            )

            self.assertEqual(summary.succeeded, 1)
            self.assertEqual(summary.failed, 0)
            self.assertFalse(output_dir.exists())
            self.assertTrue(any("compatible" in message for message in messages))

    def test_extraction_job_creates_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            source_path = root / "source.h5"
            output_dir = root / "fixtures"
            _write_minimal_source(source_path)

            summary = gui.run_extraction_job(
                source_path,
                output_dir,
                recursive=False,
                dry_run=False,
                log=lambda _message: None,
            )

            self.assertEqual(summary.created, 1)
            self.assertEqual(summary.failed, 0)
            self.assertEqual(len(list(output_dir.glob("waveform_fixture_*.h5"))), 1)


def _write_minimal_source(path: Path) -> None:
    with h5py.File(path, "w") as source:
        source.create_dataset(
            "analysis/retinal_artery_velocity_signal",
            data=np.arange(12, dtype=np.float32),
        )
        source.create_dataset(
            "analysis/retinal_vein_velocity_signal",
            data=np.arange(12, dtype=np.float32) + 20,
        )
        source.create_dataset(
            "analysis/beat_indices",
            data=np.asarray([0, 5, 10], dtype=np.int32),
        )
        source.create_dataset(
            "analysis/time_per_beat",
            data=np.asarray([0.5, 0.5], dtype=np.float32),
        )


if __name__ == "__main__":
    unittest.main()
