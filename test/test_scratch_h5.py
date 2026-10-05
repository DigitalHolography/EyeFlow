"""Tests for the shared RAM-backed HDF5 workspace."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.writers.h5 import scratch_h5  # noqa: E402
from pipelines.heartbeat_core.scratch import heartbeat_scratch_h5  # noqa: E402
from pipelines.waveform_velocity_core.scratch import velocity_scratch_h5  # noqa: E402


class ScratchH5Tests(unittest.TestCase):
    def test_shared_context_is_nonpersistent_and_accepts_purpose(self) -> None:
        with scratch_h5(purpose="diagnostic", filename_prefix="eyeflow-test") as h5file:
            h5file.create_dataset("value", data=[1, 2, 3])
            self.assertEqual("diagnostic", h5file.attrs["purpose"])
            self.assertEqual("memory", h5file.attrs["storage"])
            self.assertTrue(h5file.attrs["temporary"])
            filename = Path(h5file.filename)
            self.assertFalse(filename.exists())
        self.assertFalse(filename.exists())

    def test_pipeline_wrappers_use_shared_context(self) -> None:
        for factory, purpose in (
            (heartbeat_scratch_h5, "EyeFlow heartbeat intermediates"),
            (velocity_scratch_h5, "EyeFlow retinal velocity intermediates"),
        ):
            with self.subTest(purpose=purpose), factory(None) as h5file:
                self.assertEqual(purpose, h5file.attrs["purpose"])


if __name__ == "__main__":
    unittest.main()
