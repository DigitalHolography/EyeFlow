"""Tests for per-beat period estimation."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity.signal_analysis.per_beat.runner import (  # noqa: E402
    PerBeatAnalysisInput,
    run_per_beat_analysis,
)
from calculations.blood_flow_velocity.signal_analysis.heartbeat import (  # noqa: E402
    spectral_heartbeat_analysis,
)


class PerBeatRunnerTests(unittest.TestCase):
    def test_beat_period_uses_spectral_heartbeat_not_systole_periods(self) -> None:
        dt_seconds = 0.05
        time = np.arange(200, dtype=np.float32) * dt_seconds
        signal = (
            1.0
            + np.sin(2.0 * np.pi * 1.0 * time)
            + 0.45 * np.sin(2.0 * np.pi * 2.0 * time)
        ).astype(np.float32)
        heartbeat = spectral_heartbeat_analysis(
            signal,
            dt_seconds,
            systole_count=4,
        )
        inputs = PerBeatAnalysisInput(
            arterial_velocity_signal=signal,
            venous_velocity_signal=signal,
            cycle_boundary_indexes=np.asarray(
                [0, 80, 140, 199],
                dtype=np.int32,
            ),
            band_limited_signal_harmonic_count=4,
            heartbeat=heartbeat,
            index_base=0,
        )

        result = run_per_beat_analysis(inputs)

        np.testing.assert_allclose(result.beat_period_seconds, [1.0, 1.0, 1.0])
        self.assertIs(result.heartbeat, heartbeat)
        np.testing.assert_array_equal(
            result.cycle_boundary_indexes,
            inputs.cycle_boundary_indexes,
        )


if __name__ == "__main__":
    unittest.main()
