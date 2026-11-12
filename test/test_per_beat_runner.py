"""Tests for per-beat period estimation."""

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

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (  # noqa: E402
    spectral_cardiac_cycle_analysis,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat.runner import (  # noqa: E402
    PerBeatAnalysisInput,
    run_per_beat_analysis,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat.signal import (  # noqa: E402
    per_beat_signal_analysis,
)
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from pipelines.velocity_analysis.builder import (  # noqa: E402
    _raw_velocity_signals_for_per_beat,
    _run_velocity_per_beat_analysis,
)
from pipelines.velocity_analysis.outputs.per_beat import (  # noqa: E402
    pack_velocity_per_beat_inputs,
)


class PerBeatRunnerTests(unittest.TestCase):
    def test_waveform_analysis_reuses_canonical_segment_signals(self) -> None:
        artery_signal = np.ones((2, 3, 8), dtype=np.float32)
        vein_signal = np.full((2, 3, 8), 2.0, dtype=np.float32)
        artery_segments = SimpleNamespace(
            profile=SimpleNamespace(segment_signal=artery_signal)
        )
        vein_segments = SimpleNamespace(
            profile=SimpleNamespace(segment_signal=vein_signal)
        )
        retinal = SimpleNamespace(
            cycle_boundary_indexes=np.asarray([0, 4, 7], dtype=np.int32),
            cardiac_cycle=SimpleNamespace(spectral="cardiac-cycle"),
            continuous=lambda vessel, raw=False: np.arange(8, dtype=np.float32),
        )
        source_data = SimpleNamespace(
            source=SimpleNamespace(
                holodoppler=SimpleNamespace(
                    timing=SimpleNamespace(dt_seconds=0.1)
                )
            ),
            provenance={"beat_index_base": 0},
        )

        with patch(
            "pipelines.velocity_analysis.builder.run_per_beat_analysis",
            return_value="result",
        ) as run:
            result = _run_velocity_per_beat_analysis(
                retinal,
                source_data,
                artery_segments,
                vein_segments,
                harmonic_count=4,
            )

        self.assertEqual("result", result)
        inputs = run.call_args.args[0]
        self.assertIs(inputs.arterial_velocity_segments, artery_signal)
        self.assertIsNotNone(inputs.arterial_filtered_velocity_segments)
        self.assertIs(inputs.venous_velocity_segments, vein_signal)
        self.assertIsNotNone(inputs.venous_filtered_velocity_segments)

    def test_per_beat_input_explicitly_uses_raw_global_vessel_signals(self) -> None:
        raw_values = {
            "artery": np.asarray([100.0, 200.0], dtype=np.float32),
            "vein": np.asarray([300.0, 400.0], dtype=np.float32),
        }
        analysis = SimpleNamespace(
            continuous=lambda vessel, raw=False: raw_values[vessel]
        )

        artery, vein = _raw_velocity_signals_for_per_beat(analysis)

        np.testing.assert_array_equal(artery, [100.0, 200.0])
        np.testing.assert_array_equal(vein, [300.0, 400.0])
        self.assertEqual(artery.dtype, np.float32)
        self.assertEqual(vein.dtype, np.float32)

    def test_beat_period_uses_each_matlab_systole_interval(self) -> None:
        dt_seconds = 0.05
        time = np.arange(200, dtype=np.float32) * dt_seconds
        signal = (
            1.0
            + np.sin(2.0 * np.pi * 1.0 * time)
            + 0.45 * np.sin(2.0 * np.pi * 2.0 * time)
        ).astype(np.float32)
        cardiac_cycle = spectral_cardiac_cycle_analysis(
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
            band_limited_harmonics=4,
            cardiac_cycle=cardiac_cycle,
            dt_seconds=dt_seconds,
            index_base=0,
        )

        result = run_per_beat_analysis(inputs)

        np.testing.assert_allclose(result.beat_period_seconds, [4.0, 3.0, 2.95])
        np.testing.assert_allclose(
            result.artery.vti_per_beat,
            np.sum(result.artery.signal.velocity_signal_per_beat, axis=1)
            * dt_seconds,
        )
        self.assertIs(result.cardiac_cycle, cardiac_cycle)
        np.testing.assert_array_equal(
            result.cycle_boundary_indexes,
            inputs.cycle_boundary_indexes,
        )
        inputs = pack_velocity_per_beat_inputs(
            result,
            velocity=SimpleNamespace(
                cycle_durations_seconds=np.asarray(
                    [4.0, 3.0, 2.95],
                    dtype=np.float32,
                ),
                provenance={},
            ),
        )
        schema = EyeFlowOutputPaths.active()
        cycle_durations, attrs = inputs[
            schema.cardiac_cycle.systolic_cycle_duration_seconds
        ]
        np.testing.assert_allclose(cycle_durations, [4.0, 3.0, 2.95])
        self.assertEqual("s", attrs["unit"])


if __name__ == "__main__":
    unittest.main()
