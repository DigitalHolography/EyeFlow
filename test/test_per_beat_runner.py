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
from calculations.math import band_limited_ifft_abs  # noqa: E402
from pipelines.waveform_velocity_core.per_beat_outputs import (  # noqa: E402
    pack_velocity_per_beat_outputs,
)
from pipelines.waveform_velocity_core.runner import (  # noqa: E402
    _filtered_velocity_signals_for_per_beat,
)


class PerBeatRunnerTests(unittest.TestCase):
    def test_band_limited_reconstruction_matches_matlab_abs_ifft(self) -> None:
        phase = np.linspace(
            0.0,
            2.0 * np.pi,
            64,
            endpoint=False,
            dtype=np.float32,
        )
        waveform = (
            100.0
            + 30.0 * np.sin(phase)
            + 8.0 * np.cos(2.0 * phase)
        ).astype(np.float32)

        spectrum = np.fft.fft(waveform)
        reconstructed = band_limited_ifft_abs(
            spectrum,
            waveform.size,
            harmonic_count=3,
        )

        matlab_band_limited_spectrum = spectrum[:3] * 2.0
        matlab_band_limited_spectrum[0] = spectrum[0]
        expected = np.abs(
            np.fft.ifft(matlab_band_limited_spectrum, n=waveform.size)
        )
        np.testing.assert_allclose(reconstructed, expected, rtol=1e-6, atol=1e-5)

    def test_per_beat_input_explicitly_uses_filtered_vessel_signals(self) -> None:
        analysis = {
            "retinal_artery_velocity_signal": np.asarray([100.0, 200.0]),
            "retinal_vein_velocity_signal": np.asarray([300.0, 400.0]),
            "retinal_artery_velocity_signal_filtered": np.asarray([1.0, 2.0]),
            "retinal_vein_velocity_signal_filtered": np.asarray([3.0, 4.0]),
        }

        artery, vein = _filtered_velocity_signals_for_per_beat(analysis)

        np.testing.assert_array_equal(artery, [1.0, 2.0])
        np.testing.assert_array_equal(vein, [3.0, 4.0])
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
        self.assertIs(result.heartbeat, heartbeat)
        np.testing.assert_array_equal(
            result.cycle_boundary_indexes,
            inputs.cycle_boundary_indexes,
        )
        outputs = pack_velocity_per_beat_outputs(result)
        self.assertFalse(any("Vmax" in path for path in outputs))
        self.assertFalse(any("Vmin" in path for path in outputs))
        self.assertFalse(any("VTI" in path for path in outputs))


if __name__ == "__main__":
    unittest.main()
