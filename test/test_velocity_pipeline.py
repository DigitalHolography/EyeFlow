import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    SuspectedMissedBeatGap,
)
from calculations.topology import OpticDisc
from input_output.schema import (
    DopplerViewMetadata,
    HolodopplerMetadata,
    HolodopplerTiming,
    ImageMaps,
    PixelPitch,
    RetinalSegmentation,
    RetinalSourceData,
    VesselMasks,
)
from pipelines.velocity.runner import detect_velocity_cardiac_cycles
from pipelines.velocity.models import (
    RetinalVelocity,
    RetinalVelocityData,
    RetinalVelocityMaps,
    VesselVelocity,
    VesselVelocitySignals,
)
from pipelines.velocity.runner import _log_cardiac_cycle_warnings, velocity
from pipelines.velocity.signal_processing import build_velocity


class RetinalVelocityTests(unittest.TestCase):
    def test_detector_receives_each_workflows_unfiltered_signed_physical_velocity(self):
        for method, scale in (("doppler_moments", 1.0), ("frequency_bands", 1e-4)):
            with self.subTest(method=method):
                data = replace(
                    _retinal_velocity_data(np.ones((6, 3, 3), dtype=np.float32)),
                    artery=_vessel_signals((np.arange(6, dtype=np.float32) - 2) * scale),
                    vein=_vessel_signals((3 - np.arange(6, dtype=np.float32)) * scale),
                    provenance={"velocity_estimation_method": method},
                )
                with patch(
                    "pipelines.velocity.runner.cardiac_cycles_from_available_vessel",
                    return_value=(_cardiac_cycle(), "artery"),
                ) as detector:
                    detect_velocity_cardiac_cycles(data, dt_seconds=0.02)
                artery, vein = detector.call_args.args
                self.assertIs(artery, data.artery.velocity)
                self.assertIs(vein, data.vein.velocity)
                self.assertEqual(detector.call_args.kwargs["dt_seconds"], 0.02)

    def test_continuous_velocity_filters_raw_signal_without_double_filtering(self):
        cycle = _cardiac_cycle()
        cycle.systole.signal_filtered = np.full(6, 999.0, dtype=np.float32)
        data = _retinal_velocity_data(np.ones((6, 3, 3), dtype=np.float32))
        with patch(
            "pipelines.velocity.signal_processing._filter",
            side_effect=lambda signal, *_: np.asarray(signal) + 10.0,
        ):
            result = build_velocity(data, cycle, "artery", dt_seconds=0.02)
        np.testing.assert_array_equal(result.continuous("artery"), np.arange(6) + 10.0)
        self.assertIs(result.cardiac_cycle, cycle)

    def test_typed_result_exposes_continuous_and_per_beat_signals(self):
        result = _retinal_velocity()

        np.testing.assert_array_equal(
            result.continuous("artery", raw=True),
            np.arange(6, dtype=np.float32),
        )
        np.testing.assert_array_equal(
            result.continuous("vein"),
            np.arange(6, dtype=np.float32) + 30.0,
        )
        cycles = result.per_beat("artery", raw=True)
        self.assertEqual(2, len(cycles))
        np.testing.assert_array_equal(cycles[0], [0.0, 1.0, 2.0, 3.0])
        np.testing.assert_array_equal(cycles[1], [3.0, 4.0, 5.0])
        np.testing.assert_allclose(result.cycle_durations_seconds, [0.06, 0.04])


    def test_missed_beat_gap_is_recorded_as_a_warning(self):
        gap = SuspectedMissedBeatGap(
            start_index=100,
            stop_index=180,
            interval_samples=80,
            estimated_period_samples=40.0,
            estimated_multiple=2,
        )
        analysis = SimpleNamespace(systole=SimpleNamespace(suspected_missed_beat_gaps=(gap,)))

        with patch("pipelines.velocity.runner.Logger.log_warning") as warning:
            _log_cardiac_cycle_warnings(
                analysis,
                "artery",
                dt_seconds=0.01,
            )

        warning.assert_called_once()
        message = warning.call_args.args[0]
        self.assertIn("frames 100 to 180", message)
        self.assertIn("2x the estimated cardiac period", message)
        self.assertIn("retained", message)


def _retinal_velocity() -> RetinalVelocity:
    return RetinalVelocity(
        maps=_maps(None, shape=(2, 2)),
        artery=VesselVelocity(
            signals=_vessel_signals(np.arange(6, dtype=np.float32)),
            velocity_filtered=np.arange(6, dtype=np.float32) + 10.0,
        ),
        vein=VesselVelocity(
            signals=_vessel_signals(np.arange(6, dtype=np.float32) + 20.0),
            velocity_filtered=np.arange(6, dtype=np.float32) + 30.0,
        ),
        vessel_frms_background=np.ones(6, dtype=np.float32),
        cardiac_cycle=_cardiac_cycle(),
        cardiac_cycle_source="artery",
        dt_seconds=0.02,
    )


def _cardiac_cycle():
    spectral = SimpleNamespace(
        fundamental_hz=1.0,
        heart_rate_bpm=60.0,
        heart_rate_ste_bpm=0.0,
        period_seconds=1.0,
    )
    systole = SimpleNamespace(
        systole_indexes=np.asarray([0, 3, 5], dtype=np.int32),
        signal_filtered=np.arange(6, dtype=np.float32),
        min_peak_distance=2,
        min_peak_height=np.float32(1.0),
    )
    return SimpleNamespace(systole=systole, spectral=spectral)


def _retinal_velocity_data(velocity_map):
    signal = np.arange(velocity_map.shape[0], dtype=np.float32)
    return RetinalVelocityData(
        maps=_maps(velocity_map, shape=(3, 3)),
        artery=_vessel_signals(signal),
        vein=_vessel_signals(signal),
        vessel_frms_background=signal,
    )


def _maps(velocity, *, shape: tuple[int, int]) -> RetinalVelocityMaps:
    image = np.ones(shape, dtype=np.float32)
    return RetinalVelocityMaps(
        velocity=velocity,
        moment0_average=image,
        velocity_average_masked=image,
        frms_average=image,
        frms_background_average=image,
        delta_frms_average=image,
        section_mask=image.astype(bool),
    )


def _vessel_signals(velocity: np.ndarray) -> VesselVelocitySignals:
    values = np.ones(velocity.size, dtype=np.float32)
    return VesselVelocitySignals(
        velocity=velocity,
        frms=values,
        frms_background=values,
        delta_frms=values,
    )


def _source_data() -> RetinalSourceData:
    moment0 = np.ones((6, 3, 3), dtype=np.float32)
    return RetinalSourceData(
        image_maps=ImageMaps(moment0, np.full_like(moment0, 2.0)),
        segmentation=RetinalSegmentation(
            VesselMasks(
                np.ones((3, 3), dtype=bool),
                np.ones((3, 3), dtype=bool),
            ),
            OpticDisc(np.zeros((3, 3), dtype=bool), (1.0, 1.0), None, None),
        ),
        holodoppler=HolodopplerMetadata(
            HolodopplerTiming(50.0, 1.0),
            PixelPitch(20e-6, 20e-6),
        ),
        doppler_view=DopplerViewMetadata(1, False),
    )


if __name__ == "__main__":
    unittest.main()
