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
from pipeline_engine.context import PipelineState
from pipelines.velocity.cardiac_cycle import detect_source_cardiac_cycles
from pipelines.velocity.models import (
    RetinalVelocity,
    RetinalVelocityData,
    RetinalVelocityMaps,
    VesselVelocity,
    VesselVelocitySignals,
)
from pipelines.velocity.runner import (
    _log_cardiac_cycle_warnings,
    cardiac_cycle_indexes,
    run_velocity,
    velocity,
)
from pipelines.velocity.signal_processing import build_velocity


class RetinalVelocityTests(unittest.TestCase):
    def test_moments_detector_receives_raw_frequency(self):
        with patch(
            "pipelines.velocity.cardiac_cycle.cardiac_cycles_from_available_vessel",
            return_value=(_cardiac_cycle(), "artery"),
        ) as detector:
            detect_source_cardiac_cycles(_source_data())
        artery_frequency, vein_frequency = detector.call_args.args
        np.testing.assert_allclose(artery_frequency, np.sqrt(2.0))
        np.testing.assert_allclose(vein_frequency, np.sqrt(2.0))

    def test_bands_detector_receives_raw_calibrated_ratio(self):
        source = _source_data()
        low = np.full(source.image_maps.moment0.shape, 2.0, dtype=np.float32)
        high = np.full_like(low, 6.0)
        low[0] = 0.0
        source = replace(
            source,
            image_maps=ImageMaps(None, None, band_lf=low, band_hf=high),
            velocity_estimation_method="frequency_bands",
            band_ratio_frequency_scale_hz=7.0,
        )
        with patch(
            "pipelines.velocity.cardiac_cycle.cardiac_cycles_from_available_vessel",
            return_value=(_cardiac_cycle(), "artery"),
        ) as detector:
            detect_source_cardiac_cycles(source)
        artery_frequency, vein_frequency = detector.call_args.args
        np.testing.assert_array_equal(artery_frequency, [0.0, 21.0, 21.0, 21.0, 21.0, 21.0])
        np.testing.assert_array_equal(vein_frequency, artery_frequency)

    def test_velocity_filter_does_not_reuse_frequency_detector_signal(self):
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

    def test_pipeline_retains_velocity_map_for_analysis_consumers(self):
        source = _source_data()
        analysis = _cardiac_cycle()
        ctx = SimpleNamespace(
            state=PipelineState(),
            pipeline_scheduled=lambda name: name == "velocity_analysis",
        )
        estimated = {}

        def load_inputs(ctx):
            if ctx.velocity_estimation_method == "doppler_moments":
                return source
            raise KeyError("missing bands")

        def estimator(**kwargs):
            estimated["image_maps"] = kwargs["image_maps"]
            self.assertTrue(kwargs["retain_velocity_video"])
            video = np.full(source.image_maps.moment0.shape, 7.0, dtype=np.float32)
            data = _retinal_velocity_data(video)
            estimated["data"] = data
            return data

        with (
            patch(
                "pipelines.velocity.runner.load_velocity_inputs",
                side_effect=load_inputs,
            ),
            patch(
                "pipelines.velocity.runner.estimate_retinal_velocity",
                side_effect=estimator,
            ),
            patch(
                "pipelines.velocity.runner.detect_source_cardiac_cycles",
                return_value=(analysis, "artery"),
            ),
            patch(
                "pipelines.velocity.signal_processing._filter",
                side_effect=lambda signal, *_: np.asarray(signal, dtype=np.float32),
            ),
        ):
            result, outputs = run_velocity(ctx)

        self.assertIs(result, velocity(ctx))
        self.assertIs(estimated["image_maps"], source.image_maps)
        np.testing.assert_array_equal(cardiac_cycle_indexes(ctx), [0, 3, 5])
        self.assertIsInstance(result, RetinalVelocity)
        self.assertIs(result.maps, estimated["data"].maps)
        self.assertIs(result.artery.signals, estimated["data"].artery)
        self.assertIs(result.vein.signals, estimated["data"].vein)
        self.assertTrue(result.has_velocity_map)
        np.testing.assert_array_equal(result.maps.velocity, 7.0)
        self.assertIn("Processing/Maps/VelocityAverageMasked/value", outputs)

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
