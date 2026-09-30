import unittest
from contextlib import nullcontext
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

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
from pipelines.retinal_velocity.models import (
    RetinalVelocity,
    RetinalVelocityData,
    RetinalVelocityMaps,
    VesselVelocity,
    VesselVelocitySignals,
)
from pipelines.retinal_velocity.runner import (
    cardiac_cycle_indexes,
    retinal_velocity,
    run_retinal_velocity,
)


class RetinalVelocityTests(unittest.TestCase):
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

    def test_pipeline_retains_velocity_map_for_waveform_consumers(self):
        source = _source_data()
        analysis = _cardiac_cycle()
        ctx = SimpleNamespace(
            state=PipelineState(),
            pipeline_scheduled=lambda name: name == "waveform_velocity",
        )
        estimated = {}

        def estimator(**kwargs):
            video = kwargs["velocity_video_output"]
            video.fill(np.float32(7.0))
            data = _retinal_velocity_data(video)
            estimated["data"] = data
            return data

        with (
            patch(
                "pipelines.retinal_velocity.runner.load_retinal_velocity_inputs",
                return_value=source,
            ),
            patch(
                "pipelines.retinal_velocity.runner.retinal_velocity_scratch_h5",
                return_value=nullcontext(object()),
            ),
            patch(
                "pipelines.retinal_velocity.runner.estimate_retinal_velocity",
                side_effect=estimator,
            ),
            patch(
                "pipelines.retinal_velocity.runner.detect_cardiac_cycles",
                return_value=(analysis, "artery"),
            ),
            patch(
                "pipelines.retinal_velocity.signal_processing._filter",
                side_effect=lambda signal, *_: np.asarray(signal, dtype=np.float32),
            ),
        ):
            result, outputs = run_retinal_velocity(ctx)

        self.assertIs(result, retinal_velocity(ctx))
        np.testing.assert_array_equal(cardiac_cycle_indexes(ctx), [0, 3, 5])
        self.assertIsInstance(result, RetinalVelocity)
        self.assertIs(result.maps, estimated["data"].maps)
        self.assertIs(result.artery.signals, estimated["data"].artery)
        self.assertIs(result.vein.signals, estimated["data"].vein)
        self.assertTrue(result.has_velocity_map)
        np.testing.assert_array_equal(result.maps.velocity, 7.0)
        self.assertIn("Processing/Maps/VelocityAverage/value", outputs)
        self.assertIn("Processing/Maps/DeltaFRMSAverage/value", outputs)


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
        velocity_average=image,
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
