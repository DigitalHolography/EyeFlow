import unittest
from contextlib import nullcontext
from dataclasses import fields
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from calculations.topology import OpticDisc
from input_output.schema import (
    DopplerViewMetadata,
    ImageMaps,
    PixelPitch,
    RetinalSegmentation,
    RetinalSourceData,
    VesselMasks,
)
from input_output.schema.source_data import HolodopplerMetadata, HolodopplerTiming
from pipeline_engine.context import PipelineState
from pipelines.heartbeat_core.runner import (
    HeartbeatResult,
    cached_velocity_estimation,
    heartbeat_result,
    run_heartbeat_core,
)


class HeartbeatCoreTests(unittest.TestCase):
    def test_exposes_only_cycle_boundaries_and_index_base(self):
        inputs = _source_data(
            object(),
            object(),
            np.ones((3, 3), dtype=bool),
            np.ones((3, 3), dtype=bool),
            OpticDisc(
                np.zeros((3, 3), dtype=bool), (1.0, 1.0), None, None
            ),
            dt_seconds=0.02,
            local_background_dist=4,
        )
        analysis = SimpleNamespace(
            systole=SimpleNamespace(systole_indexes=np.asarray([2, 7, 12]))
        )
        ctx = SimpleNamespace(state=PipelineState())

        with (
            patch(
                "pipelines.heartbeat_core.runner.load_heartbeat_inputs",
                return_value=inputs,
            ),
            patch(
                "pipelines.heartbeat_core.runner.heartbeat_scratch_h5",
                return_value=nullcontext(object()),
            ),
            patch(
                "pipelines.heartbeat_core.runner.run_chunked_velocity_estimator",
                return_value={
                    "retinal_artery_velocity_signal": np.arange(15, dtype=np.float32),
                    "retinal_vein_velocity_signal": np.arange(15, dtype=np.float32),
                },
            ),
            patch(
                "pipelines.heartbeat_core.runner.heartbeat_from_available_vessel",
                return_value=(analysis, "artery"),
            ),
        ):
            result = run_heartbeat_core(ctx)

        self.assertEqual(
            [field.name for field in fields(HeartbeatResult)],
            ["cycle_boundary_indexes", "index_base"],
        )
        np.testing.assert_array_equal(result.cycle_boundary_indexes, [2, 7, 12])
        self.assertEqual(0, result.index_base)
        self.assertIs(result, heartbeat_result(ctx))

    def test_retains_reusable_velocity_video_for_waveform_core(self):
        moment0 = np.ones((4, 3, 3), dtype=np.float32)
        moment2 = np.full_like(moment0, 2.0)
        artery = np.zeros((3, 3), dtype=bool)
        vein = np.zeros_like(artery)
        artery[1, 0] = True
        vein[1, 2] = True
        inputs = _source_data(
            moment0,
            moment2,
            artery,
            vein,
            OpticDisc(
                np.zeros((3, 3), dtype=bool), (1.0, 1.0), None, None
            ),
            dt_seconds=0.02,
            local_background_dist=1,
        )
        analysis = SimpleNamespace(
            systole=SimpleNamespace(systole_indexes=np.asarray([0, 3]))
        )
        ctx = SimpleNamespace(
            state=PipelineState(),
            pipeline_scheduled=lambda name: name == "waveform_velocity_core",
        )

        def estimator(**kwargs):
            video = kwargs["velocity_video_output"]
            video.fill(np.float32(7.0))
            return {
                "velocity_map": video,
                "retinal_artery_velocity_signal": np.arange(4, dtype=np.float32),
                "retinal_vein_velocity_signal": np.arange(4, dtype=np.float32),
            }

        with (
            patch(
                "pipelines.heartbeat_core.runner.load_heartbeat_inputs",
                return_value=inputs,
            ),
            patch(
                "pipelines.heartbeat_core.runner.heartbeat_scratch_h5",
                return_value=nullcontext(SimpleNamespace()),
            ),
            patch(
                "pipelines.heartbeat_core.runner.run_chunked_velocity_estimator",
                side_effect=estimator,
            ),
            patch(
                "pipelines.heartbeat_core.runner.heartbeat_from_available_vessel",
                return_value=(analysis, "artery"),
            ),
        ):
            run_heartbeat_core(ctx)

        waveform_source = inputs
        reused = cached_velocity_estimation(ctx, waveform_source)

        self.assertIsNotNone(reused)
        self.assertIs(reused["velocity_map"], ctx.state.get(
            "heartbeat.cache.velocity_estimation"
        ).analysis["velocity_map"])
        np.testing.assert_array_equal(reused["velocity_map"], 7.0)

        waveform_source.segmentation.vessels.artery[0, 0] = True
        self.assertIsNone(cached_velocity_estimation(ctx, waveform_source))


def _source_data(
    moment0,
    moment2,
    artery,
    vein,
    optic_disc,
    *,
    dt_seconds: float,
    local_background_dist: int,
) -> RetinalSourceData:
    return RetinalSourceData(
        image_maps=ImageMaps(moment0, moment2),
        segmentation=RetinalSegmentation(
            VesselMasks(artery, vein),
            optic_disc,
        ),
        holodoppler=HolodopplerMetadata(
            HolodopplerTiming(1.0 / dt_seconds, 1.0),
            PixelPitch(20e-6, 20e-6),
        ),
        doppler_view=DopplerViewMetadata(local_background_dist, False),
    )


if __name__ == "__main__":
    unittest.main()
