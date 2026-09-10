import unittest
from contextlib import nullcontext
from dataclasses import fields
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from pipeline_engine.context import PipelineState
from pipelines.heartbeat_core.runner import (
    HeartbeatResult,
    cached_velocity_estimation,
    heartbeat_result,
    run_heartbeat_core,
)


class HeartbeatCoreTests(unittest.TestCase):
    def test_exposes_only_cycle_boundaries_and_index_base(self):
        inputs = SimpleNamespace(
            moment0=object(),
            moment2=object(),
            artery_mask=np.ones((3, 3), dtype=bool),
            vein_mask=np.ones((3, 3), dtype=bool),
            optic_disc_center=None,
            optic_disc_width=None,
            optic_disc_height=None,
            timing=SimpleNamespace(dt_seconds=0.02),
            local_background_dist=4,
            index_base=0,
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
                    "retinal_artery_velocity_signal": np.arange(15, dtype=np.float32)
                },
            ),
            patch(
                "pipelines.heartbeat_core.runner.run_heartbeat_analysis",
                return_value=analysis,
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
        inputs = SimpleNamespace(
            moment0=moment0,
            moment2=moment2,
            artery_mask=artery,
            vein_mask=vein,
            optic_disc_center=None,
            optic_disc_width=None,
            optic_disc_height=None,
            timing=SimpleNamespace(dt_seconds=0.02),
            local_background_dist=1,
            index_base=0,
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
                "pipelines.heartbeat_core.runner.run_heartbeat_analysis",
                return_value=analysis,
            ),
        ):
            run_heartbeat_core(ctx)

        waveform_source = SimpleNamespace(
            moment0=moment0,
            moment2=moment2,
            retinal_artery_mask=artery.copy(),
            retinal_vein_mask=vein.copy(),
            optic_disc_center=None,
            optic_disc_width=None,
            optic_disc_height=None,
            local_background_dist=1,
        )
        reused = cached_velocity_estimation(ctx, waveform_source)

        self.assertIsNotNone(reused)
        self.assertIs(reused["velocity_map"], ctx.state.get(
            "heartbeat.cache.velocity_estimation"
        ).analysis["velocity_map"])
        np.testing.assert_array_equal(reused["velocity_map"], 7.0)

        waveform_source.retinal_artery_mask[0, 0] = True
        self.assertIsNone(cached_velocity_estimation(ctx, waveform_source))


if __name__ == "__main__":
    unittest.main()
