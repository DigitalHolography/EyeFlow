from contextlib import nullcontext
from dataclasses import fields
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from pipeline_engine.context import PipelineState
from pipelines.heartbeat_core.runner import (
    HeartbeatResult,
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


if __name__ == "__main__":
    unittest.main()
