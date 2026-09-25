"""Tests for the independently selectable spatial-gradient pipeline."""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import pipelines
from calculations.math.spatial_gradient import spatial_gradient
from calculations.topology import OpticDisc
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines.spatial_gradient_moment0 import runner as gradient_runner


class _State:
    def __init__(self, values=None):
        self.values = dict(values or {})

    def get(self, key, default=None):
        return self.values.get(key, default)

    def set(self, key, value):
        self.values[key] = value


class SpatialGradientPipelineTests(unittest.TestCase):
    def test_spatial_gradient_uses_sobel_magnitude(self) -> None:
        ramp = np.tile(np.arange(5, dtype=np.float32), (4, 1))

        gradient = spatial_gradient(ramp)

        np.testing.assert_allclose(gradient[:, 1:-1], 8.0)
        np.testing.assert_allclose(gradient[:, (0, -1)], 4.0)

    def test_pipeline_is_optional_and_not_a_waveform_dependency(self) -> None:
        pipelines.load_pipeline_catalog()

        self.assertIn("spatial_gradient_moment0", PIPELINE_REGISTRY)
        dag = PipelineDAG(PIPELINE_REGISTRY.values())
        gradient_plan = dag.resolve_targets(["spatial_gradient_moment0"])
        self.assertLess(
            gradient_plan.names.index("heartbeat_core"),
            gradient_plan.names.index("spatial_gradient_moment0"),
        )
        self.assertLess(
            gradient_plan.names.index("topology_core"),
            gradient_plan.names.index("spatial_gradient_moment0"),
        )
        self.assertNotIn("waveform_velocity_core", gradient_plan.names)
        self.assertNotIn("waveform_velocity", gradient_plan.names)

        waveform_plan = dag.resolve_targets(["waveform_velocity"])
        self.assertNotIn("spatial_gradient_moment0", waveform_plan.names)

    def test_displacement_map_has_no_heartbeat_or_waveform_dependency(self) -> None:
        pipelines.load_pipeline_catalog()

        plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
            ["displacement_map"]
        )

        self.assertEqual(("displacement_map",), plan.names)

    def test_runner_owns_only_gradient_and_lumen_outputs(self) -> None:
        gradient_segments = SimpleNamespace(
            labels=np.asarray([1]),
            branch_ids=np.asarray([7]),
            transverse_profiles_unmasked=np.ones(
                (2, 1, 3), dtype=np.float32
            ),
        )
        inputs = SimpleNamespace(
            optic_disc=OpticDisc(None, (4.0, 4.0), 2.0, 2.0),
            timing=SimpleNamespace(dt_seconds=0.1),
        )
        topology = {"artery": object(), "vein": object()}
        state = _State()
        ctx = SimpleNamespace(
            state=state,
            output=SimpleNamespace(available=False),
        )

        with (
            patch.object(
                gradient_runner,
                "load_heartbeat_inputs",
                return_value=inputs,
            ),
            patch.object(
                gradient_runner,
                "heartbeat_result",
                return_value=SimpleNamespace(
                    cycle_boundary_indexes=(0, 3),
                    index_base=0,
                ),
            ),
            patch.object(
                gradient_runner,
                "prepared_topologies",
                return_value=topology,
            ),
            patch.object(
                gradient_runner,
                "extract_spatial_gradient_segments",
                return_value=(gradient_segments, gradient_segments),
            ) as extract,
            patch.object(
                gradient_runner,
                "pack_spatial_gradient_profile_outputs",
                return_value={"gradient": 1},
            ),
        ):
            outputs = gradient_runner.run_spatial_gradient_moment0(ctx)

        self.assertEqual({"gradient": 1}, outputs)
        self.assertIs(extract.call_args.args[2], topology)
        products = state.get(gradient_runner.SPATIAL_GRADIENT_PRODUCTS_STATE)
        self.assertEqual(outputs, products.outputs)


if __name__ == "__main__":
    unittest.main()
