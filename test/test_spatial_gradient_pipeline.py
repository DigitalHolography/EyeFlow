"""Tests for the independently selectable spatial-gradient pipeline."""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import pipelines
from calculations.math.spatial_gradient import sobel_spatial_gradient
from input_output.schema import PixelPitch
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

        gradient = sobel_spatial_gradient(ramp[None])[0]

        np.testing.assert_allclose(gradient[:, 1:-1], 8.0)
        np.testing.assert_allclose(gradient[:, (0, -1)], 4.0)

    def test_pipeline_is_optional_and_not_a_waveform_dependency(self) -> None:
        pipelines.load_pipeline_catalog()

        self.assertIn("spatial_gradient_moment0", PIPELINE_REGISTRY)
        dag = PipelineDAG(PIPELINE_REGISTRY.values())
        gradient_plan = dag.resolve_targets(["spatial_gradient_moment0"])
        self.assertLess(
            gradient_plan.names.index("velocity"),
            gradient_plan.names.index("spatial_gradient_moment0"),
        )
        self.assertLess(
            gradient_plan.names.index("topology_core"),
            gradient_plan.names.index("spatial_gradient_moment0"),
        )
        self.assertNotIn("velocity_analysis", gradient_plan.names)

        waveform_plan = dag.resolve_targets(["velocity_analysis"])
        self.assertNotIn("spatial_gradient_moment0", waveform_plan.names)

    def test_displacement_map_has_no_cardiac_cycle_or_waveform_dependency(self) -> None:
        pipelines.load_pipeline_catalog()

        dag = PipelineDAG(PIPELINE_REGISTRY.values())
        plan = dag.resolve_targets(["displacement_map"])

        self.assertEqual(("displacement_map",), plan.names)

        profile_plan = dag.resolve_targets(["displacement_profiles"])
        self.assertIn("velocity", profile_plan.names)
        self.assertIn("topology_core", profile_plan.names)
        self.assertLess(
            profile_plan.names.index("displacement_map"),
            profile_plan.names.index("displacement_profiles"),
        )
        self.assertLess(
            profile_plan.names.index("topology_core"),
            profile_plan.names.index("displacement_profiles"),
        )

    def test_runner_publishes_outputs_only_when_directly_targeted(self) -> None:
        import tempfile
        from pathlib import Path
        import h5py
        from pipeline_engine import PipelineDescriptor, run_pipelines_to_output
        from input_output.runs.layout import HoloRunLayout
        from input_output.output_manager import OutputManager

        inputs = SimpleNamespace(holodoppler=SimpleNamespace(
            pixel_pitch=PixelPitch(20e-6, 20e-6), timing=SimpleNamespace(dt_seconds=0.1)))
        for targeted in (True, False):
            with self.subTest(targeted=targeted), tempfile.TemporaryDirectory() as directory:
                observed = []
                def observe(ctx):
                    observed.append(ctx.state[gradient_runner.SPATIAL_GRADIENT_PRODUCTS_STATE])
                source = PipelineDescriptor("spatial_gradient_moment0", "gradient", True,
                    pipeline_factory=lambda: SimpleNamespace(name="spatial_gradient_moment0", run=gradient_runner.run_spatial_gradient_moment0))
                consumer = PipelineDescriptor("consumer", "consumer", True, dag_requires=(source.name,),
                    pipeline_factory=lambda: SimpleNamespace(name="consumer", run=observe))
                manager = OutputManager(HoloRunLayout.from_holo(Path(directory) / "scan.holo"))
                with (patch.object(gradient_runner, "load_vessel_topology_inputs", return_value=inputs),
                      patch.object(gradient_runner, "cardiac_cycle_indexes", return_value=np.array([0,3])),
                      patch.object(gradient_runner, "prepared_topologies", return_value={}),
                      patch.object(gradient_runner, "extract_spatial_gradient_segments", return_value=(None,None)) as extract,
                      patch.object(gradient_runner, "pack_spatial_gradient_profile_outputs", return_value={"gradient": 1})):
                    output = run_pipelines_to_output(output_manager=manager, pipelines=[source,consumer],
                        target_names=[source.name if targeted else consumer.name],
                        holodoppler_h5=None, doppler_vision_h5=None, copy_hd_passthrough=False)
                extract.assert_called_once()
                self.assertEqual({"gradient": 1}, observed[0].outputs)
                with h5py.File(output) as file:
                    self.assertEqual(targeted, "gradient" in file)


if __name__ == "__main__":
    unittest.main()
