"""Tests for the independently selectable spatial-gradient pipeline."""

from app_settings import UserSettings

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from calculations.math.spatial_gradient import sobel_spatial_gradient
from input_output.schema import PixelPitch
from pipelines.spatial_gradient_moment0 import runner as gradient_runner


pytestmark = pytest.mark.usefixtures("inline_output_service")


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
                        holodoppler_h5=None, dopplerview_h5=None, settings=UserSettings(values={'copy_hd_passthrough': False}))
                extract.assert_called_once()
                self.assertEqual({"gradient": 1}, observed[0].outputs)
                with h5py.File(output) as file:
                    self.assertEqual(targeted, "gradient" in file)


if __name__ == "__main__":
    unittest.main()
