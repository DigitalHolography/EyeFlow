"""Tests for spatial-gradient calculations and the retired full-frame export."""

import unittest

import numpy as np

import pipelines
from calculations.math.spatial_gradient import spatial_gradient
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG


class SpatialGradientPipelineTests(unittest.TestCase):
    def test_spatial_gradient_uses_sobel_magnitude(self) -> None:
        ramp = np.tile(np.arange(5, dtype=np.float32), (4, 1))

        gradient = spatial_gradient(ramp)

        np.testing.assert_allclose(gradient[:, 1:-1], 8.0)
        np.testing.assert_allclose(gradient[:, (0, -1)], 4.0)

    def test_full_frame_export_is_absent_from_catalog_and_downstream_plans(self) -> None:
        pipelines.load_pipeline_catalog()

        self.assertNotIn("spatial_gradient_moment0", PIPELINE_REGISTRY)
        dag = PipelineDAG(PIPELINE_REGISTRY.values())
        for target in (
            "waveform_velocity",
            "waveform_shape_metrics",
            "absolute_waveform_metrics",
            "pdf_report",
        ):
            with self.subTest(target=target):
                plan = dag.resolve_targets([target])
                self.assertNotIn("spatial_gradient_moment0", plan.names)
                self.assertLess(
                    plan.names.index("waveform_velocity_core"),
                    plan.names.index("waveform_velocity"),
                )


if __name__ == "__main__":
    unittest.main()
