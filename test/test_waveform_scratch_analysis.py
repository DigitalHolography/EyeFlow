from __future__ import annotations

import unittest

import numpy as np

from input_output.schema import ImageMaps
from pipelines.velocity.estimation import (
    _bounded_inpaint_result,
    _signed_rms_difference,
    estimate_retinal_velocity,
)


class ChunkedAnalysisAndSchemaTests(unittest.TestCase):
    def test_inpaint_result_is_finite_and_bounded_by_each_frame_background(self) -> None:
        source = np.asarray(
            [
                [[1.0, 2.0], [3.0, 4.0]],
                [[10.0, 20.0], [30.0, 40.0]],
            ],
            dtype=np.float32,
        )
        mask = np.asarray([[True, False], [False, False]])
        unstable = np.asarray(
            [
                [[np.inf, -100.0], [3.0, 100.0]],
                [[np.nan, -100.0], [30.0, 100.0]],
            ],
            dtype=np.float32,
        )

        actual = _bounded_inpaint_result(unstable, source, mask)

        self.assertTrue(np.all(np.isfinite(actual)))
        self.assertGreaterEqual(float(actual[0].min()), 2.0)
        self.assertLessEqual(float(actual[0].max()), 4.0)
        self.assertGreaterEqual(float(actual[1].min()), 20.0)
        self.assertLessEqual(float(actual[1].max()), 40.0)

    def test_signed_rms_difference_avoids_float32_square_overflow(self) -> None:
        foreground = np.asarray([3.0e20], dtype=np.float32)
        background = np.asarray([1.0e20], dtype=np.float32)

        with np.errstate(over="raise"):
            actual = _signed_rms_difference(foreground, background)

        self.assertTrue(np.all(np.isfinite(actual)))
        np.testing.assert_allclose(
            actual,
            np.sqrt(8.0) * np.float32(1.0e20),
            rtol=1e-6,
        )


    def test_disabled_vein_keeps_artery_background_and_returns_nan_signals(
        self,
    ) -> None:
        rng = np.random.default_rng(19)
        moment0 = (1.0 + rng.random((3, 16, 16))).astype(np.float32)
        moment2 = (2.0 + rng.random((3, 16, 16))).astype(np.float32)
        artery = np.zeros((16, 16), dtype=bool)
        source_vein = np.zeros_like(artery)
        artery[6, 6] = True
        source_vein[9, 9] = True
        background = artery | source_vein

        expected = estimate_retinal_velocity(
            image_maps=ImageMaps(moment0=moment0, moment2=moment2),
            artery_mask=artery,
            vein_mask=source_vein,
            optic_disc_center=(8.0, 8.0),
            local_background_dist=1,
            retain_velocity_video=False,
        )
        actual = estimate_retinal_velocity(
            image_maps=ImageMaps(moment0=moment0, moment2=moment2),
            artery_mask=artery,
            vein_mask=np.zeros_like(source_vein),
            background_mask=background,
            optic_disc_center=(8.0, 8.0),
            local_background_dist=1,
            retain_velocity_video=False,
        )

        np.testing.assert_array_equal(
            actual.artery.velocity,
            expected.artery.velocity,
        )
        for values in (
            actual.vein.velocity,
            actual.vein.frms,
            actual.vein.frms_background,
            actual.vein.delta_frms,
        ):
            self.assertTrue(np.all(np.isnan(values)))


if __name__ == "__main__":
    unittest.main()
