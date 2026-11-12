"""Weighted artery and vein profile-analysis tests."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from calculations.vessel_segments.profiles.fits import quadratic as fitting
from pipeline_engine.base import DatasetValue
from pipeline_engine.context import PipelineH5Output
from pipelines.velocity_analysis.outputs.profile_analysis import (
    OUTPUT_ROOT,
    SOURCE_PATHS,
    pack_velocity_profile_analysis_outputs,
)


def _analyze(y, **kwargs):
    values = np.asarray(y)[:, None, None, None, None]
    return fitting.fit_quadratic_profiles(values, **kwargs)


def _scalar(result, name):
    return result[name].item()


class VelocityProfileFittingTests(unittest.TestCase):
    def test_coefficients_fractional_roots_center_and_areas(self) -> None:
        x = np.arange(11.0)
        y = -2 * (x - 1.25) * (x - 8.75)
        result = _analyze(y)
        np.testing.assert_allclose(
            [_scalar(result, name) for name in ("a", "b", "c")],
            [-2, 20, -21.875],
            rtol=1e-6,
        )
        np.testing.assert_allclose(
            [
                _scalar(result, name)
                for name in ("index_center", "index_left_zero", "index_right_zero")
            ],
            [5, 1.25, 8.75],
            rtol=1e-6,
        )
        support = (x >= 1.25) & (x <= 8.75)
        self.assertAlmostEqual(float(y[support].sum()), _scalar(result, "Qv"), places=4)
        self.assertAlmostEqual(_scalar(result, "Qv"), _scalar(result, "Qv_fit"), places=4)
        self.assertEqual(11, _scalar(result, "n_fit_samples"))
        self.assertEqual(7, _scalar(result, "n_area_samples"))

    def test_weighted_fit_handles_missing_and_nonfinite_samples(self) -> None:
        x = np.arange(17.0)
        y = 10 - 0.3 * (x - 7) ** 2 + np.sin(x)
        y[[0, 7, 15]] = [np.nan, np.inf, np.nan]
        finite = np.isfinite(y)
        weights = fitting.border_weights(len(x))[finite]
        matrix = np.column_stack((x[finite] ** 2, x[finite], np.ones(finite.sum())))
        expected = np.linalg.lstsq(
            matrix * np.sqrt(weights)[:, None],
            y[finite] * np.sqrt(weights),
            rcond=None,
        )[0]
        result = _analyze(y)
        np.testing.assert_allclose(
            [_scalar(result, name) for name in ("a", "b", "c")],
            expected,
            rtol=1e-6,
        )
        self.assertEqual(int(finite.sum()), _scalar(result, "n_fit_samples"))
        for name in ("fit_rss", "fit_rmse", "fit_weighted_rss", "fit_weighted_rmse"):
            self.assertTrue(np.isfinite(_scalar(result, name)))

    def test_unsuitable_curvature_and_insufficient_samples(self) -> None:
        for y in (
            (np.arange(9.0) - 4) ** 2,
            np.full(9, 3.0),
            np.arange(9.0) + 1,
        ):
            result = _analyze(y)
            self.assertTrue(np.isfinite(_scalar(result, "a")))
            self.assertTrue(np.isnan(_scalar(result, "index_left_zero")))
            self.assertTrue(np.isnan(_scalar(result, "Qv")))
        result = _analyze([np.nan, 1, np.inf, 2, np.nan])
        self.assertEqual(2, _scalar(result, "n_fit_samples"))
        self.assertTrue(all(np.isnan(_scalar(result, name)) for name in fitting.FLOAT_OUTPUTS))

class VelocityProfileAnalysisOptionTests(unittest.TestCase):
    def test_both_vessels_publish_the_complete_identical_output_set(self) -> None:
        values = np.broadcast_to(
            (10 - (np.arange(9.0) - 4) ** 2)[:, None, None, None, None],
            (9, 3, 1, 2, 1),
        )
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "analysis.h5"
            with h5py.File(filename, "w") as h5:
                output = PipelineH5Output(h5)
                results = pack_velocity_profile_analysis_outputs(
                    {
                        source_path.lstrip("/"): DatasetValue(values)
                        for source_path in SOURCE_PATHS.values()
                    }
                )
                output.write_many(results)
                expected = set(fitting.FLOAT_OUTPUTS + fitting.COUNT_OUTPUTS)
                for vessel, source_path in SOURCE_PATHS.items():
                    group = h5[f"{OUTPUT_ROOT}/{vessel}"]
                    self.assertEqual(expected, set(group))
                    for name in expected:
                        dataset = group[name]["value"]
                        self.assertEqual((3, 1, 2, 1), dataset.shape)
                        dtype = np.int32 if name in fitting.COUNT_OUTPUTS else np.float32
                        self.assertEqual(np.dtype(dtype), dataset.dtype)
                        self.assertEqual(
                            ["time", "beat", "branch", "radius"],
                            list(dataset.attrs["dimDesc"]),
                        )
                        self.assertEqual(source_path, dataset.attrs["source_path"])
                        self.assertEqual(
                            fitting.DEFAULT_WEIGHT_POWER,
                            dataset.attrs["weight_power"],
                        )
                        self.assertEqual(
                            "u=x/(Nx-1); d=abs(2*u-1); w=1-d^p",
                            dataset.attrs["weight_definition"],
                        )


if __name__ == "__main__":
    unittest.main()
