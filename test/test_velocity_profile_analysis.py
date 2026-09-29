"""Weighted artery and vein profile-analysis tests."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

from pipeline_engine.context import PipelineH5Output
from pipeline_engine.base import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.velocity_profile_analysis import fitting
from pipelines.velocity_profile_analysis.runner import (
    OUTPUT_ROOT,
    SOURCE_PATHS,
    run_velocity_profile_analysis,
)


def _analyze(y, **kwargs):
    values = np.asarray(y)[:, None, None, None, None]
    return fitting.analyze_velocity_profiles(values, **kwargs)


def _scalar(result, name):
    return result[name].item()


class VelocityProfileFittingTests(unittest.TestCase):
    def test_border_weights_are_quadratic_over_complete_domain(self) -> None:
        np.testing.assert_array_equal(
            fitting.border_weights(9),
            [0, 0.4375, 0.75, 0.9375, 1, 0.9375, 0.75, 0.4375, 0],
        )
        np.testing.assert_allclose(
            fitting.border_weights(6),
            [0, 0.64, 0.96, 0.96, 0.64, 0],
        )
        np.testing.assert_array_equal(fitting.border_weights(1), [1])
        self.assertEqual(0, fitting.border_weights(0).size)

    def test_border_weight_power_is_parameterized(self) -> None:
        np.testing.assert_array_equal(
            fitting.border_weights(5, power=1),
            [0, 0.5, 1, 0.5, 0],
        )
        np.testing.assert_array_equal(
            fitting.border_weights(5, power=4),
            [0, 0.9375, 1, 0.9375, 0],
        )

    def test_border_weight_power_must_be_finite_positive_real(self) -> None:
        for power in (0, -1, np.nan, np.inf, -np.inf, True, "2", 1 + 0j):
            with self.subTest(power=power):
                with self.assertRaisesRegex(ValueError, "finite positive real"):
                    fitting.border_weights(9, power=power)
                with self.assertRaisesRegex(ValueError, "finite positive real"):
                    _analyze(np.arange(9.0), weight_power=power)

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

    def test_axis_order_time_blocking_and_shared_solves(self) -> None:
        shape = (2, 2, 3, 2)
        factors = np.arange(1, np.prod(shape) + 1).reshape(shape)
        x = np.arange(9.0)
        values = (-(x - 1.25) * (x - 6.75))[:, None, None, None, None] * factors[None]
        result = fitting.analyze_velocity_profiles(values, time_block_size=1)
        np.testing.assert_allclose(result["a"], -factors)
        self.assertTrue(all(value.shape == shape for value in result.values()))

        shared = np.repeat((-(x - 1) * (x - 7))[:, None, None, None, None], 4, axis=1)
        shared[2, :2] = np.nan
        with patch.object(np.linalg, "lstsq", wraps=np.linalg.lstsq) as solve:
            fitting.analyze_velocity_profiles(shared)
        self.assertEqual(2, solve.call_count)

    def test_bounded_dataset_reads(self) -> None:
        class SlabOnly:
            shape = (9, 7, 1, 1, 1)
            dtype = np.dtype("float32")

            def __getitem__(self, key):
                self_test.assertLessEqual(key[1].stop - key[1].start, 2)
                return np.ones((9, key[1].stop - key[1].start))

            def __array__(self, *args, **kwargs):
                raise AssertionError("full dataset read")

        self_test = self
        result = fitting.analyze_velocity_profiles(SlabOnly(), time_block_size=2)
        self.assertEqual((7, 1, 1, 1), result["a"].shape)


class VelocityProfileAnalysisPipelineTests(unittest.TestCase):
    def test_pipeline_is_visible_and_depends_on_waveform_velocity(self) -> None:
        load_pipeline_catalog()
        descriptor = PIPELINE_REGISTRY["velocity_profile_analysis"]
        self.assertEqual(("waveform_velocity",), descriptor.dag_requires)
        self.assertEqual("visible", descriptor.visibility)

    def test_both_vessels_publish_the_complete_identical_output_set(self) -> None:
        values = np.broadcast_to(
            (10 - (np.arange(9.0) - 4) ** 2)[:, None, None, None, None],
            (9, 3, 1, 2, 1),
        )
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "analysis.h5"
            with h5py.File(filename, "w") as h5:
                for source_path in SOURCE_PATHS.values():
                    h5.create_dataset(source_path, data=values)
                output = PipelineH5Output(h5)
                results = run_velocity_profile_analysis(
                    SimpleNamespace(output=SimpleNamespace(h5=output))
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

    def test_missing_source_names_the_missing_vessel_and_path(self) -> None:
        for missing in SOURCE_PATHS:
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as directory:
                with h5py.File(Path(directory) / "missing.h5", "w") as h5:
                    for vessel, path in SOURCE_PATHS.items():
                        if vessel != missing:
                            h5.create_dataset(path, data=np.ones((9, 1, 1, 1, 1)))
                    ctx = SimpleNamespace(
                        output=SimpleNamespace(h5=PipelineH5Output(h5))
                    )
                    with self.assertRaisesRegex(KeyError, f"{missing.lower()}.*{missing}"):
                        run_velocity_profile_analysis(ctx)


if __name__ == "__main__":
    unittest.main()
