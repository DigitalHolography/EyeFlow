"""Tests for two-Gaussian transverse displacement-profile fitting."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pipeline_engine.base import DatasetValue  # noqa: E402
from pipelines.waveform_velocity.displacement_gaussian import (  # noqa: E402
    fit_two_gaussian_profiles,
)
from pipelines.waveform_velocity.profiles import (  # noqa: E402
    pack_cross_section_displacement_profile_outputs,
)


class DisplacementGaussianFitTests(unittest.TestCase):
    def test_pipeline_uses_detected_peaks_for_both_vessels(self) -> None:
        x = np.arange(101, dtype=np.float32)
        curve = (
            0.1
            + 1.1 * np.exp(-np.square(x - 30.0) / (2.0 * 5.0**2))
            + 0.7 * np.exp(-np.square(x - 72.0) / (2.0 * 7.0**2))
        ).astype(np.float32)
        profiles = np.broadcast_to(curve[None, None, None, :], (1, 1, 6, 101))
        displacement = SimpleNamespace(
            transverse_displacement_profiles_unmasked=profiles,
            transverse_displacement_profiles_masked=profiles,
            longitudinal_displacement_profiles_unmasked=profiles,
            longitudinal_displacement_profiles_masked=profiles,
        )
        segments = SimpleNamespace(
            displacements={"level_set_motion": displacement}
        )

        outputs = pack_cross_section_displacement_profile_outputs(
            segments,
            segments,
            np.asarray([0, 2, 5], dtype=np.int32),
        )

        for vessel_name in ("Artery", "Vein"):
            profile_root = (
                "Processing/DisplacementProfiles/level_set_motion/"
                f"{vessel_name}/Transverse"
            )
            metrics_root = (
                "Processing/DisplacementMetrics/level_set_motion/"
                f"{vessel_name}/Transverse"
            )
            source = outputs[
                f"{profile_root}/TransverseDisplacementProfileMaskedMeaned"
            ]
            fitted = outputs[f"{profile_root}/Gaussian_Fit"]
            self.assertEqual(source.data.shape, fitted.data.shape)
            self.assertTrue(outputs[f"{metrics_root}/Gaussian_Fit_Success"].data.all())
            self.assertTrue(
                outputs[
                    f"{metrics_root}/Gaussian_Initialization_Complete"
                ].data.all()
            )
            np.testing.assert_allclose(
                outputs[f"{metrics_root}/Gaussian_Mu_L"].data,
                30.0,
                atol=1e-2,
            )
            np.testing.assert_allclose(
                outputs[f"{metrics_root}/Gaussian_Mu_R"].data,
                72.0,
                atol=1e-2,
            )

    def test_recovers_two_gaussians_and_preserves_nan_mask(self) -> None:
        x = np.arange(121, dtype=np.float32)
        expected = {
            "Gaussian_Baseline": 0.2,
            "Gaussian_A_L": 1.4,
            "Gaussian_A_R": 0.9,
            "Gaussian_Mu_L": 35.0,
            "Gaussian_Mu_R": 82.0,
            "Gaussian_Sigma_L": 6.0,
            "Gaussian_Sigma_R": 9.0,
        }
        curve = (
            expected["Gaussian_Baseline"]
            + expected["Gaussian_A_L"]
            * np.exp(
                -np.square(x - expected["Gaussian_Mu_L"])
                / (2.0 * expected["Gaussian_Sigma_L"] ** 2)
            )
            + expected["Gaussian_A_R"]
            * np.exp(
                -np.square(x - expected["Gaussian_Mu_R"])
                / (2.0 * expected["Gaussian_Sigma_R"] ** 2)
            )
        ).astype(np.float32)
        curve[:4] = np.nan
        curve[-5:] = np.nan
        profiles = np.broadcast_to(curve[:, None, None, None], (121, 2, 1, 1))
        peak_x = np.broadcast_to(
            np.asarray([35.0, 82.0], dtype=np.float32)[:, None, None, None],
            (2, 2, 1, 1),
        )
        peak_y = np.stack(
            (profiles[35], profiles[82]),
            axis=0,
        )

        fitted, metrics = fit_two_gaussian_profiles(
            _profile_value(profiles),
            _peak_value(peak_x, unit="pixels"),
            _peak_value(peak_y, unit="pixels"),
        )

        self.assertEqual(profiles.shape, fitted.data.shape)
        np.testing.assert_array_equal(np.isnan(fitted.data), np.isnan(profiles))
        self.assertTrue(metrics["Gaussian_Fit_Success"].data.all())
        self.assertTrue(metrics["Gaussian_Initialization_Complete"].data.all())
        for name, expected_value in expected.items():
            np.testing.assert_allclose(
                metrics[name].data,
                expected_value,
                rtol=2e-3,
                atol=2e-3,
            )
        np.testing.assert_allclose(
            fitted.data,
            profiles,
            rtol=2e-4,
            atol=2e-4,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            metrics["Gaussian_FWHM_L"].data,
            2.355 * metrics["Gaussian_Sigma_L"].data,
            atol=1e-6,
        )
        np.testing.assert_allclose(
            metrics["Gaussian_Area_R"].data,
            metrics["Gaussian_A_R"].data
            * metrics["Gaussian_Sigma_R"].data
            * np.sqrt(2.0 * np.pi),
            atol=1e-5,
        )
        np.testing.assert_allclose(
            metrics["Gaussian_Peak_Separation"].data,
            metrics["Gaussian_Mu_R"].data - metrics["Gaussian_Mu_L"].data,
            atol=1e-6,
        )
        self.assertLess(float(np.nanmax(metrics["Gaussian_RMSE"].data)), 1e-4)

    def test_missing_peak_marks_fit_as_failed(self) -> None:
        profiles = np.ones((20, 1, 1, 1), dtype=np.float32)
        peak_x = np.asarray([1.0, np.nan], dtype=np.float32).reshape(2, 1, 1, 1)
        peak_y = np.asarray([1.0, np.nan], dtype=np.float32).reshape(2, 1, 1, 1)

        fitted, metrics = fit_two_gaussian_profiles(
            _profile_value(profiles),
            _peak_value(peak_x, unit="pixels"),
            _peak_value(peak_y, unit="pixels"),
        )

        self.assertTrue(np.isnan(fitted.data).all())
        self.assertFalse(metrics["Gaussian_Fit_Success"].data.any())
        self.assertFalse(metrics["Gaussian_Initialization_Complete"].data.any())
        self.assertTrue(np.isnan(metrics["Gaussian_Baseline"].data).all())

    def test_missing_fwhm_intersection_uses_fallback_but_still_fits(self) -> None:
        x = np.arange(50, dtype=np.float32)
        curve = (
            1.0 * np.exp(-np.square(x - 2.0) / (2.0 * 5.0**2))
            + 0.8 * np.exp(-np.square(x - 30.0) / (2.0 * 3.0**2))
        ).astype(np.float32)
        profiles = curve.reshape(50, 1, 1, 1)
        peak_x = np.asarray([2.0, 30.0], dtype=np.float32).reshape(2, 1, 1, 1)
        peak_y = np.asarray(
            [curve[2], curve[30]],
            dtype=np.float32,
        ).reshape(2, 1, 1, 1)

        fitted, metrics = fit_two_gaussian_profiles(
            _profile_value(profiles),
            _peak_value(peak_x, unit="pixels"),
            _peak_value(peak_y, unit="pixels"),
        )

        self.assertTrue(metrics["Gaussian_Fit_Success"].data.all())
        self.assertFalse(metrics["Gaussian_Initialization_Complete"].data.any())
        self.assertTrue(np.isfinite(fitted.data).all())
        np.testing.assert_allclose(
            metrics["Gaussian_Sigma_L"].data,
            5.0,
            atol=2e-2,
        )

    def test_optimizer_failure_clears_fit_and_fitted_parameters(self) -> None:
        x = np.arange(40, dtype=np.float32)
        curve = (
            0.1
            + np.exp(-np.square(x - 10.0) / (2.0 * 3.0**2))
            + 0.8 * np.exp(-np.square(x - 28.0) / (2.0 * 4.0**2))
        ).astype(np.float32)
        profiles = curve.reshape(40, 1, 1, 1)
        peak_x = np.asarray([10.0, 28.0], dtype=np.float32).reshape(2, 1, 1, 1)
        peak_y = np.asarray(
            [curve[10], curve[28]],
            dtype=np.float32,
        ).reshape(2, 1, 1, 1)

        with patch(
            "pipelines.waveform_velocity.displacement_gaussian.least_squares",
            return_value=SimpleNamespace(success=False, x=np.ones(7)),
        ):
            fitted, metrics = fit_two_gaussian_profiles(
                _profile_value(profiles),
                _peak_value(peak_x, unit="pixels"),
                _peak_value(peak_y, unit="pixels"),
            )

        self.assertTrue(np.isnan(fitted.data).all())
        self.assertFalse(metrics["Gaussian_Fit_Success"].data.any())
        for name, metric in metrics.items():
            if name in {
                "Gaussian_Fit_Success",
                "Gaussian_Initialization_Complete",
            }:
                continue
            self.assertTrue(np.isnan(metric.data).all(), name)


def _profile_value(data: np.ndarray) -> DatasetValue:
    return DatasetValue(
        data=data,
        attrs={
            "unit": "pixels",
            "dimDesc": ["x", "beat", "branch", "radius"],
        },
    )


def _peak_value(data: np.ndarray, *, unit: str) -> DatasetValue:
    return DatasetValue(
        data=data,
        attrs={
            "unit": unit,
            "dimDesc": ["peak", "beat", "branch", "radius"],
        },
    )


if __name__ == "__main__":
    unittest.main()
