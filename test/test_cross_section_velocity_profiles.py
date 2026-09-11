"""Standard cross-section profile packing and interpolation tests."""

from __future__ import annotations

import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity.cross_section.profile_processing import (  # noqa: E402
    interpolate_velocity_profiles_per_beat,
)
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (  # noqa: E402
    CrossSectionProfileOutputs,
    _gpu_nanmean,
)
from calculations.compute_backend import optional_cupy_backend  # noqa: E402
from calculations.topology import profiles as topology_profiles  # noqa: E402
from calculations.math import rotate_image_with_nan  # noqa: E402
from input_output.output_manager import OutputType  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5 import write_value_dataset  # noqa: E402
from input_output.writers.png import FigureArtifactWriter, write_png_file  # noqa: E402
from pipelines.waveform_velocity.profiles import pack_cross_section_profile_outputs  # noqa: E402
from pipelines.waveform_velocity_core.figures.profiles import (  # noqa: E402
    _finite_median,
    _hierarchical_profile_median,
    _nanmedian,
    _positive_focused_limits,
    export_cross_section_profile_artifacts,
)


class CrossSectionProfilePackingTests(unittest.TestCase):
    def test_gpu_nanmean_avoids_unsupported_ufunc_where_argument(self) -> None:
        backend = optional_cupy_backend()
        if backend is None:
            self.skipTest("CuPy/CUDA is unavailable.")

        cupy = backend.cupy
        values = cupy.asarray(
            [
                [[1.0, cupy.nan], [3.0, cupy.nan]],
                [[2.0, 4.0], [4.0, 6.0]],
            ],
            dtype=cupy.float32,
        )

        result = cupy.asnumpy(_gpu_nanmean(values, axis=1, cupy=cupy))

        np.testing.assert_allclose(
            result,
            np.asarray([[2.0, np.nan], [3.0, 5.0]], dtype=np.float32),
            equal_nan=True,
        )

    def test_obsolete_centering_poiseuille_and_inverse_fit_interfaces_are_gone(self) -> None:
        removed = {
            "centered_velocity_profiles",
            "centered_profile_x_micrometers",
            "profile_center_micrometers",
            "profile_lumen_edges_micrometers",
            "profile_centering_fit_r_squared",
            "poiseuille_coefficients",
            "poiseuille_origin_micrometers",
            "poiseuille_roots_micrometers",
            "poiseuille_r_squared",
        }
        self.assertTrue(removed.isdisjoint(CrossSectionProfileOutputs.__dataclass_fields__))
        self.assertFalse(
            hasattr(topology_profiles, "fit_inverse_parabola_profiles_with_roots")
        )
        self.assertFalse((SRC_DIR / "pipelines/waveform_velocity/flow_asymmetry.py").exists())

    def test_h5_export_contains_only_four_standard_profiles_per_vessel(self) -> None:
        artery = _segments(radius_count=2, branch_count=1)
        vein = _segments(radius_count=2, branch_count=0)
        metrics = pack_cross_section_profile_outputs(artery, vein, [0, 2, 5])
        self.assertEqual(8, len(metrics))

        schema = EyeFlowOutputPaths.active()
        expected = set()
        for paths in (schema.artery_velocity_profiles, schema.vein_velocity_profiles):
            expected.update(
                {
                    paths.transverse_velocity_profile_unmasked,
                    paths.transverse_velocity_profile_masked,
                    paths.longitudinal_velocity_profile_unmasked,
                    paths.longitudinal_velocity_profile_masked,
                }
            )
        self.assertEqual(expected, set(metrics))
        self.assertFalse(any("Fit" in path or "FFT" in path for path in metrics))

        with h5py.File("profiles.h5", "w", driver="core", backing_store=False) as h5:
            for path, value in metrics.items():
                write_value_dataset(h5, path, value)
            transverse = h5[
                schema.artery_velocity_profiles.transverse_velocity_profile_masked
            ]
            longitudinal = h5[
                schema.artery_velocity_profiles.longitudinal_velocity_profile_masked
            ]
            self.assertEqual((181, 4, 2, 1, 2), transverse.shape)
            self.assertEqual(transverse.shape, longitudinal.shape)
            self.assertEqual(
                ["x", "time", "beat", "branch", "radius"],
                list(transverse.attrs["dimDesc"]),
            )
            self.assertEqual(
                ["y", "time", "beat", "branch", "radius"],
                list(longitudinal.attrs["dimDesc"]),
            )
            self.assertNotIn("Processing/CrossSections", h5)

    def test_interpolation_vectorizes_spatial_samples_and_skips_invalid_slots(self) -> None:
        profiles = np.empty((1, 2, 6, 3), dtype=np.float32)
        profiles[0, 0] = np.arange(18, dtype=np.float32).reshape(6, 3)
        profiles[0, 1] = 99
        from scipy.signal import resample as scipy_resample

        with patch(
            "calculations.blood_flow_velocity.cross_section.profile_processing.resample",
            wraps=scipy_resample,
        ) as resample:
            result = interpolate_velocity_profiles_per_beat(
                profiles,
                [0, 2, 5],
                valid_segments=np.array([[True, False]]),
            )

        self.assertEqual((3, 4, 2, 2, 1), result.shape)
        self.assertEqual(2, resample.call_count)
        self.assertTrue(np.all(np.isnan(result[..., 1, 0])))
        self.assertTrue(np.any(np.isfinite(result[..., 0, 0])))

    def test_nan_rotation_interpolates_only_finite_values(self) -> None:
        image = np.full((5, 5), np.nan, dtype=np.float32)
        image[1:4, 1:4] = np.arange(1, 10, dtype=np.float32).reshape(3, 3, order="F")
        rotated = rotate_image_with_nan(image, 30.0)
        self.assertTrue(np.all(np.isnan(rotated[[0, -1]])))
        self.assertAlmostEqual(5.0, float(rotated[2, 2]), places=5)


class ProfileArtifactTests(unittest.TestCase):
    def test_profile_median_handles_sparse_branches_by_radius(self) -> None:
        values = np.arange(2 * 12 * 3 * 4, dtype=np.float32).reshape(2, 12, 3, 4)
        values[0, :10, 1, 2] = np.nan
        values[1, :, 2, 3] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            expected = np.nanmedian(np.nanmedian(values, axis=1), axis=0)
        np.testing.assert_allclose(
            _hierarchical_profile_median(values), expected, equal_nan=True
        )

    def test_profile_median_has_sparse_array_fallback(self) -> None:
        values = np.arange(2 * 12 * 3, dtype=np.float32).reshape(2, 12, 3)
        expected = _finite_median(values, axis=1)
        with patch(
            "pipelines.waveform_velocity_core.figures.profiles.np.nanmedian",
            side_effect=IndexError("sparse partition failure"),
        ):
            actual = _nanmedian(values, axis=1)
        np.testing.assert_allclose(actual, expected, equal_nan=True)

    def test_profile_plot_limits_focus_positive_flow(self) -> None:
        y_min, y_max = _positive_focused_limits(
            np.array([-40.0, -12.0, 0.0, 4.0, 3.0])
        )
        self.assertAlmostEqual(-1.12, y_min)
        self.assertAlmostEqual(4.48, y_max)

    def test_only_raw_profile_diagnostics_are_written(self) -> None:
        try:
            import matplotlib  # noqa: F401
        except ModuleNotFoundError:
            self.skipTest("profile artifact dependencies are not installed")
        with tempfile.TemporaryDirectory() as temp_dir:
            output = _FakeOutput(Path(temp_dir))
            writer = FigureArtifactWriter(output, "sample")
            context = SimpleNamespace(
                source_data=SimpleNamespace(timing=SimpleNamespace(dt_seconds=0.05)),
                artery_segment_result=_segments(radius_count=1, branch_count=1),
                vein_segment_result=_segments(radius_count=1, branch_count=1),
            )
            paths = export_cross_section_profile_artifacts(writer, context)
            self.assertEqual(2, len(paths))
            self.assertTrue(all(path.is_file() for path in paths))
            self.assertTrue(all("profile_map" in path.name for path in paths))


def _segments(*, radius_count: int, branch_count: int):
    frames = 6
    width = 181
    x = np.linspace(-3.0, 3.0, width, dtype=np.float32)
    profiles = np.empty((radius_count, branch_count, frames, width), dtype=np.float32)
    for radius, branch, frame in np.ndindex(radius_count, branch_count, frames):
        profiles[radius, branch, frame] = 10 + frame - 0.8 * x**2 + 0.3 * x
    masked = profiles.copy()
    if branch_count:
        masked[..., 0] = np.nan
    return SimpleNamespace(
        topology=SimpleNamespace(
            valid_segments=np.ones((radius_count, branch_count), dtype=bool)
        ),
        branch_ids=np.arange(1, branch_count + 1, dtype=np.int32),
        velocity_profiles=profiles,
        transverse_velocity_profiles_masked=masked,
        longitudinal_velocity_profiles_unmasked=profiles + np.float32(50),
        longitudinal_velocity_profiles_masked=masked + np.float32(100),
    )


class _FakeOutput:
    available = True

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manager = SimpleNamespace(layout=SimpleNamespace(stem="sample"))

    def path_for(self, output_type: OutputType, filename: str | None = None) -> Path:
        return self.root / output_type.value / (filename or "sample")

    def write_png(self, output, filename: str | None = None) -> Path:
        return write_png_file(self.path_for(OutputType.PNG, filename), output)


if __name__ == "__main__":
    unittest.main()
