from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (  # noqa: E402
    _pack_transverse_profiles,
)
from calculations.blood_flow_velocity.cross_section.profile_processing import (  # noqa: E402
    _matlab_poiseuille_fit,
    interpolate_velocity_profiles_per_beat,
    process_velocity_profiles,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat.signal import (  # noqa: E402
    per_beat_signal_analysis,
)
from calculations.math import rotate_image_with_nan  # noqa: E402
from input_output.output_manager import OutputType  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5 import write_value_dataset  # noqa: E402
from input_output.writers.png import FigureArtifactWriter, write_png_file  # noqa: E402
from pipelines.waveform_velocity_core.figures.profiles import (  # noqa: E402
    _positive_focused_limits,
    export_cross_section_profile_artifacts,
)
from pipelines.waveform_velocity.profiles import (  # noqa: E402
    pack_cross_section_profile_outputs,
)


class TransverseProfilePackingTests(unittest.TestCase):
    def test_variable_width_profiles_are_center_padded_without_interpolation(self) -> None:
        cells = np.empty((1, 2), dtype=object)
        cells[0, 0] = np.asarray([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
        cells[0, 1] = np.asarray(
            [[10, 11, 12, 13, 14], [20, 21, 22, 23, 24]],
            dtype=np.float32,
        )

        profiles, x_pixels, counts = _pack_transverse_profiles(cells, frame_count=2)

        self.assertEqual((1, 2, 2, 5), profiles.shape)
        np.testing.assert_array_equal(x_pixels, [-2, -1, 0, 1, 2])
        np.testing.assert_array_equal(counts, [[3, 5]])
        np.testing.assert_array_equal(profiles[0, 0, :, 1:4], cells[0, 0])
        self.assertTrue(np.all(np.isnan(profiles[0, 0, :, [0, 4]])))
        np.testing.assert_array_equal(profiles[0, 1], cells[0, 1])

    def test_h5_export_contains_raw_and_interpolated_profiles(self) -> None:
        artery = _segments(radius_count=2, branch_count=1)
        vein = _segments(radius_count=2, branch_count=0)
        cycle_boundaries = np.asarray([0, 2, 5], dtype=np.int32)
        metrics = pack_cross_section_profile_outputs(
            artery,
            vein,
            cycle_boundaries,
        )
        schema = EyeFlowOutputPaths.active()
        self.assertEqual(8, len(metrics))

        with h5py.File("profiles.h5", "w", driver="core", backing_store=False) as h5:
            for path, value in metrics.items():
                write_value_dataset(h5, path, value)
            raw_dataset = h5[
                schema.artery_cross_section_profiles.raw_profile.velocity_profile
            ]
            self.assertEqual((7, 4, 2, 1, 2), raw_dataset.shape)
            self.assertEqual(
                list(raw_dataset.attrs["dimDesc"]),
                ["x", "time", "beat", "branch", "radius"],
            )
            raw_coordinate = h5[
                schema.artery_cross_section_profiles
                .raw_profile.transverse_coordinate_micrometers
            ]
            self.assertEqual((7,), raw_coordinate.shape)
            self.assertEqual("um", raw_coordinate.attrs["unit"])
            dataset = h5[
                schema.artery_cross_section_profiles
                .interpolated_profile.velocity_profile
            ]
            self.assertEqual((256, 4, 2, 1, 2), dataset.shape)
            self.assertEqual(
                list(dataset.attrs["dimDesc"]),
                ["x", "time", "beat", "branch", "radius"],
            )
            self.assertEqual("mm/s", dataset.attrs["unit"])
            coordinate = h5[
                schema.artery_cross_section_profiles
                .interpolated_profile.transverse_coordinate_micrometers
            ]
            self.assertEqual((256, 1, 2), coordinate.shape)
            self.assertEqual(
                list(coordinate.attrs["dimDesc"]),
                ["x", "branch", "radius"],
            )
            self.assertEqual("um", coordinate.attrs["unit"])
            empty_raw = h5[
                schema.vein_cross_section_profiles.raw_profile.velocity_profile
            ]
            self.assertEqual((7, 4, 2, 0, 2), empty_raw.shape)
            empty = h5[
                schema.vein_cross_section_profiles
                .interpolated_profile.velocity_profile
            ]
            self.assertEqual((256, 4, 2, 0, 2), empty.shape)

    def test_profile_time_axis_matches_standard_per_beat_interpolation(self) -> None:
        segments = _segments(radius_count=1, branch_count=1)
        cycle_boundaries = np.asarray([0, 2, 5], dtype=np.int32)

        result = interpolate_velocity_profiles_per_beat(
            segments.centered_velocity_profiles,
            cycle_boundaries,
        )
        expected = per_beat_signal_analysis(
            segments.centered_velocity_profiles[0, 0, :, 100],
            cycle_boundaries,
            1,
            index_base=0,
        ).velocity_signal_per_beat.T

        self.assertEqual((256, 4, 2, 1, 1), result.shape)
        np.testing.assert_allclose(
            result[100, :, :, 0, 0],
            expected,
            equal_nan=True,
        )

    def test_centering_anchors_edges_without_clipping_negative_velocity(self) -> None:
        raw = np.asarray(
            [[[
                [-8.0, 2.0, 8.0, -2.0, 7.0, 1.0, -9.0],
                [-8.0, 2.0, 8.0, 22.0, 7.0, 1.0, -9.0],
            ]]],
            dtype=np.float32,
        )

        result = process_velocity_profiles(
            raw,
            pixel_size_mm=0.01,
            velocity_profile_threshold=0.5,
            interpolation_points=33,
        )

        centered = result.centered_velocity[0, 0, 0]
        np.testing.assert_array_equal(result.raw_profile.velocity, raw)
        np.testing.assert_array_equal(
            result.raw_profile.x_micrometers,
            result.raw_x_micrometers,
        )
        np.testing.assert_array_equal(
            result.interpolated_profile.velocity,
            result.centered_velocity,
        )
        np.testing.assert_array_equal(
            result.interpolated_profile.x_micrometers,
            result.centered_x_micrometers,
        )
        self.assertAlmostEqual(0.0, float(centered[0]), places=5)
        self.assertAlmostEqual(0.0, float(centered[-1]), places=5)
        self.assertTrue(np.allclose(
            result.centered_x_micrometers[0, 0],
            -result.centered_x_micrometers[0, 0, ::-1],
        ))
        self.assertGreater(np.count_nonzero(centered < 0), 0)

    def test_poiseuille_fit_matches_matlab_custom_poly_functions(self) -> None:
        x_um = np.arange(-30, 31, 10, dtype=np.float32)
        profile = np.asarray([-8, 2, 8, 10, 7, 1, -9], dtype=np.float32)

        result = _matlab_poiseuille_fit(x_um, profile, 0.5)

        self.assertIsNotNone(result)
        coefficients, origin_um, roots_um, r_squared = result
        np.testing.assert_allclose(
            coefficients,
            [-0.025, -0.05, 10.0],
            rtol=1e-6,
        )
        self.assertEqual(0.0, origin_um)
        np.testing.assert_allclose(
            roots_um,
            [-15.1774468788, 13.1774468788],
            rtol=1e-6,
        )
        self.assertAlmostEqual(1.0, r_squared)

    def test_nan_rotation_matches_matlab_imrotatecustom(self) -> None:
        image = np.full((5, 5), np.nan, dtype=np.float32)
        image[1:4, 1:4] = np.arange(1, 10, dtype=np.float32).reshape(
            3,
            3,
            order="F",
        )
        expected = np.asarray(
            [
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [np.nan, 1.839746, 5.633975, 4.839746, np.nan],
                [np.nan, 1.901924, 5.0, 8.098076, np.nan],
                [np.nan, 1.5, 4.366025, 4.5, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan],
            ],
            dtype=np.float32,
        )

        np.testing.assert_allclose(
            rotate_image_with_nan(image, 30.0),
            expected,
            rtol=1e-6,
            equal_nan=True,
        )


class ProfileArtifactTests(unittest.TestCase):
    def test_profile_plot_limits_focus_positive_flow(self) -> None:
        y_min, y_max = _positive_focused_limits(
            np.asarray([-40.0, -12.0, 0.0, 4.0, 3.0], dtype=np.float32)
        )

        self.assertAlmostEqual(-1.12, y_min, places=5)
        self.assertAlmostEqual(4.48, y_max, places=5)

    def test_pngs_and_gifs_are_written_for_available_profiles(self) -> None:
        try:
            import matplotlib  # noqa: F401
            import PIL  # noqa: F401
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

            paths = export_cross_section_profile_artifacts(
                writer,
                context,
                max_gif_frames=4,
            )

            self.assertEqual(8, len(paths))
            for path in paths:
                self.assertTrue(path.is_file(), path)
                self.assertGreater(path.stat().st_size, 0, path)
            self.assertEqual(2, sum(path.suffix == ".gif" for path in paths))
            self.assertTrue(
                all("velocityProfiles" in path.parts for path in paths)
            )
            self.assertEqual(
                2,
                sum("poiseuille_profile" in path.name for path in paths),
            )


def _segments(*, radius_count: int, branch_count: int):
    frames = 6
    width = 7
    x_pixels = np.arange(-3, 4, dtype=np.float32)
    profiles = np.empty(
        (radius_count, branch_count, frames, width),
        dtype=np.float32,
    )
    for radius in range(radius_count):
        for branch in range(branch_count):
            for frame in range(frames):
                center = 10.0 + frame + radius + branch
                profiles[radius, branch, frame] = (
                    center - 0.8 * x_pixels**2 + 0.3 * x_pixels
                )
    processed = process_velocity_profiles(
        profiles,
        pixel_size_mm=0.01,
        velocity_profile_threshold=0.5,
    )
    return SimpleNamespace(
        branch_ids=np.arange(1, branch_count + 1, dtype=np.int32),
        velocity_profiles=profiles,
        profile_x_micrometers=processed.raw_x_micrometers,
        profile_sample_count=np.full(
            (radius_count, branch_count),
            width,
            dtype=np.int32,
        ),
        profile_rotation_degrees=np.zeros(
            (radius_count, branch_count),
            dtype=np.float32,
        ),
        centered_velocity_profiles=processed.centered_velocity,
        centered_profile_x_micrometers=processed.centered_x_micrometers,
        profile_center_micrometers=processed.center_micrometers,
        profile_lumen_edges_micrometers=processed.lumen_edges_micrometers,
        profile_centering_fit_r_squared=processed.centering_fit_r_squared,
        poiseuille_coefficients=processed.poiseuille_coefficients,
        poiseuille_origin_micrometers=(
            processed.poiseuille_origin_micrometers
        ),
        poiseuille_roots_micrometers=(
            processed.poiseuille_roots_micrometers
        ),
        poiseuille_r_squared=processed.poiseuille_r_squared,
        poiseuille_profile_spatial_std=np.zeros(
            (radius_count, branch_count, width),
            dtype=np.float32,
        ),
    )


class _FakeOutput:
    available = True

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manager = SimpleNamespace(layout=SimpleNamespace(stem="sample"))

    def path_for(self, output_type: OutputType, filename: str | None = None) -> Path:
        if output_type not in (OutputType.PNG, OutputType.GIF):
            raise AssertionError(output_type)
        return self.root / output_type.value / (filename or "sample")

    def write_png(self, output, filename: str | None = None) -> Path:
        return write_png_file(self.path_for(OutputType.PNG, filename), output)


if __name__ == "__main__":
    unittest.main()
