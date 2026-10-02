"""Tests for method-scoped displacement cross-section HDF5 outputs."""

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

from input_output.writers.h5 import write_value_dataset  # noqa: E402
from input_output.writers.avi import select_display_range  # noqa: E402
from pipelines.displacement_map.constants import (  # noqa: E402
    registration_method_output_name,
)
from pipelines.displacement_map.outputs import OutputCaches  # noqa: E402
from pipelines.waveform_velocity.profiles import (  # noqa: E402
    pack_cross_section_displacement_profile_outputs,
    pack_displacement_magnitude_outputs,
    pack_displacement_profile_outputs,
)
from pipelines.waveform_velocity.segment_maps import (  # noqa: E402
    pack_displacement_segment_map_outputs,
)


class DisplacementOutputTests(unittest.TestCase):
    def assert_no_value_fields(self, h5: h5py.File, *roots: str) -> None:
        value_paths: list[str] = []
        for root in roots:
            names: list[str] = []
            h5[root].visit(names.append)
            value_paths.extend(
                f"{root}/{name}"
                for name in names
                if name.rsplit("/", 1)[-1] == "value"
            )
        self.assertEqual([], value_paths)

    def test_registration_methods_use_canonical_itk_output_names(self) -> None:
        expected = {
            "classic_demons": "Demons",
            "symmetric_forces_demons": "SymmetricForcesDemons",
            "fast_symmetric_demons": "FastSymmetricForcesDemons",
            "diffeomorphic_demons": "DiffeomorphicDemons",
            "level_set_motion": "LevelSetMotion",
            "displacement_field": "DisplacementField",
        }

        self.assertEqual(
            expected,
            {
                method: registration_method_output_name(method)
                for method in expected
            },
        )

    def test_displacement_color_range_uses_temporal_median_image(self) -> None:
        frames = np.asarray(
            [
                [[0.0, 10.0], [1000.0, 50.0]],
                [[2.0, 12.0], [-1000.0, 52.0]],
                [[4.0, 14.0], [6.0, 54.0]],
            ],
            dtype=np.float32,
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            cache = OutputCaches(
                Path(temp_dir),
                frame_count=3,
                height=2,
                width=2,
                save_field=False,
                valid_mask=np.ones((2, 2), dtype=bool),
            )
            try:
                for frame in frames:
                    cache.append(
                        np.zeros((2, 2, 2), dtype=np.float32),
                        frame,
                    )
                self.assertEqual(
                    (2.0, 52.0),
                    select_display_range(
                        cache.magnitude, cache.count, cache.valid_mask,
                        "global-minmax", 1.0, 99.0, 5.0,
                    ),
                )
                expected = tuple(
                    float(value)
                    for value in np.percentile(
                        np.asarray([2.0, 6.0, 12.0, 52.0]),
                        [25.0, 75.0],
                    )
                )
                self.assertEqual(
                    expected,
                    select_display_range(
                        cache.magnitude, cache.count, cache.valid_mask,
                        "percentile", 25.0, 75.0, 5.0,
                    ),
                )
            finally:
                cache.close()

    def test_method_scoped_xy_sum_profiles_and_unmasked_xy_maps_reach_h5(
        self,
    ) -> None:
        artery = _segments()
        vein = _segments()
        boundaries = np.asarray([0, 2, 5], dtype=np.int32)

        outputs = pack_displacement_profile_outputs(
            artery,
            vein,
            boundaries,
        )
        outputs.update(
            pack_displacement_segment_map_outputs(
                artery,
                vein,
                boundaries,
            )
        )

        x_sum_profile_path = (
            "Processing/Displacement/Profiles/FastSymmetricForcesDemons/Artery/"
            "X_sum_displacement_profile"
        )
        y_sum_profile_path = (
            "Processing/Displacement/Profiles/FastSymmetricForcesDemons/Artery/"
            "Y_sum_displacement_profile"
        )
        displacement_map_path = (
            "Processing/Displacement/Map/FastSymmetricForcesDemons/Artery/"
            "PerSegment"
        )
        self.assertEqual(20, len(outputs))
        self.assertFalse(
            any(path.startswith("Processing/Debug/") for path in outputs)
        )
        self.assertFalse(any("Masked" in path for path in outputs))
        self.assertFalse(any("Transverse" in path for path in outputs))
        self.assertFalse(any("Longitudinal" in path for path in outputs))
        self.assertFalse(any("center_" in path for path in outputs))
        self.assertFalse(any("mean_" in path for path in outputs))
        self.assertFalse(
            any("/X_displacement_profile/" in path for path in outputs)
        )
        self.assertFalse(
            any("/Y_displacement_profile/" in path for path in outputs)
        )
        self.assertIn(
            "Processing/Displacement/Profiles/LevelSetMotion/Artery/"
            "X_sum_displacement_profile",
            outputs,
        )
        self.assertIn(
            "Processing/Displacement/Profiles/LevelSetMotion/Vein/"
            "Y_sum_displacement_profile",
            outputs,
        )

        amplitude_profile_path = (
            "Processing/Displacement/Profiles/LevelSetMotion/Vein/"
            "Cross_sectional_radial_movement_amplitude_profile"
        )
        asymmetry_profile_path = (
            "Processing/Displacement/Profiles/LevelSetMotion/Vein/"
            "Cross_sectional_radial_asymmetry_index_profile"
        )
        self.assertIn(amplitude_profile_path, outputs)
        self.assertIn(asymmetry_profile_path, outputs)

        with h5py.File(
            "displacement_outputs.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            self.assert_no_value_fields(
                h5,
                "Processing/Displacement/Profiles",
                "Processing/Displacement/Map",
            )

            x_sum_profile = h5[x_sum_profile_path]
            y_sum_profile = h5[y_sum_profile_path]
            displacement_map = h5[displacement_map_path]
            amplitude_profile = h5[amplitude_profile_path]
            asymmetry_profile = h5[asymmetry_profile_path]

            self.assertIsInstance(
                h5[
                    "Processing/Displacement/Map/"
                    "FastSymmetricForcesDemons/Artery"
                ],
                h5py.Group,
            )
            self.assertIsInstance(
                h5[
                    "Processing/Displacement/Map/"
                    "FastSymmetricForcesDemons/Artery/"
                    "PerSegment"
                ],
                h5py.Dataset,
            )
            self.assertNotIn("Processing/DisplacementMapPerSegment", h5)
            self.assertNotIn("Processing/DisplacementMap", h5)

            self.assertEqual((4, 2, 1, 1), x_sum_profile.shape)
            self.assertEqual((4, 2, 1, 1), y_sum_profile.shape)
            self.assertEqual((4, 3, 4, 2, 1, 1, 2), displacement_map.shape)
            self.assertEqual((4, 2, 1, 1), amplitude_profile.shape)
            self.assertEqual((4, 2, 1, 1), asymmetry_profile.shape)
            self.assertEqual("pixels", x_sum_profile.attrs["unit"])
            self.assertEqual("pixels", y_sum_profile.attrs["unit"])
            self.assertEqual(
                ["time", "beat", "branch", "radius"],
                list(x_sum_profile.attrs["dimDesc"]),
            )
            self.assertEqual(
                ["time", "beat", "branch", "radius"],
                list(y_sum_profile.attrs["dimDesc"]),
            )
            self.assertEqual(
                "rotated_segment_pixel",
                x_sum_profile.attrs["coordinate_system"],
            )
            self.assertEqual(
                "rotated_segment_local",
                x_sum_profile.attrs["component_basis"],
            )
            self.assertEqual("local_x", x_sum_profile.attrs["component"])
            self.assertEqual(
                "full_unmasked_subimage",
                x_sum_profile.attrs["spatial_region"],
            )
            self.assertEqual(
                "sum_over_valid_subimage_pixels",
                x_sum_profile.attrs["spatial_reduction"],
            )
            self.assertEqual("local_y", y_sum_profile.attrs["component"])
            self.assertEqual("pixels", displacement_map.attrs["unit"])
            self.assertEqual(
                ["local_x", "local_y"],
                list(displacement_map.attrs["components"]),
            )
            self.assertEqual("pixels", amplitude_profile.attrs["unit"])
            self.assertEqual("1", asymmetry_profile.attrs["unit"])
            self.assertEqual(
                ["time", "beat", "branch", "radius"],
                list(amplitude_profile.attrs["dimDesc"]),
            )
            self.assertEqual(
                "rotated_vessel_segment_mask",
                amplitude_profile.attrs["centerline_source"],
            )
            self.assertEqual(
                "symmetric_wall_bands_extending_outside_mask",
                asymmetry_profile.attrs["spatial_region"],
            )
            np.testing.assert_allclose(amplitude_profile[...], 6.0, atol=1e-6)
            np.testing.assert_allclose(asymmetry_profile[...], 0.25, atol=1e-6)

            self.assertGreater(
                float(np.nanmin(x_sum_profile[...])),
                0.0,
            )
            self.assertLess(
                float(np.nanmax(y_sum_profile[...])),
                0.0,
            )
            self.assertGreater(
                float(np.nanmin(displacement_map[..., 0])),
                0.0,
            )
            self.assertLess(
                float(np.nanmax(displacement_map[..., 1])),
                0.0,
            )

    def test_displacement_magnitude_is_a_per_segment_per_beat_trace(
        self,
    ) -> None:
        shape = (2, 3, 131)
        displacement = SimpleNamespace(
            x_sum_profile=np.full(
                shape,
                3.0,
                dtype=np.float32,
            ),
            y_sum_profile=np.full(
                shape,
                4.0,
                dtype=np.float32,
            ),
        )
        segments = SimpleNamespace(
            displacements={
                "fast_symmetric_demons": displacement,
                "level_set_motion": displacement,
            },
        )

        outputs = pack_displacement_magnitude_outputs(
            segments,
            segments,
            np.asarray([0, 65, 130], dtype=np.int32),
        )
        magnitude_paths = {
            "Processing/Displacement/Profiles/"
            f"{method}/{vessel}/DisplacementMagnitude"
            for method in ("FastSymmetricForcesDemons", "LevelSetMotion")
            for vessel in ("Artery", "Vein")
        }

        self.assertEqual(magnitude_paths, set(outputs))
        with h5py.File(
            "displacement_magnitude.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            self.assert_no_value_fields(
                h5,
                "Processing/Displacement/Profiles",
            )

            for path in magnitude_paths:
                dataset = h5[path]
                self.assertEqual((128, 2, 3, 2), dataset.shape)
                self.assertEqual(
                    ["time", "beat", "branch", "radius"],
                    list(dataset.attrs["dimDesc"]),
                )
                self.assertEqual("pixels", dataset.attrs["unit"])
                self.assertEqual(
                    "sqrt(x**2 + y**2)",
                    dataset.attrs["magnitude_formula"],
                )
                np.testing.assert_allclose(dataset[...], 5.0, atol=1e-6)

    def test_displacement_axis_profiles_reach_the_requested_h5_paths(
        self,
    ) -> None:
        segments = _segments()
        outputs = pack_cross_section_displacement_profile_outputs(
            segments,
            segments,
            np.asarray([0, 2, 5], dtype=np.int32),
        )
        methods = ("FastSymmetricForcesDemons", "LevelSetMotion")
        vessels = ("Artery", "Vein")
        directions = {"Longitudinal": "y", "Transverse": "x"}
        masks = ("Masked", "Unmasked")
        profile_relative_paths = {
            "tbkr/Profile",
            "bkr/Profile",
            "b/Profile",
            "tbkr/SquaredDeviation",
            "bkr/SquaredDeviation",
        }
        profile_paths = {
            "Processing/Displacement/Profiles/"
            f"{method}/{vessel}/{direction}/{mask}/{relative_path}"
            for method in methods
            for vessel in vessels
            for direction in directions
            for mask in masks
            for relative_path in profile_relative_paths
        }
        artery_root = "Processing/Displacement/Profiles/LevelSetMotion/Artery"
        artery_longitudinal = f"{artery_root}/Longitudinal"
        artery_transverse = f"{artery_root}/Transverse"
        metric_relative_paths = {
            "Position/XMax",
            "Position/YMax",
            "Position/YMaxTemporalVariance",
            "Position/XMaxMeaned",
            "Position/YMaxMeaned",
            "Position/YMaxMeanedTemporalVariance",
            "Area/Left",
            "Area/Right",
            "Area/LeftTemporalVariance",
            "Area/RightTemporalVariance",
        }
        metric_paths = {
            "Processing/Displacement/Metrics/"
            f"{method}/{vessel}/Transverse/{metric_relative_path}"
            for method in methods
            for vessel in vessels
            for metric_relative_path in metric_relative_paths
        }
        expected_paths = profile_paths | metric_paths
        self.assertEqual(expected_paths, set(outputs))
        self.assertFalse(any("Gaussian" in path for path in outputs))

        with h5py.File(
            "displacement_axis_profiles.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            self.assert_no_value_fields(
                h5,
                "Processing/Displacement/Profiles",
                "Processing/Displacement/Metrics",
            )

            metrics_group = h5[
                "Processing/Displacement/Metrics/LevelSetMotion/"
                "Artery/Transverse"
            ]
            self.assertIsInstance(metrics_group["Position"], h5py.Group)
            self.assertIsInstance(metrics_group["Area"], h5py.Group)

            for method in methods:
                for vessel in vessels:
                    for direction, spatial_axis in directions.items():
                        for mask in masks:
                            profile_root = (
                                "Processing/Displacement/Profiles/"
                                f"{method}/{vessel}/{direction}/{mask}"
                            )
                            source = h5[f"{profile_root}/tbkr/Profile"]
                            meaned = h5[f"{profile_root}/bkr/Profile"]
                            global_meaned = h5[
                                f"{profile_root}/b/Profile"
                            ]
                            power = h5[
                                f"{profile_root}/tbkr/SquaredDeviation"
                            ]
                            mean_power = h5[
                                f"{profile_root}/bkr/SquaredDeviation"
                            ]

                            self.assertEqual((181, 4, 2, 1, 1), source.shape)
                            self.assertEqual("pixels", source.attrs["unit"])
                            self.assertEqual(
                                [
                                    spatial_axis,
                                    "time",
                                    "beat",
                                    "branch",
                                    "radius",
                                ],
                                list(source.attrs["dimDesc"]),
                            )
                            self.assertEqual((181, 2, 1, 1), meaned.shape)
                            self.assertEqual(
                                [spatial_axis, "beat", "branch", "radius"],
                                list(meaned.attrs["dimDesc"]),
                            )
                            np.testing.assert_allclose(
                                meaned[...],
                                np.nanmean(source[...], axis=1),
                                atol=1e-6,
                            )
                            self.assertEqual((181, 2), global_meaned.shape)
                            self.assertEqual(
                                [spatial_axis, "beat"],
                                list(global_meaned.attrs["dimDesc"]),
                            )
                            np.testing.assert_allclose(
                                global_meaned[...],
                                np.nanmean(meaned[...], axis=(2, 3)),
                                atol=1e-6,
                            )
                            self.assertEqual(source.shape, power.shape)
                            self.assertEqual("pixels^2", power.attrs["unit"])
                            self.assertEqual(
                                "(D(t) - mean_t(D(t)))**2",
                                power.attrs["formula"],
                            )
                            np.testing.assert_allclose(
                                power[...],
                                (source[...] - meaned[...][:, None, ...]) ** 2,
                                atol=1e-6,
                            )
                            self.assertEqual(meaned.shape, mean_power.shape)
                            np.testing.assert_allclose(
                                mean_power[...],
                                np.nanmean(power[...], axis=1),
                                atol=1e-6,
                            )
            transverse_unmasked = h5[
                f"{artery_transverse}/Unmasked/tbkr/Profile"
            ][...]
            np.testing.assert_allclose(
                h5[f"{artery_transverse}/Masked/tbkr/Profile"][...],
                2.0 * transverse_unmasked,
                atol=1e-5,
            )
            np.testing.assert_allclose(
                h5[f"{artery_longitudinal}/Unmasked/tbkr/Profile"][...],
                3.0 * transverse_unmasked,
                atol=1e-5,
            )
            np.testing.assert_allclose(
                h5[f"{artery_longitudinal}/Masked/tbkr/Profile"][...],
                4.0 * transverse_unmasked,
                atol=1e-5,
            )

    def test_global_displacement_profiles_average_all_valid_segments(
        self,
    ) -> None:
        segments = _segments()
        displacement = segments.displacements["level_set_motion"]
        segment_scales = np.asarray(
            [[1.0, 2.0], [3.0, np.nan]],
            dtype=np.float32,
        )[:, :, None, None]
        for field in (
            "transverse_profiles_unmasked",
            "transverse_profiles_masked",
            "longitudinal_profiles_unmasked",
            "longitudinal_profiles_masked",
        ):
            base_profile = getattr(displacement, field)
            setattr(displacement, field, base_profile * segment_scales)
        segments.displacements = {"level_set_motion": displacement}

        outputs = pack_cross_section_displacement_profile_outputs(
            segments,
            None,
            np.asarray([0, 2, 5], dtype=np.int32),
        )

        for direction in ("Longitudinal", "Transverse"):
            for mask in ("Masked", "Unmasked"):
                root = (
                    "Processing/Displacement/Profiles/LevelSetMotion/Artery/"
                    f"{direction}/{mask}"
                )
                source = outputs[f"{root}/bkr/Profile"]
                global_meaned = outputs[f"{root}/b/Profile"]
                self.assertEqual((181, 2), global_meaned.data.shape)
                np.testing.assert_allclose(
                    global_meaned.data,
                    np.nanmean(source.data, axis=(2, 3)),
                    atol=1e-6,
                )
                np.testing.assert_allclose(
                    global_meaned.data,
                    2.0 * source.data[..., 0, 0],
                    atol=1e-6,
                )

    def test_transverse_peak_metrics_use_first_two_curve_peaks(
        self,
    ) -> None:
        curves = np.asarray(
            [
                [
                    [0, 1, 0, 5, 0, 3, 0, 0, 0],
                    [0, 2, 2, 0, 3, 3, 0, 0, 0],
                ],
                [
                    [0, 0, 1, 4, 1, 0, 0, 0, 0],
                    [0, 1, 2, 3, 4, 5, 6, 7, 8],
                ],
            ],
            dtype=np.float32,
        )
        frame_scales = np.arange(1, 7, dtype=np.float32)[None, None, :, None]
        profiles = curves[:, :, None, :] * frame_scales
        displacement = SimpleNamespace(
            transverse_profiles_unmasked=profiles,
            transverse_profiles_masked=profiles,
            longitudinal_profiles_unmasked=profiles,
            longitudinal_profiles_masked=profiles,
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
            metrics_root = (
                "Processing/Displacement/Metrics/LevelSetMotion/"
                f"{vessel_name}/Transverse"
            )
            position_root = f"{metrics_root}/Position"
            area_root = f"{metrics_root}/Area"
            max_x = outputs[f"{position_root}/XMax"]
            max_y = outputs[f"{position_root}/YMax"]
            diff_y = outputs[f"{position_root}/YMaxTemporalVariance"]
            mean_x = outputs[f"{position_root}/XMaxMeaned"]
            mean_y = outputs[f"{position_root}/YMaxMeaned"]
            mean_diff_y = outputs[
                f"{position_root}/YMaxMeanedTemporalVariance"
            ]
            area_l = outputs[f"{area_root}/Left"]
            area_r = outputs[f"{area_root}/Right"]
            diff_area_l = outputs[f"{area_root}/LeftTemporalVariance"]
            diff_area_r = outputs[f"{area_root}/RightTemporalVariance"]

            self.assertEqual((2, 2, 2, 2), max_x.data.shape)
            self.assertEqual(
                ["peak", "beat", "branch", "radius"],
                max_x.attrs["dimDesc"],
            )
            self.assertEqual(
                ["first_peak_x", "second_peak_x"], max_x.attrs["value_order"]
            )
            np.testing.assert_array_equal(
                max_x.data[:, :, 0, 0],
                np.asarray([[1, 1], [3, 3]], dtype=np.float32),
            )
            np.testing.assert_array_equal(
                max_x.data[:, :, 1, 0],
                np.asarray([[1, 1], [4, 4]], dtype=np.float32),
            )
            np.testing.assert_array_equal(max_x.data[0, :, 0, 1], 3.0)
            self.assertTrue(np.isnan(max_x.data[1, :, 0, 1]).all())
            self.assertTrue(np.isnan(max_x.data[:, :, 1, 1]).all())

            self.assertEqual(max_x.data.shape, max_y.data.shape)
            self.assertEqual(
                ["peak", "beat", "branch", "radius"],
                max_y.attrs["dimDesc"],
            )
            self.assertEqual(max_y.data.shape, diff_y.data.shape)
            self.assertEqual("pixels", max_y.attrs["unit"])
            self.assertEqual("pixels^2", diff_y.attrs["unit"])

            self.assertEqual((2, 2), mean_x.data.shape)
            self.assertEqual(["peak", "beat"], mean_x.attrs["dimDesc"])
            self.assertEqual(mean_x.data.shape, mean_y.data.shape)
            self.assertEqual(mean_x.data.shape, mean_diff_y.data.shape)
            np.testing.assert_allclose(
                mean_x.data,
                np.nanmean(max_x.data, axis=(2, 3)),
                atol=1e-6,
            )
            np.testing.assert_allclose(
                mean_y.data,
                np.nanmean(max_y.data, axis=(2, 3)),
                atol=1e-6,
            )
            np.testing.assert_allclose(
                mean_diff_y.data,
                np.nanmean(diff_y.data, axis=(2, 3)),
                atol=1e-6,
            )
            for area in (area_l, area_r, diff_area_l, diff_area_r):
                self.assertEqual((2, 2, 2), area.data.shape)
                self.assertEqual(
                    ["beat", "branch", "radius"],
                    area.attrs["dimDesc"],
                )
            self.assertEqual("pixels^2", area_l.attrs["unit"])
            self.assertEqual("pixels^2", area_r.attrs["unit"])
            self.assertEqual("pixels^3", diff_area_l.attrs["unit"])
            self.assertEqual("pixels^3", diff_area_r.attrs["unit"])

        artery_metrics_root = (
            "Processing/Displacement/Metrics/LevelSetMotion/"
            "Artery/Transverse"
        )
        transverse_root = (
            "Processing/Displacement/Profiles/LevelSetMotion/"
            "Artery/Transverse"
        )
        position_root = f"{artery_metrics_root}/Position"
        area_root = f"{artery_metrics_root}/Area"
        max_x = outputs[f"{position_root}/XMax"].data
        max_y = outputs[f"{position_root}/YMax"].data
        diff_y = outputs[f"{position_root}/YMaxTemporalVariance"].data
        meaned = outputs[f"{transverse_root}/Masked/bkr/Profile"].data
        mean_power = outputs[
            f"{transverse_root}/Masked/bkr/SquaredDeviation"
        ].data
        area_l = outputs[f"{area_root}/Left"].data
        area_r = outputs[f"{area_root}/Right"].data
        diff_area_l = outputs[f"{area_root}/LeftTemporalVariance"].data
        diff_area_r = outputs[f"{area_root}/RightTemporalVariance"].data
        for peak_index, beat_index, branch_index, radius_index in np.ndindex(
            max_x.shape
        ):
            position = max_x[
                peak_index, beat_index, branch_index, radius_index
            ]
            if not np.isfinite(position):
                continue
            spatial_index = int(position)
            np.testing.assert_allclose(
                max_y[peak_index, beat_index, branch_index, radius_index],
                meaned[spatial_index, beat_index, branch_index, radius_index],
                atol=1e-6,
            )
            np.testing.assert_allclose(
                diff_y[peak_index, beat_index, branch_index, radius_index],
                mean_power[
                    spatial_index, beat_index, branch_index, radius_index
                ],
                atol=1e-6,
            )

        for beat_index, branch_index in np.ndindex((2, 2)):
            curve = meaned[:, beat_index, branch_index, 0]
            expected_area = np.sum((curve[:-1] + curve[1:]) / 2.0)
            np.testing.assert_allclose(
                area_l[beat_index, branch_index, 0]
                + area_r[beat_index, branch_index, 0],
                expected_area,
                atol=1e-6,
            )
            difference_curve = mean_power[:, beat_index, branch_index, 0]
            expected_diff_area = np.sum(
                (difference_curve[:-1] + difference_curve[1:]) / 2.0
            )
            np.testing.assert_allclose(
                diff_area_l[beat_index, branch_index, 0]
                + diff_area_r[beat_index, branch_index, 0],
                expected_diff_area,
                atol=1e-6,
            )

    def test_velocity_only_packing_emits_no_displacement_keys(self) -> None:
        segments = SimpleNamespace(displacements={})
        boundaries = np.asarray([0, 2, 5], dtype=np.int32)
        self.assertEqual(
            {},
            pack_displacement_profile_outputs(segments, segments, boundaries),
        )
        self.assertEqual(
            {},
            pack_displacement_segment_map_outputs(segments, segments, boundaries),
        )

def _segments():
    scalar_maps = np.full((1, 1, 6, 3, 4), -2.0, dtype=np.float32)
    profile_shape = (1, 1, 6, 181)
    profile_time = np.broadcast_to(
        np.arange(1, 7, dtype=np.float32)[None, None, :, None],
        profile_shape,
    )
    x_sum_profile = np.full((1, 1, 6), 12.0, dtype=np.float32)
    y_sum_profile = np.full((1, 1, 6), -24.0, dtype=np.float32)
    radial_amplitude = np.full((1, 1, 6), 3.0, dtype=np.float32)
    radial_asymmetry = np.full((1, 1, 6), 0.25, dtype=np.float32)
    vector_maps = np.stack(
        (
            np.full_like(scalar_maps, 1.0),
            scalar_maps,
        ),
        axis=-1,
    )

    def result(scale: float):
        return SimpleNamespace(
            maps=vector_maps * scale,
            transverse_profiles_unmasked=(
                profile_time * scale
            ),
            transverse_profiles_masked=(
                profile_time * np.float32(2.0) * scale
            ),
            longitudinal_profiles_unmasked=(
                profile_time * np.float32(3.0) * scale
            ),
            longitudinal_profiles_masked=(
                profile_time * np.float32(4.0) * scale
            ),
            x_sum_profile=x_sum_profile * scale,
            y_sum_profile=y_sum_profile * scale,
            radial_movement_amplitude=(
                radial_amplitude * scale
            ),
            radial_asymmetry_index=radial_asymmetry,
        )

    return SimpleNamespace(
        displacements={
            "fast_symmetric_demons": result(1.0),
            "level_set_motion": result(2.0),
        }
    )


if __name__ == "__main__":
    unittest.main()
