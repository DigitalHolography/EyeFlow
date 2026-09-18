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
from pipelines.displacement_map.outputs import (  # noqa: E402
    OutputCaches,
    select_display_range,
)
from pipelines.waveform_velocity.profiles import (  # noqa: E402
    _combined_displacement_magnitude_dataset,
    pack_cross_section_displacement_profile_outputs,
    pack_displacement_magnitude_outputs,
    pack_displacement_profile_outputs,
)
from pipelines.waveform_velocity.segment_maps import (  # noqa: E402
    pack_displacement_segment_map_outputs,
)
from pipelines.waveform_velocity_core.runner import (  # noqa: E402
    _load_displacement_maps,
)
from pipelines.displacement_map.runner import (  # noqa: E402
    DISPLACEMENT_MAP_STATE,
    DisplacementMapArtifacts,
)


class DisplacementOutputTests(unittest.TestCase):
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
                    select_display_range(cache, "global-minmax", 1.0, 99.0, 5.0),
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
                    select_display_range(cache, "percentile", 25.0, 75.0, 5.0),
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
            "Processing/DisplacementProfiles/fast_symmetric_demons/Artery/"
            "X_sum_displacement_profile/value"
        )
        y_sum_profile_path = (
            "Processing/DisplacementProfiles/fast_symmetric_demons/Artery/"
            "Y_sum_displacement_profile/value"
        )
        displacement_map_path = (
            "Processing/DisplacementMapPerSegment/fast_symmetric_demons/Artery"
        )
        self.assertEqual(24, len(outputs))
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
            "Processing/DisplacementProfiles/level_set_motion/Artery/"
            "X_sum_displacement_profile/value",
            outputs,
        )
        self.assertIn(
            "Processing/DisplacementProfiles/level_set_motion/Vein/"
            "Y_sum_displacement_profile/value",
            outputs,
        )

        magnitude_profile_path = (
            "Processing/DisplacementProfiles/level_set_motion/Vein/"
            "Magnitude_displacement_profile/value"
        )
        amplitude_profile_path = (
            "Processing/DisplacementProfiles/level_set_motion/Vein/"
            "Cross_sectional_radial_movement_amplitude_profile/value"
        )
        asymmetry_profile_path = (
            "Processing/DisplacementProfiles/level_set_motion/Vein/"
            "Cross_sectional_radial_asymmetry_index_profile/value"
        )
        self.assertIn(magnitude_profile_path, outputs)
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

            x_sum_profile = h5[x_sum_profile_path]
            y_sum_profile = h5[y_sum_profile_path]
            displacement_map = h5[displacement_map_path]
            magnitude_profile = h5[magnitude_profile_path]
            amplitude_profile = h5[amplitude_profile_path]
            asymmetry_profile = h5[asymmetry_profile_path]

            self.assertEqual((4, 2, 1, 1), x_sum_profile.shape)
            self.assertEqual((4, 2, 1, 1), y_sum_profile.shape)
            self.assertEqual((4, 3, 4, 2, 1, 1, 2), displacement_map.shape)
            self.assertEqual((4, 2), magnitude_profile.shape)
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
            self.assertEqual("pixels", magnitude_profile.attrs["unit"])
            self.assertEqual(
                ["time", "beat"],
                list(magnitude_profile.attrs["dimDesc"]),
            )
            self.assertEqual(
                "sum_of_segment_vector_magnitudes",
                magnitude_profile.attrs["spatial_reduction"],
            )
            np.testing.assert_allclose(
                magnitude_profile[...],
                np.hypot(24.0, -48.0),
                atol=1e-5,
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

    def test_combined_magnitude_adds_segments_without_directional_cancellation(
        self,
    ) -> None:
        x_sums = np.asarray([[[3.0] * 6, [-3.0] * 6]], dtype=np.float32)
        y_sums = np.asarray([[[4.0] * 6, [-4.0] * 6]], dtype=np.float32)

        dataset = _combined_displacement_magnitude_dataset(
            x_sums,
            y_sums,
            np.asarray([0, 2, 5], dtype=np.int32),
            index_base=0,
        )

        self.assertEqual((4, 2), dataset.data.shape)
        np.testing.assert_allclose(dataset.data, 10.0, atol=1e-6)

    def test_displacement_magnitude_is_a_per_segment_per_beat_trace(
        self,
    ) -> None:
        shape = (2, 3, 131)
        displacement = SimpleNamespace(
            x_sum_displacement_profile=np.full(
                shape,
                3.0,
                dtype=np.float32,
            ),
            y_sum_displacement_profile=np.full(
                shape,
                4.0,
                dtype=np.float32,
            ),
        )
        segments = SimpleNamespace(
            displacements={"level_set_motion": displacement},
        )

        outputs = pack_displacement_magnitude_outputs(
            segments,
            segments,
            np.asarray([0, 65, 130], dtype=np.int32),
        )
        artery_path = (
            "Processing/DisplacementProfiles/level_set_motion/Artery/"
            "displacement_magnitude"
        )
        vein_path = (
            "Processing/DisplacementProfiles/level_set_motion/Vein/"
            "displacement_magnitude"
        )

        self.assertEqual({artery_path, vein_path}, set(outputs))
        with h5py.File(
            "displacement_magnitude.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            for path in (artery_path, vein_path):
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
        artery_root = (
            "Processing/DisplacementProfiles/level_set_motion/Artery"
        )
        vein_root = "Processing/DisplacementProfiles/level_set_motion/Vein"
        artery_longitudinal = f"{artery_root}/Longitudinal"
        artery_transverse = f"{artery_root}/Transverse"
        vein_longitudinal = f"{vein_root}/Longitudinal"
        vein_transverse = f"{vein_root}/Transverse"
        profile_paths = {
            f"{artery_longitudinal}/LongitudinalDisplacementProfileMasked",
            f"{artery_longitudinal}/LongitudinalDisplacementProfileUnmasked",
            f"{artery_transverse}/TransverseDisplacementProfileMasked",
            f"{artery_transverse}/TransverseDisplacementProfileUnmasked",
            f"{vein_longitudinal}/LongitudinalDisplacementProfileMasked",
            f"{vein_longitudinal}/LongitudinalDisplacementProfileUnmasked",
        }
        meaned_paths = {
            f"{artery_longitudinal}/LongitudinalDisplacementProfileMaskedMeaned",
            f"{artery_transverse}/TransverseDisplacementProfileMaskedMeaned",
            f"{vein_transverse}/TransverseDisplacementProfileMaskedMeaned",
        }
        global_meaned_paths = {
            f"{artery_longitudinal}/"
            "LongitudinalDisplacementProfileMaskedGlobalMeaned",
            f"{artery_transverse}/"
            "TransverseDisplacementProfileMaskedGlobalMeaned",
        }
        mean_power_paths = {
            f"{artery_transverse}/Mean_P_D_transverse",
            f"{vein_transverse}/Mean_P_D_transverse",
        }
        gaussian_profile_paths = {
            f"{artery_transverse}/Gaussian_Fit",
            f"{vein_transverse}/Gaussian_Fit",
        }
        power_paths = {
            f"{artery_longitudinal}/P_D_longitudinal",
            f"{artery_transverse}/P_D_transverse",
        }
        metric_names = {
            "Max_X_Position",
            "Max_Y_Position",
            "Diff_Y_Value",
            "Mean_X_Position",
            "Mean_Y_Position",
            "Mean_Diff_Y_Value",
            "Area_L",
            "Area_R",
            "Diff_Area_L",
            "Diff_Area_R",
            "Gaussian_Baseline",
            "Gaussian_A_L",
            "Gaussian_A_R",
            "Gaussian_Mu_L",
            "Gaussian_Mu_R",
            "Gaussian_Sigma_L",
            "Gaussian_Sigma_R",
            "Gaussian_FWHM_L",
            "Gaussian_FWHM_R",
            "Gaussian_Area_L",
            "Gaussian_Area_R",
            "Gaussian_Peak_Separation",
            "Gaussian_RMSE",
            "Gaussian_Fit_Success",
            "Gaussian_Initialization_Complete",
        }
        metric_paths = {
            "Processing/DisplacementMetrics/level_set_motion/"
            f"{vessel_name}/Transverse/{metric_name}"
            for vessel_name in ("Artery", "Vein")
            for metric_name in metric_names
        }
        level_set_paths = (
            profile_paths
            | meaned_paths
            | global_meaned_paths
            | mean_power_paths
            | gaussian_profile_paths
            | power_paths
        )
        other_method_paths = {
            path.replace(
                "/level_set_motion/",
                "/fast_symmetric_demons/",
            )
            for path in level_set_paths
        }
        other_method_metric_paths = {
            path.replace(
                "/level_set_motion/",
                "/fast_symmetric_demons/",
            )
            for path in metric_paths
        }
        expected_paths = (
            level_set_paths
            | other_method_paths
            | metric_paths
            | other_method_metric_paths
        )
        self.assertEqual(expected_paths, set(outputs))

        with h5py.File(
            "displacement_axis_profiles.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            for path in profile_paths | {
                path.replace(
                    "/level_set_motion/",
                    "/fast_symmetric_demons/",
                )
                for path in profile_paths
            }:
                dataset = h5[path]
                self.assertEqual((181, 4, 2, 1, 1), dataset.shape)
                self.assertEqual("pixels", dataset.attrs["unit"])
                expected_axis = (
                    "x" if "Transverse" in path else "y"
                )
                self.assertEqual(
                    [expected_axis, "time", "beat", "branch", "radius"],
                    list(dataset.attrs["dimDesc"]),
                )
            transverse_unmasked = h5[
                f"{artery_transverse}/TransverseDisplacementProfileUnmasked"
            ][...]
            np.testing.assert_allclose(
                h5[
                    f"{artery_transverse}/TransverseDisplacementProfileMasked"
                ][...],
                2.0 * transverse_unmasked,
                atol=1e-5,
            )
            np.testing.assert_allclose(
                h5[
                    f"{artery_longitudinal}/"
                    "LongitudinalDisplacementProfileUnmasked"
                ][...],
                3.0 * transverse_unmasked,
                atol=1e-5,
            )
            np.testing.assert_allclose(
                h5[
                    f"{artery_longitudinal}/"
                    "LongitudinalDisplacementProfileMasked"
                ][...],
                4.0 * transverse_unmasked,
                atol=1e-5,
            )
            for profile_name in (
                "LongitudinalDisplacementProfileMasked",
                "TransverseDisplacementProfileMasked",
            ):
                direction = (
                    "transverse"
                    if "Transverse" in profile_name
                    else "longitudinal"
                )
                direction_root = (
                    artery_transverse
                    if direction == "transverse"
                    else artery_longitudinal
                )
                source = h5[f"{direction_root}/{profile_name}"]
                meaned = h5[f"{direction_root}/{profile_name}Meaned"]
                self.assertEqual((181, 2, 1, 1), meaned.shape)
                expected_axis = (
                    "x" if "Transverse" in profile_name else "y"
                )
                self.assertEqual(
                    [expected_axis, "beat", "branch", "radius"],
                    list(meaned.attrs["dimDesc"]),
                )
                self.assertEqual(
                    "mean_over_interpolated_beat_time",
                    meaned.attrs["temporal_reduction"],
                )
                np.testing.assert_allclose(
                    meaned[...],
                    np.nanmean(source[...], axis=1),
                    atol=1e-6,
                )
                global_meaned = h5[
                    f"{direction_root}/{profile_name}GlobalMeaned"
                ]
                self.assertEqual((181, 2), global_meaned.shape)
                self.assertEqual(
                    [expected_axis, "beat"],
                    list(global_meaned.attrs["dimDesc"]),
                )
                self.assertEqual(
                    "mean_over_valid_branch_radius_segments",
                    global_meaned.attrs["segment_reduction"],
                )
                np.testing.assert_allclose(
                    global_meaned[...],
                    np.nanmean(meaned[...], axis=(2, 3)),
                    atol=1e-6,
                )
                power = h5[f"{direction_root}/P_D_{direction}"]
                self.assertEqual(source.shape, power.shape)
                self.assertEqual("pixels^2", power.attrs["unit"])
                self.assertEqual(
                    ["x" if direction == "transverse" else "y",
                     "time", "beat", "branch", "radius"],
                    list(power.attrs["dimDesc"]),
                )
                self.assertEqual(
                    "(D(t) - mean_t(D(t)))**2",
                    power.attrs["formula"],
                )
                np.testing.assert_allclose(
                    power[...],
                    (source[...] - meaned[...][:, None, ...]) ** 2,
                    atol=1e-6,
                )
                self.assertGreater(float(np.nanmax(power[...])), 0.0)
                if direction == "transverse":
                    mean_power = h5[
                        f"{direction_root}/Mean_P_D_transverse"
                    ]
                    self.assertEqual((181, 2, 1, 1), mean_power.shape)
                    self.assertEqual(
                        ["x", "beat", "branch", "radius"],
                        list(mean_power.attrs["dimDesc"]),
                    )
                    np.testing.assert_allclose(
                        mean_power[...],
                        np.nanmean(power[...], axis=1),
                        atol=1e-6,
                    )

            np.testing.assert_allclose(
                h5[f"{vein_root}/Transverse/Mean_P_D_transverse"][...],
                h5[f"{artery_transverse}/Mean_P_D_transverse"][...],
                atol=1e-6,
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
            "transverse_displacement_profiles_unmasked",
            "transverse_displacement_profiles_masked",
            "longitudinal_displacement_profiles_unmasked",
            "longitudinal_displacement_profiles_masked",
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
            root = (
                "Processing/DisplacementProfiles/level_set_motion/Artery/"
                f"{direction}/{direction}DisplacementProfileMasked"
            )
            source = outputs[f"{root}Meaned"]
            global_meaned = outputs[f"{root}GlobalMeaned"]
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
            metrics_root = (
                "Processing/DisplacementMetrics/level_set_motion/"
                f"{vessel_name}/Transverse"
            )
            max_x = outputs[f"{metrics_root}/Max_X_Position"]
            max_y = outputs[f"{metrics_root}/Max_Y_Position"]
            diff_y = outputs[f"{metrics_root}/Diff_Y_Value"]
            mean_x = outputs[f"{metrics_root}/Mean_X_Position"]
            mean_y = outputs[f"{metrics_root}/Mean_Y_Position"]
            mean_diff_y = outputs[f"{metrics_root}/Mean_Diff_Y_Value"]
            area_l = outputs[f"{metrics_root}/Area_L"]
            area_r = outputs[f"{metrics_root}/Area_R"]
            diff_area_l = outputs[f"{metrics_root}/Diff_Area_L"]
            diff_area_r = outputs[f"{metrics_root}/Diff_Area_R"]

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
            "Processing/DisplacementMetrics/level_set_motion/"
            "Artery/Transverse"
        )
        transverse_root = (
            "Processing/DisplacementProfiles/level_set_motion/"
            "Artery/Transverse"
        )
        max_x = outputs[f"{artery_metrics_root}/Max_X_Position"].data
        max_y = outputs[f"{artery_metrics_root}/Max_Y_Position"].data
        diff_y = outputs[f"{artery_metrics_root}/Diff_Y_Value"].data
        meaned = outputs[
            f"{transverse_root}/TransverseDisplacementProfileMaskedMeaned"
        ].data
        mean_power = outputs[f"{transverse_root}/Mean_P_D_transverse"].data
        area_l = outputs[f"{artery_metrics_root}/Area_L"].data
        area_r = outputs[f"{artery_metrics_root}/Area_R"].data
        diff_area_l = outputs[f"{artery_metrics_root}/Diff_Area_L"].data
        diff_area_r = outputs[f"{artery_metrics_root}/Diff_Area_R"].data
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

    def test_displacement_map_loading_is_optional_and_method_aware(self) -> None:
        unscheduled = SimpleNamespace(
            pipeline_scheduled=lambda name: False,
        )
        self.assertEqual({}, _load_displacement_maps(unscheduled))

        with tempfile.TemporaryDirectory() as temp_dir:
            artery_path = Path(temp_dir) / "artery.npy"
            vein_path = Path(temp_dir) / "vein.npy"
            np.save(artery_path, np.zeros((2, 3, 4, 2), dtype=np.float32))
            np.save(vein_path, np.ones((2, 3, 4, 2), dtype=np.float32))
            artifacts = DisplacementMapArtifacts(
                registration_method="diffeomorphic_demons",
                field_paths_by_vessel={
                    "artery": artery_path,
                    "vein": vein_path,
                },
                temporary_directory=SimpleNamespace(cleanup=lambda: None),
            )
            scheduled = SimpleNamespace(
                pipeline_scheduled=lambda name: name == "displacement_map",
                state=SimpleNamespace(
                    get=lambda key: (
                        artifacts if key == DISPLACEMENT_MAP_STATE else None
                    )
                ),
            )
            loaded = _load_displacement_maps(scheduled)
            try:
                self.assertEqual({"artery", "vein"}, set(loaded))
                self.assertEqual(
                    ["diffeomorphic_demons"],
                    list(loaded["artery"]),
                )
                np.testing.assert_array_equal(
                    loaded["artery"]["diffeomorphic_demons"],
                    0.0,
                )
                np.testing.assert_array_equal(
                    loaded["vein"]["diffeomorphic_demons"],
                    1.0,
                )
            finally:
                for maps_for_vessel in loaded.values():
                    for displacement_map in maps_for_vessel.values():
                        displacement_map._mmap.close()


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
            displacement_maps_per_segment=vector_maps * scale,
            transverse_displacement_profiles_unmasked=(
                profile_time * scale
            ),
            transverse_displacement_profiles_masked=(
                profile_time * np.float32(2.0) * scale
            ),
            longitudinal_displacement_profiles_unmasked=(
                profile_time * np.float32(3.0) * scale
            ),
            longitudinal_displacement_profiles_masked=(
                profile_time * np.float32(4.0) * scale
            ),
            x_sum_displacement_profile=x_sum_profile * scale,
            y_sum_displacement_profile=y_sum_profile * scale,
            cross_sectional_radial_movement_amplitude=(
                radial_amplitude * scale
            ),
            cross_sectional_radial_asymmetry_index=radial_asymmetry,
        )

    return SimpleNamespace(
        displacements={
            "fast_symmetric_demons": result(1.0),
            "level_set_motion": result(2.0),
        }
    )


if __name__ == "__main__":
    unittest.main()
