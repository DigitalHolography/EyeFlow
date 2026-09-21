"""Tests for annular vessel profiles of the moment0ff spatial gradient."""

from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from calculations.blood_flow_velocity import CrossSectionSignalSettings
from pipelines.waveform_velocity import spatial_gradient_profiles as profile_module
from pipelines.waveform_velocity.spatial_gradient_profiles import (
    SPATIAL_GRADIENT_METRICS_ROOT,
    SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES,
    SPATIAL_GRADIENT_PROFILE_ROOT,
    extract_spatial_gradient_segments,
    pack_spatial_gradient_profile_outputs,
)


class SpatialGradientProfileTests(unittest.TestCase):
    def test_extracts_both_vessels_using_annular_cross_section_engine(self) -> None:
        source = SimpleNamespace(
            optic_disc_mask=None,
            optic_disc_width=2.0,
            optic_disc_height=2.0,
            optic_disc_center=np.asarray([4.0, 4.0]),
            retinal_artery_mask=np.ones((8, 8), dtype=bool),
            retinal_vein_mask=np.eye(8, dtype=bool),
            cross_section_settings=CrossSectionSignalSettings(0.01),
        )
        moment0ff = np.full((3, 8, 8), 42.0, dtype=np.float32)
        ctx = SimpleNamespace(
            inputs=SimpleNamespace(
                hd=SimpleNamespace(
                    as_holodoppler=lambda: SimpleNamespace(
                        moment0_flat_field_dataset=lambda: moment0ff
                    )
                )
            ),
        )
        waveform_context = SimpleNamespace(
            source_data=source,
            attrs={"number_of_radii_in_FOV": 4},
            artery_segment_result=SimpleNamespace(
                topology=SimpleNamespace(prepared_topology="artery topology")
            ),
            vein_segment_result=SimpleNamespace(
                topology=SimpleNamespace(prepared_topology="vein topology")
            ),
        )

        with patch.object(
            profile_module,
            "analyze_velocity_segments",
            return_value={"artery": "artery", "vein": "vein"},
        ) as extract:
            result = extract_spatial_gradient_segments(ctx, waveform_context)

        self.assertEqual(("artery", "vein"), result)
        args, kwargs = extract.call_args
        self.assertEqual((3, 8, 8), args[0].shape)
        self.assertIs(moment0ff, args[0])
        self.assertIs(source.cross_section_settings, args[4])
        np.testing.assert_array_equal(source.retinal_artery_mask, args[1]["artery"])
        np.testing.assert_array_equal(source.retinal_vein_mask, args[1]["vein"])
        disc = kwargs["optic_disc_mask"]
        self.assertEqual((8, 8), disc.shape)
        self.assertTrue(disc[4, 4])
        self.assertFalse(disc[0, 0])
        self.assertEqual(
            5,
            kwargs["transverse_mask_dilation_pixels"],
        )
        self.assertEqual("staged", kwargs["transform_mode"])
        self.assertEqual(6, kwargs["temporal_halo"])
        self.assertEqual(
            {"artery": "artery topology", "vein": "vein topology"},
            kwargs["prepared_topologies"],
        )
        self.assertFalse(kwargs["retain_displacement_maps"])

    def test_packs_requested_profiles_for_arteries_and_veins(self) -> None:
        unmasked = np.arange(30, dtype=np.float32).reshape(2, 1, 5, 3)
        masked = unmasked.copy()
        masked[..., 0] = np.nan
        segments = SimpleNamespace(
            velocity_profiles=unmasked,
            transverse_velocity_profiles_masked=masked,
        )

        outputs = pack_spatial_gradient_profile_outputs(
            segments,
            segments,
            np.asarray([0, 2, 4], dtype=np.int32),
        )

        expected_paths = set()
        for vessel in ("Artery", "Vein"):
            root = f"{SPATIAL_GRADIENT_PROFILE_ROOT}/{vessel}/Transverse"
            metrics_root = (
                f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel}/Transverse"
            )
            expected_paths.update(
                {
                    f"{root}/Masked/SpatialGradientProfile/value",
                    f"{root}/Masked/SpatialGradientProfileMeaned/value",
                    f"{root}/Unmasked/SpatialGradientProfile/value",
                }
            )
            for mask_name in ("Masked", "Unmasked"):
                for dimension_tag in ("tbkr", "bkr", "bk", "kr", "k"):
                    expected_paths.update(
                        {
                            f"{metrics_root}/{mask_name}/{dimension_tag}/left_edge_index",
                            f"{metrics_root}/{mask_name}/{dimension_tag}/right_edge_index",
                            f"{metrics_root}/{mask_name}/{dimension_tag}/lumen/size",
                            f"{metrics_root}/{mask_name}/{dimension_tag}/lumen/size_qc",
                        }
                    )
            expected_paths.update(
                {
                    f"{metrics_root}/Masked/tbkr/lumen/size_median",
                    f"{metrics_root}/Masked/tbkr/lumen/size_std",
                    f"{metrics_root}/Masked/tk/lumen_size",
                }
            )
        self.assertEqual(expected_paths, set(outputs))
        for path, value in outputs.items():
            if "SpatialGradientProfiles" in path:
                self.assertEqual("/moment0ff", value.attrs["source_dataset"])
                self.assertEqual(7, value.attrs["temporal_moving_average_window"])
                self.assertEqual(2, value.attrs["temporal_moving_average_passes"])
                self.assertEqual(
                    "centered_pixelwise_moving_average", value.attrs["temporal_filter"]
                )
                self.assertEqual("truncated_window", value.attrs["temporal_boundary_mode"])
                self.assertEqual("propagate", value.attrs["temporal_nan_policy"])
                self.assertEqual("ImageJ", value.attrs["gaussian_blur_algorithm"])
                self.assertEqual(6.0, value.attrs["gaussian_blur_radius_pixels"])
                self.assertEqual("propagate", value.attrs["gaussian_blur_nan_policy"])
                self.assertEqual(
                    "nearest_extension_zero_output_border",
                    value.attrs["gaussian_blur_boundary_mode"],
                )
                self.assertEqual("ImageJ", value.attrs["unsharp_mask_algorithm"])
                self.assertEqual(8.0, value.attrs["unsharp_mask_radius_pixels"])
                self.assertAlmostEqual(0.6, value.attrs["unsharp_mask_weight"])
                self.assertEqual("propagate", value.attrs["unsharp_mask_nan_policy"])
                self.assertEqual("clip_to_zero", value.attrs["unsharp_mask_negative_values"])
                self.assertEqual(
                    "spatial_interpolation, moving_average, sobel_magnitude, "
                    "gaussian2d_blur, unsharp_mask, moving_average, segment_rotation",
                    value.attrs["processing_order"],
                )
                self.assertEqual(
                    "3x3 Sobel magnitude",
                    value.attrs["spatial_operator"],
                )
                if "/SpatialGradientProfileMeaned/" in path:
                    self.assertEqual(
                        ["x", "beat", "branch", "radius"],
                        value.attrs["dimDesc"],
                    )
                    self.assertEqual((3, 2, 1, 2), value.data.shape)
                else:
                    self.assertEqual(
                        ["x", "time", "beat", "branch", "radius"],
                        value.attrs["dimDesc"],
                    )
                    self.assertEqual((3, 2, 2, 1, 2), value.data.shape)

        metric_dimensions = {}
        for mask_name in ("Masked", "Unmasked"):
            for dimension_tag, dimensions in {
                "tbkr": ["time", "beat", "branch", "radius"],
                "bkr": ["beat", "branch", "radius"],
                "bk": ["beat", "branch"],
                "kr": ["branch", "radius"],
                "k": ["branch"],
            }.items():
                metric_dimensions.update(
                    {
                        f"{mask_name}/{dimension_tag}/left_edge_index": dimensions,
                        f"{mask_name}/{dimension_tag}/right_edge_index": dimensions,
                        f"{mask_name}/{dimension_tag}/lumen/size": dimensions,
                        f"{mask_name}/{dimension_tag}/lumen/size_qc": dimensions,
                    }
                )
        metric_dimensions.update(
            {
                "Masked/tbkr/lumen/size_median": [],
                "Masked/tbkr/lumen/size_std": [],
                "Masked/tk/lumen_size": ["time", "branch"],
            }
        )
        for vessel in ("Artery", "Vein"):
            metrics_root = f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel}/Transverse"
            for name, dimensions in metric_dimensions.items():
                self.assertEqual(
                    dimensions,
                    outputs[f"{metrics_root}/{name}"].attrs["dimDesc"],
                )
            for mask_name in ("Masked", "Unmasked"):
                kr_qc = outputs[
                    f"{metrics_root}/{mask_name}/kr/lumen/size_qc"
                ]
                self.assertEqual(np.float32, kr_qc.data.dtype)
                self.assertEqual((1, 2), kr_qc.data.shape)
                self.assertTrue(np.all((kr_qc.data >= 0.0) & (kr_qc.data <= 1.0)))
                self.assertEqual("fraction", kr_qc.attrs["unit"])
                self.assertEqual(
                    ["branch", "radius"],
                    kr_qc.attrs["dimDesc"],
                )

    def test_gradient_peak_metrics_select_highest_separated_values(self) -> None:
        masked_profile = np.asarray(
            [0.0, 1.0, 8.0, 7.5, 0.0, 0.0, 1.0, 6.0, 7.0, 0.0, 0.0],
            dtype=np.float32,
        )
        unmasked_profile = np.asarray(
            [9.0, 0.0, 1.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0],
            dtype=np.float32,
        )
        masked_profiles = np.broadcast_to(masked_profile, (1, 1, 5, 11)).copy()
        unmasked_profiles = np.broadcast_to(unmasked_profile, (1, 1, 5, 11)).copy()
        segments = SimpleNamespace(
            velocity_profiles=unmasked_profiles,
            transverse_velocity_profiles_masked=masked_profiles,
        )

        outputs = pack_spatial_gradient_profile_outputs(
            segments,
            segments,
            np.asarray([0, 2, 4], dtype=np.int32),
        )

        for vessel in ("Artery", "Vein"):
            root = f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel}/Transverse"
            masked_left_index = outputs[f"{root}/Masked/tbkr/left_edge_index"]
            masked_right_index = outputs[f"{root}/Masked/tbkr/right_edge_index"]
            np.testing.assert_allclose(
                masked_left_index.data[:, 0, 0, 0],
                [2.4333334, 2.4333334],
            )
            np.testing.assert_allclose(
                masked_right_index.data[:, 0, 0, 0],
                [7.625, 7.625],
            )
            self.assertEqual("pixels", masked_left_index.attrs["unit"])
            unmasked_left_index = outputs[
                f"{root}/Unmasked/tbkr/left_edge_index"
            ]
            unmasked_right_index = outputs[
                f"{root}/Unmasked/tbkr/right_edge_index"
            ]
            np.testing.assert_array_equal(
                unmasked_left_index.data[:, 0, 0, 0],
                [0, 0],
            )
            np.testing.assert_array_equal(
                unmasked_right_index.data[:, 0, 0, 0],
                [6, 6],
            )
            self.assertEqual(
                "Unmasked/SpatialGradientProfile/value",
                unmasked_left_index.attrs["source_profile"],
            )
            self.assertEqual(
                "right_edge_index - left_edge_index",
                outputs[f"{root}/Masked/tbkr/lumen/size"].attrs["definition"],
            )
            np.testing.assert_allclose(
                outputs[f"{root}/Masked/tbkr/lumen/size"].data[:, 0, 0, 0],
                [5.1916666, 5.1916666],
            )
            self.assertEqual(
                SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES,
                masked_left_index.attrs["minimum_peak_gap_samples"],
            )

    def test_peak_indexes_use_quadratic_fractional_refinement(self) -> None:
        x = np.arange(7, dtype=np.float32)
        profile = -((x - np.float32(2.25)) ** 2)

        self.assertAlmostEqual(
            2.25,
            float(profile_module._fractional_peak_index(profile, 2)),
        )
        self.assertEqual(
            0.0,
            float(profile_module._fractional_peak_index(profile, 0)),
        )

        profile_with_gap = profile.copy()
        profile_with_gap[3] = np.nan
        self.assertEqual(
            2.0,
            float(profile_module._fractional_peak_index(profile_with_gap, 2)),
        )

        non_maximum = np.asarray([0.0, 1.0, 2.0], dtype=np.float32)
        self.assertEqual(
            1.0,
            float(profile_module._fractional_peak_index(non_maximum, 1)),
        )

    def test_edge_index_hierarchy_uses_nanmedian_in_requested_order(self) -> None:
        left_tbkr = np.asarray(
            [
                [
                    [[0.0, np.nan, 4.0], [3.0, 8.0, 2.0]],
                    [[9.0, 1.0, 5.0], [6.0, 7.0, 2.0]],
                ],
                [
                    [[2.0, 6.0, 8.0], [5.0, np.nan, 4.0]],
                    [[7.0, 3.0, 1.0], [8.0, 3.0, 9.0]],
                ],
            ],
            dtype=np.float32,
        )
        right_tbkr = left_tbkr + np.float32(10.0)

        metrics = profile_module._spatial_gradient_edge_index_metrics(
            left_tbkr,
            right_tbkr,
            mask_name="Masked",
            source_profile="profile",
            common_attrs={},
        )

        expected_hierarchies = {}
        for side, tbkr in (("left", left_tbkr), ("right", right_tbkr)):
            bkr = np.nanmedian(tbkr, axis=0)
            kr = np.nanmedian(bkr, axis=0)
            expected_hierarchies[side] = {
                "tbkr": tbkr,
                "bkr": bkr,
                "bk": np.nanmedian(bkr, axis=-1),
                "kr": kr,
                "k": np.nanmedian(kr, axis=-1),
            }

        for side, expected in expected_hierarchies.items():
            for dimension_tag, values in expected.items():
                edge_path = f"Masked/{dimension_tag}"
                np.testing.assert_array_equal(
                    metrics[f"{edge_path}/{side}_edge_index"].data,
                    values,
                )
                if side == "right":
                    lumen_size = values - expected_hierarchies["left"][dimension_tag]
                    np.testing.assert_array_equal(
                        metrics[f"{edge_path}/lumen/size"].data,
                        lumen_size,
                    )
                    if dimension_tag == "tbkr":
                        distribution_qc, _, _ = (
                            profile_module._lumen_size_standard_deviation_quality_control(
                                lumen_size
                            )
                        )
                        expected_qc = (
                            profile_module._tbkr_lumen_size_quality_control(
                                profile_module._kr_lumen_size_quality_control(
                                    distribution_qc
                                ),
                                lumen_size.shape,
                                threshold=(
                                    profile_module.TBKR_LUMEN_SIZE_QC_THRESHOLD
                                ),
                            )
                        )
                    elif dimension_tag == "kr":
                        tbkr_lumen_size = (
                            expected_hierarchies["right"]["tbkr"]
                            - expected_hierarchies["left"]["tbkr"]
                        )
                        tbkr_qc, _, _ = (
                            profile_module._lumen_size_standard_deviation_quality_control(
                                tbkr_lumen_size
                            )
                        )
                        expected_qc = (
                            profile_module._kr_lumen_size_quality_control(tbkr_qc)
                        )
                    else:
                        expected_qc, _, _ = (
                            profile_module._lumen_size_quality_control(lumen_size)
                        )
                    np.testing.assert_array_equal(
                        metrics[f"{edge_path}/lumen/size_qc"].data,
                        expected_qc,
                    )
        masked_lumen_size = right_tbkr - left_tbkr
        self.assertEqual(
            np.nanmedian(masked_lumen_size),
            metrics["Masked/tbkr/lumen/size_median"].data,
        )
        self.assertEqual(
            np.nanstd(masked_lumen_size),
            metrics["Masked/tbkr/lumen/size_std"].data,
        )

    def test_lumen_size_qc_uses_inclusive_percentiles_and_rejects_nan(self) -> None:
        lumen_size = np.concatenate(
            (
                np.arange(201, dtype=np.float32),
                np.asarray([np.nan], dtype=np.float32),
            )
        )

        qc, lower_limit, upper_limit = (
            profile_module._lumen_size_quality_control(lumen_size)
        )

        self.assertEqual(np.uint8, qc.dtype)
        self.assertEqual(1.0, lower_limit)
        self.assertEqual(199.0, upper_limit)
        self.assertEqual(0, qc[0])
        self.assertEqual(1, qc[1])
        self.assertEqual(1, qc[199])
        self.assertEqual(0, qc[200])
        self.assertEqual(0, qc[201])

    def test_lumen_size_distribution_statistics_ignore_nan(self) -> None:
        values = np.asarray([1.0, 2.0, 6.0, np.nan], dtype=np.float32)

        median, standard_deviation = (
            profile_module._lumen_size_distribution_statistics(values)
        )

        self.assertEqual(2.0, median)
        self.assertAlmostEqual(float(np.std([1.0, 2.0, 6.0])), standard_deviation)

        median, standard_deviation = (
            profile_module._lumen_size_distribution_statistics(
                np.asarray([np.nan], dtype=np.float32)
            )
        )
        self.assertTrue(np.isnan(median))
        self.assertTrue(np.isnan(standard_deviation))

    def test_tbkr_lumen_size_qc_uses_median_plus_or_minus_std(self) -> None:
        lumen_size = np.asarray(
            [0.0, 1.0, 2.0, 3.0, 4.0, 100.0, np.nan],
            dtype=np.float32,
        )

        qc, median, standard_deviation = (
            profile_module._lumen_size_standard_deviation_quality_control(
                lumen_size
            )
        )

        self.assertEqual(2.5, median)
        self.assertAlmostEqual(
            float(np.std([0.0, 1.0, 2.0, 3.0, 4.0, 100.0])),
            standard_deviation,
        )
        np.testing.assert_array_equal(qc, [1, 1, 1, 1, 1, 0, 0])

    def test_kr_lumen_size_qc_is_fractional_mean_of_tbkr_qc(self) -> None:
        tbkr_qc = np.asarray(
            [
                [
                    [[0, 1, 1], [1, 1, 0]],
                    [[0, 1, 1], [1, 0, 0]],
                ],
                [
                    [[1, 1, 1], [1, 0, 0]],
                    [[0, 0, 1], [0, 0, 0]],
                ],
            ],
            dtype=np.uint8,
        )

        kr_qc = profile_module._kr_lumen_size_quality_control(tbkr_qc)

        self.assertEqual(np.float32, kr_qc.dtype)
        self.assertEqual((2, 3), kr_qc.shape)
        np.testing.assert_array_equal(
            kr_qc,
            [[0.25, 0.75, 1.0], [0.75, 0.25, 0.0]],
        )

        empty = profile_module._kr_lumen_size_quality_control(
            np.empty((0, 2, 3, 4), dtype=np.uint8)
        )
        self.assertEqual(np.float32, empty.dtype)
        np.testing.assert_array_equal(
            empty,
            np.zeros((3, 4), dtype=np.float32),
        )

    def test_tbkr_lumen_size_qc_uses_half_threshold(self) -> None:
        kr_qc = np.asarray([[0.5, 0.6]], dtype=np.float32)

        actual = profile_module._tbkr_lumen_size_quality_control(
            kr_qc,
            (2, 1, 1, 2),
            threshold=profile_module.TBKR_LUMEN_SIZE_QC_THRESHOLD,
        )

        self.assertEqual(0.5, profile_module.TBKR_LUMEN_SIZE_QC_THRESHOLD)
        np.testing.assert_array_equal(
            actual,
            np.asarray([[[[0, 1]]], [[[0, 1]]]], dtype=np.uint8),
        )


if __name__ == "__main__":
    unittest.main()
