"""Tests for annular vessel profiles of the moment0ff spatial gradient."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from pipelines.spatial_gradient_moment0.runner import (
    STATE_KEY,
    SpatialGradientMoment0Artifacts,
)
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
        temporary_directory = tempfile.TemporaryDirectory()
        gradient_path = Path(temporary_directory.name) / "gradient.npy"
        np.save(gradient_path, np.ones((3, 8, 8), dtype=np.float32))
        artifacts = SpatialGradientMoment0Artifacts(
            avi_path=Path("gradient.avi"),
            mean_png_path=Path("gradient.png"),
            gradient_path=gradient_path,
            frame_count=3,
            display_maximum=1.0,
            temporary_directory=temporary_directory,
        )
        source = SimpleNamespace(
            optic_disc_width=2.0,
            optic_disc_height=2.0,
            optic_disc_center=np.asarray([4.0, 4.0]),
            retinal_artery_mask=np.ones((8, 8), dtype=bool),
            retinal_vein_mask=np.eye(8, dtype=bool),
            cross_section_settings="settings",
        )
        ctx = SimpleNamespace(
            state=SimpleNamespace(
                get=lambda key: artifacts if key == STATE_KEY else None
            )
        )
        waveform_context = SimpleNamespace(
            source_data=source,
            attrs={"number_of_radii_in_FOV": 4},
        )

        with patch.object(
            profile_module,
            "segment_velocity_results",
            return_value=("artery", "vein"),
        ) as extract:
            result = extract_spatial_gradient_segments(ctx, waveform_context)

        self.assertEqual(("artery", "vein"), result)
        args, kwargs = extract.call_args
        self.assertEqual((3, 8, 8), args[0].shape)
        np.testing.assert_array_equal(source.retinal_artery_mask, args[1])
        np.testing.assert_array_equal(source.retinal_vein_mask, args[2])
        self.assertFalse(kwargs["retain_displacement_maps"])
        self.assertFalse(gradient_path.exists())

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
                    f"{root}/TransverseSpatialGradientProfileMasked",
                    f"{root}/TransverseSpatialGradientProfileMaskedMeaned",
                    f"{root}/TransverseSpatialGradientProfileUnmasked",
                    f"{metrics_root}/index_left_max",
                    f"{metrics_root}/index_right_max",
                    f"{metrics_root}/peak_value_left_max",
                    f"{metrics_root}/peak_value_right_max",
                    f"{metrics_root}/lumen_diameter",
                }
            )
        self.assertEqual(expected_paths, set(outputs))
        for path, value in outputs.items():
            if "SpatialGradientProfiles" not in path:
                self.assertEqual(
                    ["beat", "branch", "radius"],
                    value.attrs["dimDesc"],
                )
                self.assertEqual((2, 1, 2), value.data.shape)
                continue
            self.assertEqual("/moment0ff", value.attrs["source_dataset"])
            self.assertEqual("3x3 Sobel magnitude", value.attrs["spatial_operator"])
            if path.endswith("Meaned"):
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

    def test_gradient_peak_metrics_select_highest_separated_values(self) -> None:
        profile = np.asarray(
            [0.0, 1.0, 8.0, 7.5, 0.0, 0.0, 1.0, 6.0, 7.0, 0.0, 0.0],
            dtype=np.float32,
        )
        profiles = np.broadcast_to(profile, (1, 1, 5, 11)).copy()
        segments = SimpleNamespace(
            velocity_profiles=profiles,
            transverse_velocity_profiles_masked=profiles,
        )

        outputs = pack_spatial_gradient_profile_outputs(
            segments,
            segments,
            np.asarray([0, 2, 4], dtype=np.int32),
        )

        for vessel in ("Artery", "Vein"):
            root = f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel}/Transverse"
            left_index = outputs[f"{root}/index_left_max"]
            right_index = outputs[f"{root}/index_right_max"]
            left_value = outputs[f"{root}/peak_value_left_max"]
            right_value = outputs[f"{root}/peak_value_right_max"]
            lumen_diameter = outputs[f"{root}/lumen_diameter"]
            np.testing.assert_array_equal(left_index.data[:, 0, 0], [2, 2])
            np.testing.assert_array_equal(right_index.data[:, 0, 0], [8, 8])
            np.testing.assert_allclose(left_value.data[:, 0, 0], [8.0, 8.0])
            np.testing.assert_allclose(right_value.data[:, 0, 0], [7.0, 7.0])
            np.testing.assert_array_equal(lumen_diameter.data[:, 0, 0], [6, 6])
            self.assertEqual("pixels", left_index.attrs["unit"])
            self.assertEqual("a.u.", left_value.attrs["unit"])
            self.assertEqual("pixels", lumen_diameter.attrs["unit"])
            self.assertEqual(
                "index_right_max - index_left_max",
                lumen_diameter.attrs["definition"],
            )
            self.assertEqual(
                SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES,
                left_index.attrs["minimum_peak_gap_samples"],
            )


if __name__ == "__main__":
    unittest.main()
