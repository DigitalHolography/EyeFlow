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
        self.assertEqual(12, len(outputs))
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

            self.assertEqual((4, 2, 1, 1), x_sum_profile.shape)
            self.assertEqual((4, 2, 1, 1), y_sum_profile.shape)
            self.assertEqual((4, 3, 4, 2, 1, 1, 2), displacement_map.shape)
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
    x_sum_profile = np.full((1, 1, 6), 12.0, dtype=np.float32)
    y_sum_profile = np.full((1, 1, 6), -24.0, dtype=np.float32)
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
            x_sum_displacement_profile=x_sum_profile * scale,
            y_sum_displacement_profile=y_sum_profile * scale,
        )

    return SimpleNamespace(
        displacements={
            "fast_symmetric_demons": result(1.0),
            "level_set_motion": result(2.0),
        }
    )


if __name__ == "__main__":
    unittest.main()
