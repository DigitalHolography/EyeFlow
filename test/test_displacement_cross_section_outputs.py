"""Tests for method-scoped displacement cross-section HDF5 outputs."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.writers.h5 import write_value_dataset  # noqa: E402
from pipelines.waveform_velocity.profiles import (  # noqa: E402
    pack_displacement_profile_outputs,
)
from pipelines.waveform_velocity.segment_maps import (  # noqa: E402
    pack_displacement_segment_map_outputs,
)
from pipelines.waveform_velocity_core.runner import (  # noqa: E402
    _load_displacement_maps,
)


class DisplacementOutputTests(unittest.TestCase):
    def test_method_scoped_scalar_and_debug_vector_outputs_reach_h5(self) -> None:
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

        scalar_profile_path = (
            "Processing/DisplacementProfiles/fast_symmetric_demons/Artery/"
            "TransverseDisplacementProfileUnmasked/value"
        )
        vector_profile_path = (
            "Processing/Debug/DisplacementProfiles/fast_symmetric_demons/"
            "Artery/TransverseDisplacementProfileUnmasked/value"
        )
        scalar_map_path = (
            "Processing/DisplacementMapPerSegment/fast_symmetric_demons/Artery"
        )
        vector_map_path = (
            "Processing/Debug/DisplacementMapPerSegment/"
            "fast_symmetric_demons/Artery"
        )
        self.assertEqual(40, len(outputs))

        with h5py.File(
            "displacement_outputs.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            scalar_profile = h5[scalar_profile_path]
            vector_profile = h5[vector_profile_path]
            scalar_map = h5[scalar_map_path]
            vector_map = h5[vector_map_path]

            self.assertEqual((4, 4, 2, 1, 1), scalar_profile.shape)
            self.assertEqual((4, 4, 2, 1, 1, 2), vector_profile.shape)
            self.assertEqual((4, 3, 4, 2, 1, 1), scalar_map.shape)
            self.assertEqual((4, 3, 4, 2, 1, 1, 2), vector_map.shape)
            self.assertEqual("pixels", scalar_profile.attrs["unit"])
            self.assertEqual("pixels", scalar_map.attrs["unit"])
            self.assertEqual(
                ["local_x", "local_y"],
                list(vector_profile.attrs["components"]),
            )
            self.assertEqual(
                "signed_local_y",
                scalar_map.attrs["component"],
            )
            self.assertLess(float(np.nanmax(scalar_profile[...])), 0.0)
            self.assertLess(float(np.nanmax(scalar_map[...])), 0.0)

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

        with h5py.File(
            "displacement_input.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            dataset = h5.create_dataset(
                "Processing/DisplacementMap",
                data=np.zeros((2, 3, 4, 2), dtype=np.float32),
            )
            dataset.attrs["registration_method"] = "diffeomorphic_demons"
            scheduled = SimpleNamespace(
                pipeline_scheduled=lambda name: name == "displacement_map",
                runtime=SimpleNamespace(work_h5=h5),
            )
            loaded = _load_displacement_maps(scheduled)

            self.assertEqual(["diffeomorphic_demons"], list(loaded))
            self.assertEqual(
                dataset.name,
                loaded["diffeomorphic_demons"].name,
            )


def _segments():
    scalar_profiles = np.full((1, 1, 6, 4), -2.0, dtype=np.float32)
    vector_profiles = np.stack(
        (
            np.full_like(scalar_profiles, 1.0),
            scalar_profiles,
        ),
        axis=-1,
    )
    scalar_maps = np.full((1, 1, 6, 3, 4), -2.0, dtype=np.float32)
    vector_maps = np.stack(
        (
            np.full_like(scalar_maps, 1.0),
            scalar_maps,
        ),
        axis=-1,
    )

    def result(scale: float):
        return SimpleNamespace(
            displacement_maps_per_segment=scalar_maps * scale,
            displacement_vectors_per_segment=vector_maps * scale,
            displacement_profiles=scalar_profiles * scale,
            transverse_displacement_profiles_masked=scalar_profiles * scale,
            longitudinal_displacement_profiles_unmasked=scalar_profiles * scale,
            longitudinal_displacement_profiles_masked=scalar_profiles * scale,
            displacement_vector_profiles=vector_profiles * scale,
            transverse_displacement_vector_profiles_masked=(
                vector_profiles * scale
            ),
            longitudinal_displacement_vector_profiles_unmasked=(
                vector_profiles * scale
            ),
            longitudinal_displacement_vector_profiles_masked=(
                vector_profiles * scale
            ),
        )

    return SimpleNamespace(
        displacements={
            "fast_symmetric_demons": result(1.0),
            "level_set_motion": result(2.0),
        }
    )


if __name__ == "__main__":
    unittest.main()
