"""Tests for transverse profiles of per-segment velocity FFT magnitudes."""

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

from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5 import write_value_dataset  # noqa: E402
from pipelines.waveform_velocity.profiles import (  # noqa: E402
    pack_velocity_profile_fft_outputs,
    velocity_fft_transverse_profiles,
)
from pipelines.waveform_velocity.segment_maps import (  # noqa: E402
    interpolate_velocity_maps_per_beat,
)


class VelocityFFTProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        y_values = np.arange(51, dtype=np.float32)
        maps = np.broadcast_to(
            y_values[None, None, None, :, None],
            (1, 1, 6, 51, 2),
        ).copy()
        masks = np.zeros((1, 1, 51, 2), dtype=bool)
        masks[0, 0, 10, :] = True
        self.segments = SimpleNamespace(
            velocity_maps_per_segment=maps,
            segment_masks=masks,
        )
        self.boundaries = np.asarray([0, 2, 5], dtype=np.int32)
        self.maps_per_beat = interpolate_velocity_maps_per_beat(
            maps,
            self.boundaries,
        )

    def test_fft_then_nanmean_uses_dilated_mask(self) -> None:
        unmasked, masked = velocity_fft_transverse_profiles(
            self.maps_per_beat,
            self.segments.segment_masks,
        )

        self.assertEqual((2, 4, 2, 1, 1), unmasked.shape)
        self.assertEqual(unmasked.shape, masked.shape)
        np.testing.assert_allclose(unmasked[:, 0, :, 0, 0], 100.0, atol=1e-5)
        np.testing.assert_allclose(masked[:, 0, :, 0, 0], 60.0, atol=1e-5)
        np.testing.assert_allclose(unmasked[:, 1:, :, 0, 0], 0.0, atol=1e-5)
        np.testing.assert_allclose(masked[:, 1:, :, 0, 0], 0.0, atol=1e-5)

    def test_packer_writes_exact_artery_paths_and_metadata(self) -> None:
        outputs = pack_velocity_profile_fft_outputs(
            self.segments,
            self.segments,
            self.maps_per_beat,
            self.maps_per_beat,
        )
        schema = EyeFlowOutputPaths.active()
        paths = schema.artery_velocity_profiles
        expected_paths = {
            "Processing/VelocityProfilesFFT/Artery/"
            "TransverseVelocityProfileUnmasked",
            "Processing/VelocityProfilesFFT/Artery/"
            "TransverseVelocityProfileMasked",
            "Processing/VelocityProfilesFFT/Vein/"
            "TransverseVelocityProfileUnmasked",
            "Processing/VelocityProfilesFFT/Vein/"
            "TransverseVelocityProfileMasked",
        }

        self.assertEqual(expected_paths, set(outputs))
        self.assertEqual(
            "Processing/VelocityProfilesFFT/Artery/"
            "TransverseVelocityProfileUnmasked",
            paths.transverse_velocity_profile_fft_unmasked,
        )
        self.assertEqual(
            "Processing/VelocityProfilesFFT/Artery/"
            "TransverseVelocityProfileMasked",
            paths.transverse_velocity_profile_fft_masked,
        )
        self.assertEqual(
            "Processing/VelocityProfilesFFT/Vein/"
            "TransverseVelocityProfileUnmasked",
            schema.vein_velocity_profiles.transverse_velocity_profile_fft_unmasked,
        )
        self.assertEqual(
            "Processing/VelocityProfilesFFT/Vein/TransverseVelocityProfileMasked",
            schema.vein_velocity_profiles.transverse_velocity_profile_fft_masked,
        )

        with h5py.File(
            "velocity-fft-profiles.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5:
            for path, value in outputs.items():
                write_value_dataset(h5, path, value)

            unmasked = h5[paths.transverse_velocity_profile_fft_unmasked]
            masked = h5[paths.transverse_velocity_profile_fft_masked]
            self.assertEqual((2, 4, 2, 1, 1), unmasked.shape)
            self.assertEqual(
                ["x", "frequency", "beat", "branch", "radius"],
                list(masked.attrs["dimDesc"]),
            )
            self.assertEqual("full", masked.attrs["fft_spectrum"])
            self.assertEqual(20, masked.attrs["mask_dilation_iterations"])
            self.assertEqual(0, unmasked.attrs["mask_dilation_iterations"])
            self.assertEqual("gzip", masked.compression)


if __name__ == "__main__":
    unittest.main()
