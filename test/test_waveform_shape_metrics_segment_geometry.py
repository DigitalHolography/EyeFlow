"""Tests for spatial segment geometry used to build segment waveforms."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.vessel_segments.measurement import SegmentMeasurementSettings  # noqa: E402
from calculations.topology import OpticDisc, section_masks  # noqa: E402
from pipelines.velocity_analysis.analysis.segments import (  # noqa: E402
    analyze_velocity_segment_profiles,
)
class SegmentCenterTests(unittest.TestCase):
    def test_segment_velocity_results_runs_with_prepared_topology(self) -> None:
        velocity = np.ones((3, 61, 61), dtype=np.float32)
        artery_mask = np.zeros((61, 61), dtype=bool)
        artery_mask[27:34, 5:56] = True
        vein_mask = np.zeros_like(artery_mask)
        vein_mask[5:56, 37:44] = True
        optic_disc_mask = np.zeros_like(artery_mask)
        optic_disc_mask[27:34, 27:34] = True
        ring_settings = SimpleNamespace(
            inner_radius_frac=0.1,
            outer_radius_frac=0.7,
            ring_width_frac=0.25,
            ring_count=2,
            segment_length_frac=None,
        )
        settings = SegmentMeasurementSettings(0.01, 1.0)

        results = analyze_velocity_segment_profiles(
            velocity,
            {
                "artery": artery_mask,
                "vein": vein_mask,
            },
            OpticDisc(optic_disc_mask, (30.0, 30.0), None, None),
            ring_settings,
            settings,
            retain_velocity_maps=False,
            cycle_boundary_indexes=np.asarray([0, 2], dtype=np.int32),
            velocity_profile_fft=True,
        )

        for result in results.values():
            profile = result.profile
            self.assertGreater(profile.topology.native.branch_ids.size, 0)
            self.assertIsNone(profile.maps)
            fft = result.require_fft()
            self.assertEqual(
                (181, 2, 1, profile.topology.native.branch_ids.size, 2),
                fft.unmasked.shape,
            )
            self.assertEqual(
                fft.unmasked.shape,
                fft.masked.shape,
            )
            valid = profile.topology.valid_segments
            self.assertTrue(np.any(valid))
            np.testing.assert_allclose(profile.segment_signal[valid], 1.0)

    def test_annulus_starts_at_optic_disc_and_uses_fixed_fov_extent(self) -> None:
        center = np.asarray([267.0, 230.0])
        settings = OpticDisc(None, center, 55, 69).annulus_geometry((512, 512))
        mask = section_masks(
            (512, 512),
            center,
            settings,
        )[0]

        # The first annulus starts at ceil(min(width, height) / 2) = 28 px.
        self.assertFalse(mask[230, 294])
        self.assertTrue(mask[230, 295])
        self.assertFalse(mask[257, 267])
        self.assertTrue(mask[258, 267])

        # The first annulus has the configured 0.04 radial spacing.
        self.assertTrue(mask[230, 305])
        self.assertTrue(mask[230, 315])
        self.assertFalse(mask[230, 316])
        self.assertFalse(mask[279, 267])

if __name__ == "__main__":
    unittest.main()
