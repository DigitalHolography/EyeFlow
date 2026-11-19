"""Tests for transverse profiles of per-segment velocity FFT magnitudes."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.vessel_segments.measurement import MaskedArrays  # noqa: E402
from calculations.vessel_segments.profiles.fft import (  # noqa: E402
    SegmentFftAccumulator,
)
from calculations.vessel_segments.sampling.transforms import dilate_segment_masks  # noqa: E402
from pipelines.velocity_analysis.analysis.segment_maps import (  # noqa: E402
    interpolate_velocity_maps_per_beat,
)


def _retained_fft_reference(maps_per_beat, segment_masks):
    """Small full-array reference for the streamed FFT accumulator."""
    maps = np.asarray(maps_per_beat, dtype=np.float32)
    masks = dilate_segment_masks(
        segment_masks, iterations=10, horizontal_only=True,
    )
    shape = (maps.shape[0], maps.shape[2], maps.shape[3], maps.shape[4], maps.shape[5])
    unmasked = np.full(shape, np.nan, dtype=np.float32)
    masked = np.full_like(unmasked, np.nan)
    for radius in range(maps.shape[5]):
        for branch in range(maps.shape[4]):
            magnitude = np.abs(np.fft.fft(maps[..., branch, radius], axis=2)).astype(np.float32)
            mask = masks[radius, branch].T[:, :, None, None]
            for output, values in (
                (unmasked, magnitude),
                (masked, np.where(mask, magnitude, np.nan)),
            ):
                finite = np.isfinite(values)
                count = finite.sum(axis=1)
                total = np.where(finite, values, 0).sum(axis=1, dtype=np.float32)
                np.divide(
                    total, count,
                    out=output[..., branch, radius], where=count > 0,
                )
    return unmasked, masked


class VelocityFFTProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        y_values = np.arange(51, dtype=np.float32)
        maps = np.broadcast_to(
            y_values[None, None, None, :, None],
            (1, 1, 6, 51, 2),
        ).copy()
        masks = np.zeros((1, 1, 51, 2), dtype=bool)
        masks[0, 0, 10, :] = True
        self.boundaries = np.asarray([0, 2, 5], dtype=np.int32)
        self.maps_per_beat = interpolate_velocity_maps_per_beat(
            maps,
            self.boundaries,
        )
        accumulator = SegmentFftAccumulator(
            frame_count=maps.shape[2],
            ring_count=1,
            branch_count=1,
            canvas_side=maps.shape[-1],
            cycle_boundary_indexes=self.boundaries,
            index_base=0,
        )
        accumulator.observe(
            0,
            0,
            maps[0, 0],
            dilate_segment_masks(
                masks, iterations=10, horizontal_only=True
            )[0, 0],
        )
        self.segments = SimpleNamespace(
            segment_maps=maps,
            segment_masks=masks,
            transverse_fft=MaskedArrays(accumulator.unmasked, accumulator.masked),
        )

    def test_fft_then_nanmean_uses_dilated_mask(self) -> None:
        unmasked, masked = _retained_fft_reference(
            self.maps_per_beat,
            self.segments.segment_masks,
        )

        self.assertEqual((2, 4, 2, 1, 1), unmasked.shape)
        self.assertEqual(unmasked.shape, masked.shape)
        np.testing.assert_allclose(unmasked[:, 0, :, 0, 0], 100.0, atol=1e-5)
        np.testing.assert_allclose(masked[:, 0, :, 0, 0], 40.0, atol=1e-5)
        np.testing.assert_allclose(unmasked[:, 1:, :, 0, 0], 0.0, atol=1e-5)
        np.testing.assert_allclose(masked[:, 1:, :, 0, 0], 0.0, atol=1e-5)

    def test_streamed_accumulator_matches_legacy_for_temporal_data(self) -> None:
        rng = np.random.default_rng(3401)
        maps = rng.normal(size=(1, 1, 7, 51, 5)).astype(np.float32)
        maps[0, 0, :, :4, 0] = np.nan
        masks = np.zeros((1, 1, 51, 5), dtype=bool)
        masks[0, 0, 25, 1:4] = True
        boundaries = np.asarray([0, 3, 6], dtype=np.int32)
        maps_per_beat = interpolate_velocity_maps_per_beat(maps, boundaries)
        expected_unmasked, expected_masked = _retained_fft_reference(
            maps_per_beat,
            masks,
        )
        accumulator = SegmentFftAccumulator(
            frame_count=maps.shape[2],
            ring_count=1,
            branch_count=1,
            canvas_side=maps.shape[-1],
            cycle_boundary_indexes=boundaries,
            index_base=0,
        )

        accumulator.observe(
            0,
            0,
            maps[0, 0],
            dilate_segment_masks(
                masks, iterations=10, horizontal_only=True
            )[0, 0],
        )

        np.testing.assert_allclose(
            accumulator.unmasked,
            expected_unmasked,
            rtol=1e-6,
            atol=1e-6,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            accumulator.masked,
            expected_masked,
            rtol=1e-6,
            atol=1e-6,
            equal_nan=True,
        )

if __name__ == "__main__":
    unittest.main()
