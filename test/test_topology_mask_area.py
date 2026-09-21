"""Regression tests for reusable exact-annulus geometry and diagnostics."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from calculations.topology import SegmentRingSettings, circle_pixel_coverage
from calculations.topology.mask_area import annulus_widths_pixels
from pipelines.waveform_velocity.profiles import (
    pack_mask_detection_blood_volume_rate_outputs,
)


def _segments():
    labels = np.ones((41, 41), dtype=np.int32)
    yy, xx = np.indices(labels.shape, dtype=np.float32)
    radius_scale = np.hypot(20.0, 20.0)
    radius_sq = (yy - 20.0) ** 2 + (xx - 20.0) ** 2
    middle_radius = 0.25 * radius_scale
    outer_radius = 0.5 * radius_scale
    sections = np.asarray(
        [
            radius_sq <= middle_radius**2,
            (radius_sq > middle_radius**2) & (radius_sq <= outer_radius**2),
        ]
    )
    topology = SimpleNamespace(
        labels=labels,
        branch_ids=np.asarray([1], dtype=np.int32),
        section_masks=sections,
        ring_settings=SegmentRingSettings(0.0, 0.5, 0.25, 2, 0.25),
    )
    return SimpleNamespace(
        velocity=np.full((2, 1, 3), 2.0, dtype=np.float32),
        topology=topology,
    )


def test_circle_pixel_coverage_conserves_exact_area() -> None:
    radius = 12.4
    coverage = circle_pixel_coverage(
        (61, 57),
        cy=28.3,
        cx=27.7,
        radius_pixels=radius,
    )
    assert float(coverage.min()) >= 0.0
    assert float(coverage.max()) <= 1.0
    assert np.isclose(float(coverage.sum()), np.pi * radius**2, rtol=0, atol=1e-9)


def test_annulus_width_uses_clipped_last_ring() -> None:
    widths = annulus_widths_pixels(
        (7, 9),
        SegmentRingSettings(0.0, 0.6, 0.5, 2, 0.25),
        2,
    )
    np.testing.assert_allclose(widths, [1.25, 0.5], rtol=1e-6)


def test_mask_detection_circular_area_can_be_disabled_explicitly() -> None:
    segments = _segments()
    outputs = pack_mask_detection_blood_volume_rate_outputs(
        segments,
        segments,
        np.asarray([0, 2], dtype=np.int32),
        optic_disc_center=(20.0, 20.0),
        pixel_size_mm=0.1,
        apply_circular_area=False,
    )
    artery = outputs["Processing/BloodVolumeRate/Artery/maskDetection/value"]
    np.testing.assert_allclose(artery.data, 2.0)
    assert artery.attrs["unit"] == "mm/s"
    assert artery.attrs["quantity"] == "mean_velocity"
    assert artery.attrs["circular_area_scaling_applied"] == 0
    assert "radial_width_pixels" not in artery.attrs


def test_mask_detection_circular_area_is_enabled_by_default() -> None:
    segments = _segments()
    outputs = pack_mask_detection_blood_volume_rate_outputs(
        segments,
        segments,
        np.asarray([0, 2], dtype=np.int32),
        optic_disc_center=(20.0, 20.0),
        pixel_size_mm=0.1,
    )
    artery = outputs["Processing/BloodVolumeRate/Artery/maskDetection/value"]
    radius_scale = np.hypot(20.0, 20.0)
    diameter_mm = np.pi * radius_scale * np.asarray([0.25, 0.75]) * 0.1
    expected = 2.0 * np.pi / 4.0 * diameter_mm**2
    np.testing.assert_allclose(
        artery.data[:, 0, 0, :],
        np.broadcast_to(expected, (2, 2)),
        rtol=0.01,
    )
    assert artery.attrs["unit"] == "mm^3/s"
    assert artery.attrs["circular_area_scaling_applied"] == 1
