"""Regression tests for reusable exact-annulus geometry and diagnostics."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from calculations.blood_volume_rate import (
    circular_lumen_flow,
    mask_derived_lumen_geometry,
)
from calculations.topology import AnnulusGeometry
from calculations.topology.mask_area import annulus_widths_pixels


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
        annulus_masks=sections,
        ring_settings=AnnulusGeometry(0.0, 0.5, 0.25, 2, 0.25),
    )
    return SimpleNamespace(
        projected_signal=np.full((2, 1, 1, 2), 2.0, dtype=np.float32),
        topology=topology,
    )


def test_annulus_width_uses_clipped_last_ring() -> None:
    widths = annulus_widths_pixels(
        (7, 9),
        AnnulusGeometry(0.0, 0.6, 0.5, 2, 0.25),
        2,
    )
    np.testing.assert_allclose(widths, [1.25, 0.5], rtol=1e-6)


def test_mask_derived_flow_keeps_equivalent_diameter_model() -> None:
    segments = _segments()
    diameters, areas, radial_widths = mask_derived_lumen_geometry(
        (segments.topology,),
        pixel_size_mm=0.1,
    )
    artery = circular_lumen_flow(segments.projected_signal, diameters[0])
    diameter_mm = areas[0][0] * 0.1 / radial_widths
    area = np.float32(np.pi / 4.0) * diameter_mm**2
    expected = np.float32(2.0) / np.sqrt(
        diameter_mm / np.float32(0.08)
    ) * area
    np.testing.assert_allclose(
        artery[:, 0, 0, :],
        np.broadcast_to(expected, (2, 2)),
        rtol=0.01,
    )
    assert areas[0].shape == (1, 2)
    assert radial_widths.shape == (2,)
