"""Regression tests for reusable exact-annulus geometry and diagnostics."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from calculations.blood_volume_rate import (
    mask_derived_lumen_geometry,
    masked_edges_flow,
)
from calculations.topology import AnnulusGeometry, OpticDisc, prepare_topology
from calculations.topology.mask_area import (
    annulus_widths_pixels,
    segment_mask_areas_pixels,
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
        annulus_masks=sections,
        ring_settings=AnnulusGeometry(0.0, 0.5, 0.25, 2, 0.25),
    )
    return SimpleNamespace(
        velocity=np.full((2, 1, 1, 2), 2.0, dtype=np.float32),
        topology=topology,
    )


def test_annulus_width_uses_clipped_last_ring() -> None:
    widths = annulus_widths_pixels(
        (7, 9),
        AnnulusGeometry(0.0, 0.6, 0.5, 2, 0.25),
        2,
    )
    np.testing.assert_allclose(widths, [1.25, 0.5], rtol=1e-6)


def test_canonical_prepared_topology_provides_native_mask_areas() -> None:
    vessel = np.zeros((41, 41), dtype=bool)
    vessel[18:23, 4:37] = True
    prepared = prepare_topology(
        vessel,
        OpticDisc(None, (20.0, 20.0), 6.0, 6.0),
        AnnulusGeometry(0.1, 0.7, 0.2, 3, 0.2),
        window_size_percentile_kept=1.0,
    )

    areas = segment_mask_areas_pixels(prepared.topology)
    diameters, geometry_areas, widths = mask_derived_lumen_geometry(
        (prepared,),
        pixel_size_mm=0.01,
    )

    assert areas.shape == (
        prepared.topology.branch_ids.size,
        prepared.topology.annulus_masks.shape[0],
    )
    np.testing.assert_array_equal(geometry_areas[0], areas)
    assert diameters[0].shape == areas.shape
    assert widths.shape == (prepared.topology.annulus_masks.shape[0],)


def test_mask_derived_flow_keeps_equivalent_diameter_model() -> None:
    segments = _segments()
    diameters, areas, radial_widths = mask_derived_lumen_geometry(
        (segments.topology,),
        pixel_size_mm=0.1,
    )
    artery = masked_edges_flow(segments.velocity, diameters[0])
    diameter_mm = areas[0][0] * 0.1 / radial_widths
    expected = 2.0 * np.pi / 4.0 * diameter_mm**2
    np.testing.assert_allclose(
        artery[:, 0, 0, :],
        np.broadcast_to(expected, (2, 2)),
        rtol=0.01,
    )
    assert areas[0].shape == (1, 2)
    assert radial_widths.shape == (2,)
