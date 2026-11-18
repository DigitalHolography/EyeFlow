"""Regression tests for reusable exact-annulus geometry and diagnostics."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

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


@pytest.mark.parametrize("count_based", [False, True])
def test_overlapping_measurements_use_thickness_and_add_histogram_samples(tmp_path, count_based):
    import h5py
    from input_output.schema import EyeFlowOutputPaths, PixelPitch
    from input_output.writers.h5_output import PipelineH5Output
    from pipelines.topology_core.outputs import pack_topology_outputs
    from pipelines.topology_core.annulus_visualization import annulus_outline_image
    from input_output.schema.serialization import serialize_spatial_image

    shape = (101, 101)
    vessel = np.zeros(shape, dtype=bool)
    vessel[48:53, 54:98] = True
    disc = OpticDisc(None, (50.0, 50.0), 8, 8)
    legacy_geometry = disc.annulus_geometry(shape, number_of_radii_in_fov=10)
    geometry = (disc.annulus_geometry(shape, number_of_radii_in_half_diagonal=32)
                if count_based else disc.annulus_geometry(shape, annulus_size=10.1))
    size = geometry.ring_width_frac * np.hypot(50, 50)
    legacy = prepare_sampling_plan(vessel, disc, legacy_geometry)
    prepared = prepare_sampling_plan(vessel, disc, geometry)
    diameters, areas, widths = mask_derived_lumen_geometry((prepared,), pixel_size_mm=0.01)
    legacy_diameters, _, _ = mask_derived_lumen_geometry((legacy,), pixel_size_mm=0.01)
    bounds = geometry.metadata(shape)["annulus_bounds_pixels"]
    np.testing.assert_allclose(widths, np.diff(bounds, axis=1).ravel(), rtol=1e-6)
    assert np.allclose(widths[:-2], size)  # Measurement thickness, never stride.
    if count_based:
        np.testing.assert_allclose(widths, size, rtol=1e-6)
    valid = areas[0] > 0
    np.testing.assert_allclose(diameters[0][valid], (areas[0] * 0.01 / widths)[valid], rtol=1e-6)
    assert np.isfinite(diameters[0]).sum() > np.isfinite(legacy_diameters[0]).sum()
    histogram = _lumen_diameter_distribution_figure(diameters[0].ravel() * 1000)
    assert sum(bar.get_height() for bar in histogram.axes[0].patches) == valid.sum()

    outputs = pack_topology_outputs(vessel, np.zeros(shape, dtype=bool), disc,
                                    {"artery": prepared}, PixelPitch(10e-6, 10e-6))
    schema = EyeFlowOutputPaths.active().segmentation
    with h5py.File(tmp_path / "geometry.h5", "w") as file:
        PipelineH5Output(file).write_many(outputs)
        assert file[schema.annuli.size][()] == pytest.approx(size)
        assert file[schema.annuli.stride][()] == pytest.approx(size / 2)
        assert file[schema.annuli.count][()] == geometry.ring_count
        np.testing.assert_allclose(file[schema.annuli.bounds_pixels][:], bounds)
        assert file[schema.artery.segment_mask_area].shape[-1] == geometry.ring_count
        np.testing.assert_allclose(file[schema.artery.lumen_diameter][:][valid],
                                   diameters[0][valid] / 0.01, rtol=1e-6)
        assert file[schema.artery.delta_radius].shape == (geometry.ring_count,)
        outline = file[schema.annuli.outline_intensity]
        assert outline.dtype == np.uint8
        assert outline.shape == shape[::-1]
        expected_outline = annulus_outline_image(shape, disc.center, geometry)
        np.testing.assert_array_equal(outline[:], serialize_spatial_image(expected_outline))
        np.testing.assert_array_equal(outline.attrs["annulus_gray_levels"], [255, 128])
        np.testing.assert_array_equal(outline.attrs["dimDesc"], ["x", "y"])
        assert outline.attrs["image_origin"] == "lower_left"
        assert outline.attrs["shared_boundary_policy"] == "earliest_annulus_wins"
