"""Segment velocity arrays from CrossSection/generateCrossSectionSignals.m."""

from __future__ import annotations

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from utils.logger import Logger

from .generate_cross_section_signals import (
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    _fixed_substack_side_pixels,
    _generate_cross_section_signals_from_geometry,
    _prepare_cross_section_geometry,
)
from .segment_geometry import SegmentRingSettings


def segment_velocity_inputs(
    velocity,
    artery_mask,
    vein_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    artery, vein = segment_velocity_results(
        velocity,
        artery_mask,
        vein_mask,
        optic_disc_center,
        ring_settings,
        cross_section_settings,
    )
    return artery.velocity, vein.velocity


def segment_velocity_results(
    velocity,
    artery_mask,
    vein_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings | None = None,
) -> tuple[CrossSectionSignalResult, CrossSectionSignalResult]:
    settings = _cross_section_settings(cross_section_settings)
    backend = optional_cupy_backend()
    Logger.log(
        "Cross-section compute backend: "
        + ("CuPy/CUDA" if backend is not None else "CPU with parallel segments")
        + "."
    )
    artery_geometry = _prepare_cross_section_geometry(
        artery_mask,
        optic_disc_center,
        ring_settings,
    )
    vein_geometry = _prepare_cross_section_geometry(
        vein_mask,
        optic_disc_center,
        ring_settings,
    )
    substack_side_pixels = _fixed_substack_side_pixels(
        (artery_geometry, vein_geometry),
        settings.submask_size_percentile_kept,
    )
    return (
        _generate_cross_section_signals_from_geometry(
            velocity,
            artery_geometry,
            optic_disc_center,
            ring_settings,
            settings,
            substack_side_pixels,
        ),
        _generate_cross_section_signals_from_geometry(
            velocity,
            vein_geometry,
            optic_disc_center,
            ring_settings,
            settings,
            substack_side_pixels,
        ),
    )


def _cross_section_settings(value: CrossSectionSignalSettings | None):
    if value is not None:
        return value
    return CrossSectionSignalSettings(
        hydrodynamic_diameters=True,
        velocity_profile_threshold=0.5,
        rotate_from_mask=False,
        pixel_size_mm=0.0191,
    )
