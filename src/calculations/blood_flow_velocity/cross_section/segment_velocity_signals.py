"""Segment velocity arrays from CrossSection/generateCrossSectionSignals.m."""

from __future__ import annotations

import numpy as np

from .generate_cross_section_signals import (
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    generate_cross_section_signals,
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
    velocity_array = np.asarray(velocity, dtype=np.float32)
    settings = _cross_section_settings(cross_section_settings)
    return (
        _segment_velocity(velocity_array, artery_mask, optic_disc_center, ring_settings, settings),
        _segment_velocity(velocity_array, vein_mask, optic_disc_center, ring_settings, settings),
    )


def _segment_velocity(
    velocity: np.ndarray,
    vessel_mask,
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
) -> CrossSectionSignalResult:
    return generate_cross_section_signals(
        velocity,
        np.asarray(vessel_mask, dtype=bool),
        optic_disc_center,
        ring_settings,
        cross_section_settings,
    )


def _cross_section_settings(value: CrossSectionSignalSettings | None):
    if value is not None:
        return value
    return CrossSectionSignalSettings(
        scale_factor_width=3.0,
        hydrodynamic_diameters=True,
        velocity_profile_threshold=0.5,
        rotate_from_mask=False,
        pixel_size_mm=0.0191,
    )
