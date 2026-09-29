"""Tests for generic segment-profile data models."""

from dataclasses import fields, replace

import numpy as np
import pytest

from calculations.segment_profiles import SegmentProfileResult, SegmentProfileSettings
from pipelines.waveform_velocity_core.models import VelocitySegmentResult


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("pixel_size_mm", 0),
        ("pixel_size_mm", np.nan),
        ("working_memory_mb", -1),
        ("submask_size_percentile_kept", np.nan),
    ],
)
def test_invalid_settings_are_rejected(field, value):
    with pytest.raises(ValueError):
        replace(SegmentProfileSettings(0.01), **{field: value})


def test_velocity_result_only_adds_pipeline_specific_products():
    base_fields = {item.name for item in fields(SegmentProfileResult)}
    velocity_fields = {item.name for item in fields(VelocitySegmentResult)}

    assert issubclass(VelocitySegmentResult, SegmentProfileResult)
    assert velocity_fields - base_fields == {
        "displacements",
        "transverse_fft_profiles_masked",
        "transverse_fft_profiles_unmasked",
    }
