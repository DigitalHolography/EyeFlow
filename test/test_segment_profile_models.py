"""Tests for generic segment-profile data models."""

from dataclasses import fields, replace

import numpy as np
import pytest

from calculations.vessel_segments.measurement import (
    CompactSegmentMaps,
    SegmentMeasurements,
    SegmentMeasurementSettings,
)
from pipelines.velocity_analysis.models import VelocitySegmentResult


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("pixel_size_mm", 0),
    ],
)
def test_invalid_settings_are_rejected(field, value):
    with pytest.raises(ValueError):
        replace(SegmentMeasurementSettings(0.01), **{field: value})


def test_velocity_result_composes_generic_profiles_with_pipeline_products():
    base_fields = {item.name for item in fields(SegmentMeasurements)}
    velocity_fields = {item.name for item in fields(VelocitySegmentResult)}

    assert not issubclass(VelocitySegmentResult, SegmentMeasurements)
    assert base_fields == {
        "topology",
        "segment_signal",
        "transverse",
        "longitudinal",
        "mean_images",
        "sample_spacing_mm",
        "maps",
    }
    assert velocity_fields == {"profile", "transverse_fft", "filtered_signal"}


def test_compact_segment_maps_resolve_rows_and_expand_to_the_segment_grid():
    values = np.arange(8, dtype=np.float32).reshape(2, 1, 2, 2)
    maps = CompactSegmentMaps(
        values=values,
        indexes=np.asarray([[0, 1], [1, 0]], dtype=np.int32),
    )

    assert maps.row_for(0, 1) == 0
    assert maps.row_for(1, 1) is None
    dense = maps.to_dense((2, 2))
    np.testing.assert_array_equal(dense[0, 1], values[0])
    np.testing.assert_array_equal(dense[1, 0], values[1])
    assert np.all(np.isnan(dense[0, 0]))
