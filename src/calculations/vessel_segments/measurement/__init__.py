"""Measure signals, profiles, mean images, and optional maps from sampled segments."""

from .runner import analyze_segment_profiles
from .models import (
    CompactSegmentMaps,
    MaskedArrays,
    SegmentMeasurements,
    SegmentMeasurementSettings,
)

# Transitional names for callers of the former segment-profile API.
SegmentProfileResult = SegmentMeasurements
SegmentProfileSettings = SegmentMeasurementSettings

__all__ = [
    "CompactSegmentMaps",
    "MaskedArrays",
    "SegmentMeasurements",
    "SegmentMeasurementSettings",
    "SegmentProfileResult",
    "SegmentProfileSettings",
    "analyze_segment_profiles",
]
