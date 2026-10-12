"""Compatibility imports for generic vessel-segment measurements."""

from calculations.vessel_segments.measurement import (
    CompactSegmentMaps,
    MaskedArrays,
    SegmentMeasurements,
    SegmentMeasurementSettings,
    analyze_segment_profiles,
)

SegmentProfileResult = SegmentMeasurements
SegmentProfileSettings = SegmentMeasurementSettings

__all__ = [
    "CompactSegmentMaps",
    "MaskedArrays",
    "SegmentProfileResult",
    "SegmentProfileSettings",
    "analyze_segment_profiles",
]
