"""Spatial plan preparation, patch extraction, and bounded segment sampling."""

from .models import (
    SampledSegment,
    SampledSegmentChunk,
    SampledSegmentChunks,
    SampledSegments,
    SegmentSamplingPlan,
)
from .preparation import (
    prepare_sampling_plan,
    prepare_sampling_plans,
    resolve_segment_rotations,
)
from .extraction import extract_segment, extract_segments
from .streaming import prepare_segment_chunks, prepare_segments

# Transitional names available from the sampling layer only.
PreparedTopology = SegmentSamplingPlan
PreparedSegment = SampledSegment
PreparedSegmentChunk = SampledSegmentChunk
PreparedSegments = SampledSegments
PreparedSegmentChunks = SampledSegmentChunks
prepare_topology = prepare_sampling_plan
prepare_topologies = prepare_sampling_plans

__all__ = [
    "SegmentSamplingPlan",
    "SampledSegment",
    "SampledSegmentChunk",
    "SampledSegmentChunks",
    "SampledSegments",
    "extract_segment",
    "extract_segments",
    "prepare_sampling_plan",
    "prepare_sampling_plans",
    "prepare_segment_chunks",
    "prepare_segments",
    "resolve_segment_rotations",
    "PreparedTopology",
    "PreparedSegment",
    "PreparedSegmentChunk",
    "PreparedSegments",
    "PreparedSegmentChunks",
    "prepare_topology",
    "prepare_topologies",
]
