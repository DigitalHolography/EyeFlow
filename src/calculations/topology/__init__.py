"""General-purpose topology calculations for retinal maps."""

from .branch_identity import (
    BranchIdentityResult,
    BranchIdentityStages,
    label_vessel_branches,
)
from .cache import (
    TOPOLOGY_CACHE_STATE,
    TopologyCacheKey,
    run_topology_cache,
    topology_cache_key,
    topology_source_id,
)
from .geometry import (
    AnnulusGeometry,
    annulus_mask,
    image_half_diagonal,
    ring_masks,
    section_masks,
)
from .mask_area import annulus_widths_pixels, segment_mask_areas_pixels
from .optic_disc import OpticDisc
from .profile_interpolation import interpolate_profiles_per_beat
from .profiles import (
    longitudinal_profiles,
    mean_profiles,
    profile_deviation_power,
    transverse_profiles,
)
from .quadrants import QUADRANT_NAMES, quadrant_membership
from .segments import (
    SegmentTopology,
    build_segment_topology,
    competing_segment_masks,
    extract_segment,
    extract_segments,
    resize_segment_topology_windows,
)
from .transforms import (
    determine_segment_rotations,
    dilate_segment_masks,
    interpolate_segment_masks,
    interpolate_segments,
    resample_rotate_segment,
    rotate_segment_masks,
    rotate_segments,
)
from .workflow import (
    PreparedSegment,
    PreparedSegmentChunk,
    PreparedSegmentChunks,
    PreparedSegments,
    PreparedTopology,
    prepare_segment_chunks,
    prepare_segments,
    prepare_topologies,
    prepare_topology,
    resolve_segment_rotations,
)

__all__ = [
    "TOPOLOGY_CACHE_STATE",
    "AnnulusGeometry",
    "BranchIdentityResult",
    "BranchIdentityStages",
    "OpticDisc",
    "PreparedSegment",
    "PreparedSegmentChunk",
    "PreparedSegmentChunks",
    "PreparedSegments",
    "PreparedTopology",
    "QUADRANT_NAMES",
    "SegmentTopology",
    "TopologyCacheKey",
    "annulus_mask",
    "annulus_widths_pixels",
    "build_segment_topology",
    "competing_segment_masks",
    "determine_segment_rotations",
    "dilate_segment_masks",
    "extract_segment",
    "extract_segments",
    "image_half_diagonal",
    "interpolate_profiles_per_beat",
    "interpolate_segment_masks",
    "interpolate_segments",
    "label_vessel_branches",
    "longitudinal_profiles",
    "mean_profiles",
    "prepare_segment_chunks",
    "prepare_segments",
    "prepare_topologies",
    "prepare_topology",
    "profile_deviation_power",
    "quadrant_membership",
    "resample_rotate_segment",
    "resize_segment_topology_windows",
    "resolve_segment_rotations",
    "ring_masks",
    "rotate_segment_masks",
    "rotate_segments",
    "run_topology_cache",
    "section_masks",
    "segment_mask_areas_pixels",
    "topology_cache_key",
    "topology_source_id",
    "transverse_profiles",
]
