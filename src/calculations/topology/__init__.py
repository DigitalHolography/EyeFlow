"""General-purpose topology calculations for retinal maps."""

from .cache import (
    TOPOLOGY_CACHE_STATE,
    TopologyCacheKey,
    run_topology_cache,
    topology_cache_key,
    topology_source_id,
)
from .branch_identity import (
    BranchIdentityResult,
    BranchIdentityStages,
    label_vessel_branches,
)
from .geometry import (
    AnnulusGeometry,
    annulus_mask,
    image_half_diagonal,
    ring_masks,
    section_masks,
)
from .optic_disc import OpticDisc
from .profiles import (
    longitudinal_profiles,
    mean_profiles,
    profile_deviation_power,
    transverse_profiles,
)
from .profile_interpolation import interpolate_profiles_per_beat
from .mask_area import annulus_widths_pixels, segment_mask_areas_pixels
from .scale import retinal_pixel_size_mm
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
    prepare_segments,
    prepare_segment_chunks,
    prepare_topologies,
    prepare_topology,
    resolve_segment_rotations,
)

__all__ = [
    "BranchIdentityResult",
    "BranchIdentityStages",
    "PreparedSegments",
    "PreparedSegment",
    "PreparedSegmentChunk",
    "PreparedSegmentChunks",
    "PreparedTopology",
    "AnnulusGeometry",
    "OpticDisc",
    "SegmentTopology",
    "TOPOLOGY_CACHE_STATE",
    "TopologyCacheKey",
    "annulus_mask",
    "build_segment_topology",
    "competing_segment_masks",
    "determine_segment_rotations",
    "dilate_segment_masks",
    "extract_segments",
    "annulus_widths_pixels",
    "extract_segment",
    "resize_segment_topology_windows",
    "image_half_diagonal",
    "interpolate_segment_masks",
    "interpolate_segments",
    "interpolate_profiles_per_beat",
    "resample_rotate_segment",
    "label_vessel_branches",
    "longitudinal_profiles",
    "mean_profiles",
    "prepare_segments",
    "prepare_segment_chunks",
    "prepare_topologies",
    "prepare_topology",
    "resolve_segment_rotations",
    "profile_deviation_power",
    "retinal_pixel_size_mm",
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
