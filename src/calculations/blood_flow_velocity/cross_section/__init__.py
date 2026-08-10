"""Cross-section geometry and velocity-signal calculations."""

from .branch_identity import (
    BranchIdentityResult,
    BranchIdentityStages,
    label_vessel_branches,
)
from .generate_cross_section_signals import (
    CrossSectionProfileOutputs,
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    generate_cross_section_signals,
)
from .reusable_cross_section_signals import (
    CrossSectionProjectionPlan,
    MultiCubeCrossSectionSignalResult,
    fit_cross_section_plan,
    generate_cross_section_signals_for_cubes,
    project_cross_section_cube,
    project_cross_section_cubes,
)
from .segment_geometry import SegmentRingSettings
from .segment_velocity_signals import segment_velocity_inputs, segment_velocity_results

__all__ = [
    "BranchIdentityResult",
    "BranchIdentityStages",
    "CrossSectionProfileOutputs",
    "CrossSectionProjectionPlan",
    "CrossSectionSignalResult",
    "CrossSectionSignalSettings",
    "MultiCubeCrossSectionSignalResult",
    "SegmentRingSettings",
    "fit_cross_section_plan",
    "generate_cross_section_signals",
    "generate_cross_section_signals_for_cubes",
    "label_vessel_branches",
    "project_cross_section_cube",
    "project_cross_section_cubes",
    "segment_velocity_inputs",
    "segment_velocity_results",
]
