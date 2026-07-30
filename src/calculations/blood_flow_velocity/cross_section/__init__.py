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
from .segment_geometry import SegmentRingSettings
from .segment_velocity_signals import segment_velocity_inputs, segment_velocity_results

__all__ = [
    "BranchIdentityResult",
    "BranchIdentityStages",
    "CrossSectionSignalResult",
    "CrossSectionProfileOutputs",
    "CrossSectionSignalSettings",
    "SegmentRingSettings",
    "generate_cross_section_signals",
    "label_vessel_branches",
    "segment_velocity_inputs",
    "segment_velocity_results",
]
