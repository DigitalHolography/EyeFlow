"""Velocity-specific extensions to generic segment profile results."""

from __future__ import annotations

from dataclasses import dataclass, field, fields

import numpy as np

from calculations.segment_profiles import SegmentProfileResult
from pipelines.displacement_map.segments import DisplacementSegmentResult


@dataclass(frozen=True, kw_only=True)
class VelocitySegmentResult(SegmentProfileResult):
    """Generic segment profiles with optional velocity pipeline products."""

    displacements: dict[str, DisplacementSegmentResult] = field(default_factory=dict)
    transverse_fft_profiles_unmasked: np.ndarray | None = None
    transverse_fft_profiles_masked: np.ndarray | None = None

    @classmethod
    def from_profile_result(
        cls,
        profiles: SegmentProfileResult,
        *,
        transverse_fft_profiles_unmasked: np.ndarray | None = None,
        transverse_fft_profiles_masked: np.ndarray | None = None,
    ) -> VelocitySegmentResult:
        values = {
            item.name: getattr(profiles, item.name)
            for item in fields(SegmentProfileResult)
        }
        return cls(
            **values,
            transverse_fft_profiles_unmasked=transverse_fft_profiles_unmasked,
            transverse_fft_profiles_masked=transverse_fft_profiles_masked,
        )


__all__ = ["VelocitySegmentResult"]
