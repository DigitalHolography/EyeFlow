"""Velocity-specific extensions to generic segment profile results."""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from calculations.segment_profiles import SegmentProfileResult


@dataclass(frozen=True, kw_only=True)
class VelocitySegmentResult(SegmentProfileResult):
    """Generic segment profiles with optional velocity pipeline products."""

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
