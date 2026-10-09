"""Typed state produced by the velocity-analysis pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisResult
from calculations.topology.segment_profiles import MaskedArrays, SegmentProfileResult
from pipelines.velocity.models import RetinalVelocity

from .sources import VelocityAnalysisSourceData

VesselName = Literal["artery", "vein"]


@dataclass(frozen=True, kw_only=True)
class VelocitySegmentResult:
    """A generic segment-profile result plus velocity-only products."""

    profile: SegmentProfileResult
    transverse_fft: MaskedArrays | None = None

    @classmethod
    def from_profile_result(
        cls,
        profiles: SegmentProfileResult,
        *,
        transverse_fft_profiles_unmasked: np.ndarray | None = None,
        transverse_fft_profiles_masked: np.ndarray | None = None,
    ) -> VelocitySegmentResult:
        if (transverse_fft_profiles_unmasked is None) != (
            transverse_fft_profiles_masked is None
        ):
            raise ValueError("masked and unmasked FFT profiles must be supplied together.")
        fft = (
            None
            if transverse_fft_profiles_unmasked is None
            else MaskedArrays(
                unmasked=transverse_fft_profiles_unmasked,
                masked=transverse_fft_profiles_masked,
            )
        )
        return cls(profile=profiles, transverse_fft=fft)

    def require_fft(self) -> MaskedArrays:
        if self.transverse_fft is None:
            raise RuntimeError("Velocity FFT profiles were not requested for this run.")
        return self.transverse_fft


@dataclass(frozen=True)
class VelocityAnalysis:
    """Canonical temporal and spatial state shared with downstream pipelines."""

    velocity: RetinalVelocity
    source_data: VelocityAnalysisSourceData
    artery_segments: VelocitySegmentResult | None
    vein_segments: VelocitySegmentResult | None
    per_beat_result: PerBeatAnalysisResult
    attrs: dict[str, object]

    def segments(self, vessel: VesselName) -> VelocitySegmentResult | None:
        """Return the extracted spatial segments for one vessel class."""

        if vessel == "artery":
            return self.artery_segments
        if vessel == "vein":
            return self.vein_segments
        raise ValueError("vessel must be 'artery' or 'vein'.")

    @property
    def cycle_boundary_indexes(self) -> np.ndarray:
        return self.per_beat_result.cycle_boundary_indexes


__all__ = ["VelocityAnalysis", "VelocitySegmentResult"]
