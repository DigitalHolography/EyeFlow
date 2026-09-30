"""Typed state produced by the waveform-velocity pipeline."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Literal

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisResult
from calculations.segment_profiles import SegmentProfileResult
from pipelines.retinal_velocity.models import RetinalVelocity

from .sources import WaveformVelocitySourceData

VesselName = Literal["artery", "vein"]


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


@dataclass(frozen=True)
class WaveformVelocity:
    """Canonical waveform and spatial state shared with downstream pipelines."""

    retinal_velocity: RetinalVelocity
    source_data: WaveformVelocitySourceData
    artery_segments: VelocitySegmentResult | None
    vein_segments: VelocitySegmentResult | None
    per_beat_result: PerBeatAnalysisResult | None
    attrs: dict[str, object]

    def segments(self, vessel: VesselName) -> VelocitySegmentResult | None:
        """Return the extracted spatial segments for one vessel class."""

        if vessel == "artery":
            return self.artery_segments
        if vessel == "vein":
            return self.vein_segments
        raise ValueError("vessel must be 'artery' or 'vein'.")

    def require_per_beat(self) -> PerBeatAnalysisResult:
        """Return per-beat waveforms or explain why they are unavailable."""

        if self.per_beat_result is None:
            raise RuntimeError(
                "Per-beat waveform velocity was not requested for this run."
            )
        return self.per_beat_result

    @property
    def cycle_boundary_indexes(self) -> np.ndarray:
        if self.per_beat_result is not None:
            return self.per_beat_result.cycle_boundary_indexes
        return self.retinal_velocity.cycle_boundary_indexes


__all__ = ["VelocitySegmentResult", "WaveformVelocity"]
