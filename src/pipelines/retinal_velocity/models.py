"""Typed retinal-velocity calculation and pipeline data."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from calculations.math.cycle_boundaries import normalize_cycle_boundaries

VesselName = Literal["artery", "vein"]


@dataclass(frozen=True)
class RetinalVelocityMaps:
    """Spatial maps shared throughout retinal-velocity processing."""

    velocity: object | None
    moment0_average: np.ndarray
    velocity_average: np.ndarray
    frms_average: np.ndarray
    frms_background_average: np.ndarray
    delta_frms_average: np.ndarray
    section_mask: np.ndarray


@dataclass(frozen=True)
class VesselVelocitySignals:
    """Unfiltered signals calculated for one retinal vessel class."""

    velocity: np.ndarray
    frms: np.ndarray
    frms_background: np.ndarray
    delta_frms: np.ndarray


@dataclass(frozen=True)
class RetinalVelocityData:
    """Neutral components available before cardiac-cycle processing."""

    maps: RetinalVelocityMaps
    artery: VesselVelocitySignals
    vein: VesselVelocitySignals
    vessel_frms_background: np.ndarray


@dataclass(frozen=True)
class VesselVelocity:
    """Raw vessel signals paired with their filtered velocity."""

    signals: VesselVelocitySignals
    velocity_filtered: np.ndarray

    def continuous(self, *, raw: bool = False) -> np.ndarray:
        return self.signals.velocity if raw else self.velocity_filtered

    @property
    def frms(self) -> np.ndarray:
        return self.signals.frms

    @property
    def frms_background(self) -> np.ndarray:
        return self.signals.frms_background

    @property
    def delta_frms(self) -> np.ndarray:
        return self.signals.delta_frms


@dataclass(frozen=True)
class RetinalVelocity:
    """Base retinal velocity, frequency maps, and cardiac-cycle timing."""

    maps: RetinalVelocityMaps
    artery: VesselVelocity
    vein: VesselVelocity
    vessel_frms_background: np.ndarray
    cardiac_cycle: CardiacCycleAnalysis
    cardiac_cycle_source: str
    dt_seconds: float
    index_base: int = 0

    def continuous(
        self,
        vessel: VesselName,
        *,
        raw: bool = False,
    ) -> np.ndarray:
        """Return one whole-vessel velocity signal."""

        return self.vessel(vessel).continuous(raw=raw)

    def vessel(self, vessel: VesselName) -> VesselVelocity:
        """Return the composed signals for one vessel class."""

        return getattr(self, _validate_vessel(vessel))

    def per_beat(
        self,
        vessel: VesselName,
        *,
        raw: bool = False,
    ) -> tuple[np.ndarray, ...]:
        """Slice a continuous signal between successive cardiac boundaries.

        Cycle lengths are intentionally retained. Waveform pipelines may
        interpolate or Fourier-normalize these arrays when a dense matrix is
        required.
        """

        signal = self.continuous(vessel, raw=raw)
        boundaries = normalize_cycle_boundaries(
            self.cycle_boundary_indexes,
            signal.size,
            index_base=0,
        )
        return tuple(
            signal[int(start) : int(stop) + 1]
            for start, stop in zip(boundaries[:-1], boundaries[1:])
        )

    def derivative(
        self,
        vessel: VesselName,
        *,
        raw: bool = False,
    ) -> np.ndarray:
        """Return the temporal derivative of a whole-vessel signal."""

        return np.gradient(
            self.continuous(vessel, raw=raw),
            np.float32(self.dt_seconds),
        ).astype(np.float32)

    def per_beat_matrix(
        self,
        vessel: VesselName,
        *,
        raw: bool = False,
        sample_count: int = 128,
    ) -> np.ndarray:
        """Return cardiac cycles interpolated to a common sample count."""

        if sample_count < 1:
            raise ValueError("sample_count must be positive.")
        return _interpolate_cycles(
            self.per_beat(vessel, raw=raw),
            sample_count=sample_count,
        )

    @property
    def cycle_boundary_indexes(self) -> np.ndarray:
        indexes = np.asarray(
            self.cardiac_cycle.systole.systole_indexes,
            dtype=np.int32,
        )
        return (indexes - int(self.index_base)).astype(np.int32, copy=False)

    @property
    def cycle_durations_seconds(self) -> np.ndarray:
        return (
            np.diff(self.cycle_boundary_indexes).astype(np.float32)
            * np.float32(self.dt_seconds)
        ).astype(np.float32, copy=False)

    @property
    def cycle_count(self) -> int:
        return max(0, int(self.cycle_boundary_indexes.size) - 1)

    @property
    def has_velocity_map(self) -> bool:
        return self.maps.velocity is not None


def _validate_vessel(vessel: str) -> VesselName:
    if vessel not in {"artery", "vein"}:
        raise ValueError("vessel must be 'artery' or 'vein'.")
    return vessel


def _interpolate_cycles(
    cycles: tuple[np.ndarray, ...],
    *,
    sample_count: int,
) -> np.ndarray:
    output = np.empty((len(cycles), sample_count), dtype=np.float32)
    target = np.linspace(0.0, 1.0, sample_count, dtype=np.float32)
    for index, cycle in enumerate(cycles):
        source = np.linspace(0.0, 1.0, cycle.size, dtype=np.float32)
        output[index] = np.interp(target, source, cycle).astype(np.float32)
    return output


__all__ = [
    "RetinalVelocity",
    "RetinalVelocityData",
    "RetinalVelocityMaps",
    "VesselName",
    "VesselVelocity",
    "VesselVelocitySignals",
]
