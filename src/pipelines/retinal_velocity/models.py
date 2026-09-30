"""Canonical retinal-velocity data shared by downstream pipelines."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import ClassVar, Literal

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from calculations.math.cycle_boundaries import normalize_cycle_boundaries

VesselName = Literal["artery", "vein"]


@dataclass(frozen=True)
class RetinalVelocity(Mapping[str, object]):
    """Base retinal velocity, frequency maps, and cardiac-cycle timing.

    The mapping interface is a temporary bridge for figure code that still
    consumes the former analysis dictionary. New code should use fields and
    methods directly.
    """

    artery_velocity_raw: np.ndarray
    vein_velocity_raw: np.ndarray
    artery_velocity_filtered: np.ndarray
    vein_velocity_filtered: np.ndarray
    velocity_map: object | None
    moment0_average: np.ndarray
    velocity_average: np.ndarray
    frms_average: np.ndarray
    frms_background_average: np.ndarray
    delta_frms_average: np.ndarray
    velocity_section_mask: np.ndarray
    artery_frms: np.ndarray
    vein_frms: np.ndarray
    artery_frms_background: np.ndarray
    vein_frms_background: np.ndarray
    vessel_frms_background: np.ndarray
    artery_delta_frms: np.ndarray
    vein_delta_frms: np.ndarray
    cardiac_cycle: CardiacCycleAnalysis
    cardiac_cycle_source: str
    dt_seconds: float
    provenance: Mapping[str, object] = field(default_factory=dict)
    index_base: int = 0

    _COMPATIBILITY_KEYS: ClassVar[tuple[str, ...]] = (
        "fRMS",
        "fRMS_bkg",
        "deltafRMS",
        "velocity_map",
        "retinal_vessel_velocity",
        "moment0_avg",
        "velocity_map_avg",
        "fRMS_avg",
        "fRMS_bkg_avg",
        "deltafRMS_avg",
        "velocity_section_mask",
        "velocity_section_geometry",
        "retinal_artery_velocity_signal",
        "retinal_vein_velocity_signal",
        "retinal_artery_velocity_signal_filtered",
        "retinal_vein_velocity_signal_filtered",
        "retinal_artery_velocity_signal_derivative",
        "retinal_vein_velocity_signal_derivative",
        "retinal_artery_velocity_signal_filtered_perbeat",
        "retinal_artery_fRMS_signal",
        "retinal_vein_fRMS_signal",
        "retinal_artery_fRMS_bkg_signal",
        "retinal_vein_fRMS_bkg_signal",
        "retinal_vessel_fRMS_bkg_signal",
        "retinal_artery_deltafRMS_signal",
        "retinal_vein_deltafRMS_signal",
        "beat_indices",
        "time_per_beat",
        "beat_detection_min_peak_distance",
        "beat_detection_min_peak_height",
        "cardiac_cycle_detection_source",
        "_cardiac_cycle_analysis",
        "provenance",
        "velocity_estimation_method",
        "velocity_quantity",
        "velocity_unit",
    )

    def continuous(
        self,
        vessel: VesselName,
        *,
        raw: bool = False,
    ) -> np.ndarray:
        """Return one whole-vessel velocity signal."""

        vessel_name = _validate_vessel(vessel)
        suffix = "raw" if raw else "filtered"
        return getattr(self, f"{vessel_name}_velocity_{suffix}")

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
        return self.velocity_map is not None

    def __getitem__(self, key: str) -> object:
        detection = self.cardiac_cycle.systole
        values = {
            "fRMS": lambda: None,
            "fRMS_bkg": lambda: None,
            "deltafRMS": lambda: None,
            "velocity_map": lambda: self.velocity_map,
            "retinal_vessel_velocity": lambda: self.velocity_map,
            "moment0_avg": lambda: self.moment0_average,
            "velocity_map_avg": lambda: self.velocity_average,
            "fRMS_avg": lambda: self.frms_average,
            "fRMS_bkg_avg": lambda: self.frms_background_average,
            "deltafRMS_avg": lambda: self.delta_frms_average,
            "velocity_section_mask": lambda: self.velocity_section_mask,
            "velocity_section_geometry": lambda: (
                "optic_disc_centered_frame_fraction"
            ),
            "retinal_artery_velocity_signal": lambda: self.artery_velocity_raw,
            "retinal_vein_velocity_signal": lambda: self.vein_velocity_raw,
            "retinal_artery_velocity_signal_filtered": lambda: (
                self.artery_velocity_filtered
            ),
            "retinal_vein_velocity_signal_filtered": lambda: (
                self.vein_velocity_filtered
            ),
            "retinal_artery_velocity_signal_derivative": lambda: self.derivative(
                "artery"
            ),
            "retinal_vein_velocity_signal_derivative": lambda: self.derivative(
                "vein"
            ),
            "retinal_artery_velocity_signal_filtered_perbeat": lambda: (
                _interpolate_cycles(self.per_beat("artery"), sample_count=128)
            ),
            "retinal_artery_fRMS_signal": lambda: self.artery_frms,
            "retinal_vein_fRMS_signal": lambda: self.vein_frms,
            "retinal_artery_fRMS_bkg_signal": lambda: self.artery_frms_background,
            "retinal_vein_fRMS_bkg_signal": lambda: self.vein_frms_background,
            "retinal_vessel_fRMS_bkg_signal": lambda: self.vessel_frms_background,
            "retinal_artery_deltafRMS_signal": lambda: self.artery_delta_frms,
            "retinal_vein_deltafRMS_signal": lambda: self.vein_delta_frms,
            "beat_indices": lambda: self.cycle_boundary_indexes,
            "time_per_beat": lambda: self.cycle_durations_seconds,
            "beat_detection_min_peak_distance": lambda: detection.min_peak_distance,
            "beat_detection_min_peak_height": lambda: detection.min_peak_height,
            "cardiac_cycle_detection_source": lambda: self.cardiac_cycle_source,
            "_cardiac_cycle_analysis": lambda: self.cardiac_cycle,
            "provenance": lambda: self.provenance,
            "velocity_estimation_method": lambda: self.provenance.get(
                "velocity_estimation_method",
                "doppler_moments",
            ),
            "velocity_quantity": lambda: self.provenance.get(
                "velocity_quantity",
                "physical_velocity",
            ),
            "velocity_unit": lambda: self.provenance.get(
                "velocity_unit",
                "mm/s",
            ),
        }
        if key in self.provenance:
            return self.provenance[key]
        try:
            return values[key]()
        except KeyError as exc:
            raise KeyError(key) from exc

    def __iter__(self) -> Iterator[str]:
        return iter(dict.fromkeys((*self._COMPATIBILITY_KEYS, *self.provenance)))

    def __len__(self) -> int:
        return len(tuple(iter(self)))


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


__all__ = ["RetinalVelocity", "VesselName"]
