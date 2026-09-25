"""Source assembly for shared waveform velocity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from calculations.blood_flow_velocity import CrossSectionSignalSettings
from calculations.topology import OpticDisc, retinal_pixel_size_mm
from input_output.schema import DopplerViewSource, HolodopplerSource, HolodopplerTiming

from .constants import (
    CROSS_SECTION_SUBMASK_SIZE_PERCENTILE_KEPT,
)

if TYPE_CHECKING:
    from pipeline_engine import PipelineContext


MOMENT_AXES = ("frame", "y", "x")
MASK_AXES = ("y", "x")
BEAT_INDEX_BASE = 0


@dataclass(frozen=True)
class WaveformVelocitySourceData:
    """Resolved source data with explicit axis contracts for waveform metrics."""

    moment0: object
    moment2: object
    retinal_artery_mask: np.ndarray
    retinal_vein_mask: np.ndarray
    retinal_labeled_vessels: np.ndarray | None
    optic_disc: OpticDisc
    timing: HolodopplerTiming
    local_background_dist: int
    cross_section_settings: CrossSectionSignalSettings
    provenance: dict[str, object]


@dataclass(frozen=True)
class WaveformVelocitySources:
    """Typed input adapters needed by the waveform velocity core."""

    hd: HolodopplerSource
    dv: DopplerViewSource

    @classmethod
    def from_context(cls, ctx: PipelineContext) -> "WaveformVelocitySources":
        ctx.require_inputs("hd", "dv")
        return cls(
            hd=ctx.inputs.hd.as_holodoppler(),
            dv=ctx.inputs.dv.as_dopplerview(),
        )

    def load(self) -> WaveformVelocitySourceData:
        moment0, moment2 = _load_moment_pair(self.hd)
        spatial_shape = moment0.shape[-2:]
        timing = self.hd.timing()
        artery_mask, artery_axes_swapped = _align_spatial_array(
            self.dv.retinal_artery_mask(),
            spatial_shape,
            "retinal_artery_mask",
        )
        vein_mask = _apply_spatial_alignment(
            self.dv.retinal_vein_mask(),
            spatial_shape,
            "retinal_vein_mask",
            swapped=artery_axes_swapped,
        )
        spatial_axes_swapped = artery_axes_swapped
        labeled_vessels = _apply_optional_spatial_alignment(
            self.dv.retinal_labeled_vessels(),
            spatial_shape,
            "retinal_labeled_vessels",
            swapped=spatial_axes_swapped,
        )
        optic_disc = self.dv.optic_disc()
        if spatial_axes_swapped:
            optic_disc = optic_disc.transposed()
        if optic_disc.mask is not None and optic_disc.mask.shape != tuple(spatial_shape):
            raise ValueError(
                "optic-disc mask must share the DopplerView vessel-mask frame; "
                f"expected {tuple(spatial_shape)}, got {optic_disc.mask.shape}."
            )
        return WaveformVelocitySourceData(
            moment0=moment0,
            moment2=moment2,
            retinal_artery_mask=np.asarray(artery_mask, dtype=bool),
            retinal_vein_mask=np.asarray(vein_mask, dtype=bool),
            retinal_labeled_vessels=labeled_vessels,
            optic_disc=optic_disc,
            timing=timing,
            local_background_dist=self.dv.local_background_dist(),
            cross_section_settings=self._cross_section_settings(optic_disc),
            provenance=_source_provenance(
                self.hd,
                self.dv,
                labeled_vessels,
                optic_disc,
                spatial_axes_swapped=spatial_axes_swapped,
            ),
        )

    def _cross_section_settings(self, optic_disc: OpticDisc):
        return CrossSectionSignalSettings(
            pixel_size_mm=self._pixel_size(optic_disc),
            submask_size_percentile_kept=(
                CROSS_SECTION_SUBMASK_SIZE_PERCENTILE_KEPT
            ),
        )

    def _pixel_size(self, optic_disc: OpticDisc) -> float:
        return retinal_pixel_size_mm(optic_disc)


def load_waveform_velocity_source_data(ctx: PipelineContext) -> WaveformVelocitySourceData:
    return WaveformVelocitySources.from_context(ctx).load()


def _load_moment_pair(hd: HolodopplerSource):
    return hd.moment0_dataset(), hd.moment2_dataset()


def _source_provenance(
    hd,
    dv,
    labeled_vessels,
    optic_disc: OpticDisc,
    *,
    spatial_axes_swapped: bool,
) -> dict[str, object]:
    return {
        "hd_source_file": str(hd.filename or ""),
        "dv_source_file": str(dv.filename or ""),
        "has_retinal_labeled_vessels": labeled_vessels is not None,
        "has_optic_disc_mask": optic_disc.mask is not None,
        "has_optic_disc_center": True,
        "dv_spatial_axes_swapped_to_match_hd": spatial_axes_swapped,
        "beat_index_base": BEAT_INDEX_BASE,
        "moment_axes": list(MOMENT_AXES),
        "mask_axes": list(MASK_AXES),
    }


def _apply_optional_spatial_alignment(
    array,
    spatial_shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
):
    if array is None:
        return None
    return _apply_spatial_alignment(array, spatial_shape, name, swapped=swapped)


def _align_spatial_array(
    array,
    spatial_shape: tuple[int, int],
    name: str,
):
    value = np.asarray(array)
    if value.ndim < 2:
        raise ValueError(
            f"{name} must include two spatial axes, got shape {value.shape}."
        )
    expected = tuple(int(size) for size in spatial_shape)
    actual = tuple(int(size) for size in value.shape[-2:])
    if actual == expected:
        return value, False
    if actual == expected[::-1]:
        return np.swapaxes(value, -1, -2), True
    raise ValueError(
        f"{name} spatial shape {actual} does not match HD spatial shape "
        f"{expected}, even after transposition."
    )


def _apply_spatial_alignment(
    array,
    spatial_shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
):
    value = np.asarray(array)
    if value.ndim < 2:
        raise ValueError(
            f"{name} must include two spatial axes, got shape {value.shape}."
        )
    aligned = np.swapaxes(value, -1, -2) if swapped else value
    expected = tuple(int(size) for size in spatial_shape)
    if tuple(int(size) for size in aligned.shape[-2:]) != expected:
        raise ValueError(
            f"{name} does not share the selected DopplerView spatial orientation; "
            f"aligned shape {aligned.shape[-2:]} does not match {expected}."
        )
    return aligned
