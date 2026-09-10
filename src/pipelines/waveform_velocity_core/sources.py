"""Source assembly for shared waveform velocity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from calculations.blood_flow_velocity import CrossSectionSignalSettings
from input_output.schema import DopplerViewSource, HolodopplerSource, HolodopplerTiming

from .constants import (
    CROSS_SECTION_SUBMASK_SIZE_PERCENTILE_KEPT,
    DEFAULT_PIXEL_SIZE_MM,
    REFERENCE_OPTIC_DISC_DIAMETER_MM,
    SPATIAL_INTERPOLATION_FACTOR,
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
    optic_disc_mask: np.ndarray | None
    optic_disc_center: np.ndarray | None
    optic_disc_width: np.ndarray | None
    optic_disc_height: np.ndarray | None
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
        vein_mask, vein_axes_swapped = _align_spatial_array(
            self.dv.retinal_vein_mask(),
            spatial_shape,
            "retinal_vein_mask",
        )
        spatial_axes_swapped = artery_axes_swapped or vein_axes_swapped
        labeled_vessels, labeled_axes_swapped = _align_optional_spatial_array(
            self.dv.retinal_labeled_vessels(),
            spatial_shape,
            "retinal_labeled_vessels",
        )
        spatial_axes_swapped = spatial_axes_swapped or labeled_axes_swapped
        optic_disc_mask, optic_disc_axes_swapped = _align_optional_spatial_array(
            self.dv.optic_disc_mask(),
            spatial_shape,
            "optic_disc_mask",
        )
        spatial_axes_swapped = spatial_axes_swapped or optic_disc_axes_swapped
        optic_disc_center = _align_spatial_pair(
            self.dv.optic_disc_center(),
            swapped=spatial_axes_swapped,
        )
        optic_disc_width, optic_disc_height = _align_spatial_sizes(
            self.dv.optic_disc_width(),
            self.dv.optic_disc_height(),
            swapped=spatial_axes_swapped,
        )
        return WaveformVelocitySourceData(
            moment0=moment0,
            moment2=moment2,
            retinal_artery_mask=np.asarray(artery_mask, dtype=bool),
            retinal_vein_mask=np.asarray(vein_mask, dtype=bool),
            retinal_labeled_vessels=labeled_vessels,
            optic_disc_mask=optic_disc_mask,
            optic_disc_center=optic_disc_center,
            optic_disc_width=optic_disc_width,
            optic_disc_height=optic_disc_height,
            timing=timing,
            local_background_dist=self.dv.local_background_dist(),
            cross_section_settings=self._cross_section_settings(
                optic_disc_width,
                optic_disc_height,
            ),
            provenance=_source_provenance(
                self.hd,
                self.dv,
                labeled_vessels,
                optic_disc_mask,
                optic_disc_center,
                spatial_axes_swapped=spatial_axes_swapped,
            ),
        )

    def _cross_section_settings(self, optic_disc_width, optic_disc_height):
        return CrossSectionSignalSettings(
            pixel_size_mm=self._pixel_size(optic_disc_width, optic_disc_height),
            submask_size_percentile_kept=(
                CROSS_SECTION_SUBMASK_SIZE_PERCENTILE_KEPT
            ),
        )

    def _pixel_size(self, optic_disc_width, optic_disc_height) -> float:
        diameter = _mean_pair(optic_disc_width, optic_disc_height)
        if diameter is not None:
            return REFERENCE_OPTIC_DISC_DIAMETER_MM / diameter
        return DEFAULT_PIXEL_SIZE_MM / (2.0**SPATIAL_INTERPOLATION_FACTOR)


def load_waveform_velocity_source_data(ctx: PipelineContext) -> WaveformVelocitySourceData:
    return WaveformVelocitySources.from_context(ctx).load()


def _load_moment_pair(hd: HolodopplerSource):
    return hd.moment0_dataset(), hd.moment2_dataset()


def _source_provenance(
    hd,
    dv,
    labeled_vessels,
    optic_disc_mask,
    optic_disc_center,
    *,
    spatial_axes_swapped: bool,
) -> dict[str, object]:
    return {
        "hd_source_file": str(hd.filename or ""),
        "dv_source_file": str(dv.filename or ""),
        "has_retinal_labeled_vessels": labeled_vessels is not None,
        "has_optic_disc_mask": optic_disc_mask is not None,
        "has_optic_disc_center": optic_disc_center is not None,
        "dv_spatial_axes_swapped_to_match_hd": spatial_axes_swapped,
        "beat_index_base": BEAT_INDEX_BASE,
        "moment_axes": list(MOMENT_AXES),
        "mask_axes": list(MASK_AXES),
    }


def _align_optional_spatial_array(
    array,
    spatial_shape: tuple[int, int],
    name: str,
):
    if array is None:
        return None, False
    return _align_spatial_array(array, spatial_shape, name)


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


def _align_spatial_pair(value, *, swapped: bool):
    if value is None or not swapped:
        return value
    array = np.asarray(value).copy()
    flat = array.reshape(-1)
    if flat.size >= 2:
        flat[0], flat[1] = flat[1].copy(), flat[0].copy()
    return array


def _align_spatial_sizes(width, height, *, swapped: bool):
    if swapped:
        return height, width
    return width, height


def _mean_pair(first, second) -> float | None:
    if first is None or second is None:
        return None
    values = np.asarray([first, second], dtype=np.float32).reshape(-1)
    if values.size < 2 or not np.all(np.isfinite(values[:2])):
        return None
    return float(np.mean(values[:2], dtype=np.float32))
