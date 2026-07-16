"""Source assembly for the waveform-shape metrics pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from calculations.blood_flow_velocity import CrossSectionSignalSettings
from input_output.schema import DopplerViewSource, HolodopplerSource, HolodopplerTiming

if TYPE_CHECKING:
    from pipeline_engine import PipelineContext


MOMENT_AXES = ("frame", "y", "x")
MASK_AXES = ("y", "x")
DOPPLERVIEW_BEAT_INDEX_BASE = 0


@dataclass(frozen=True)
class WaveformShapeSourceData:
    """Resolved source data with explicit axis contracts for waveform metrics."""

    moment0: np.ndarray
    moment2: np.ndarray
    retinal_artery_mask: np.ndarray
    retinal_vein_mask: np.ndarray
    retinal_labeled_vessels: np.ndarray | None
    optic_disc_mask: np.ndarray | None
    optic_disc_center: np.ndarray | None
    optic_disc_width: np.ndarray | None
    optic_disc_height: np.ndarray | None
    timing: HolodopplerTiming
    local_background_dist: int
    peripapillary_inner_radius: np.float32
    peripapillary_ring_count: np.int32
    peripapillary_ring_width: np.float32
    segment_length: np.float32
    cross_section_settings: CrossSectionSignalSettings
    dopplerview_analysis: dict[str, object] | None
    provenance: dict[str, object]

    def dopplerview_cache(self) -> dict[str, object]:
        return {
            "moment0": self.moment0,
            "moment2": self.moment2,
            "retinal_artery_mask": self.retinal_artery_mask,
            "retinal_vein_mask": self.retinal_vein_mask,
        }


@dataclass(frozen=True)
class WaveformShapeSources:
    """Typed input adapters needed by the waveform-shape pipeline."""

    hd: HolodopplerSource
    dv: DopplerViewSource

    @classmethod
    def from_context(cls, ctx: PipelineContext) -> "WaveformShapeSources":
        ctx.require_inputs("hd", "dv")
        return cls(
            hd=ctx.inputs.hd.as_holodoppler(),
            dv=ctx.inputs.dv.as_dopplerview(),
        )

    def load(self) -> WaveformShapeSourceData:
        moment0 = self.hd.moment0()
        moment2 = self.hd.moment2()
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
        dopplerview_analysis = _align_analysis_spatial_arrays(
            self.dv.analysis(default=None),
            spatial_shape,
        )
        return WaveformShapeSourceData(
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
            peripapillary_inner_radius=self._inner_radius(),
            peripapillary_ring_count=self._circle_count(),
            peripapillary_ring_width=self._ring_width(),
            segment_length=self._segment_length(),
            cross_section_settings=self._cross_section_settings(
                optic_disc_width,
                optic_disc_height,
            ),
            dopplerview_analysis=dopplerview_analysis,
            provenance=_source_provenance(
                self.hd,
                self.dv,
                labeled_vessels,
                optic_disc_mask,
                optic_disc_center,
                spatial_axes_swapped=spatial_axes_swapped,
                dopplerview_analysis_available=dopplerview_analysis is not None,
            ),
        )

    def _inner_radius(self) -> np.float32:
        return np.float32(
            _config_float(
                self.hd,
                "SizeOfField",
                "SmallRadiusRatio",
                float(self.dv.peripapillary_inner_radius()),
            )
        )

    def _circle_count(self) -> np.int32:
        return np.int32(
            _config_int(
                self.hd,
                "generateCrossSectionSignals",
                "NumberOfCircles",
                int(self.dv.peripapillary_ring_count()),
            )
        )

    def _ring_width(self) -> np.float32:
        return np.float32(self.dv.peripapillary_ring_width())

    def _segment_length(self) -> np.float32:
        return np.float32(
            _config_float(
                self.hd,
                "generateCrossSectionSignals",
                "SegmentsLength",
                -1.0,
            )
        )

    def _cross_section_settings(self, optic_disc_width, optic_disc_height):
        return CrossSectionSignalSettings(
            scale_factor_width=_config_float(
                self.hd,
                "generateCrossSectionSignals",
                "ScaleFactorWidth",
                3.0,
            ),
            hydrodynamic_diameters=_config_bool(
                self.hd,
                "generateCrossSectionSignals",
                "HydrodynamicDiameters",
                True,
            ),
            velocity_profile_threshold=_config_float(
                self.hd,
                "generateCrossSectionSignals",
                "velocityProfileThreshold",
                0.5,
            ),
            rotate_from_mask=_config_bool(
                self.hd,
                "generateCrossSectionSignals",
                "RotateFromMask",
                False,
            ),
            pixel_size_mm=self._pixel_size(optic_disc_width, optic_disc_height),
        )

    def _pixel_size(self, optic_disc_width, optic_disc_height) -> float:
        diameter = _mean_pair(optic_disc_width, optic_disc_height)
        if diameter is not None:
            ref = _config_float(self.hd, "generateCrossSectionSignals", "RefPapillaSize", 1.91)
            return ref / diameter
        default = _config_float(self.hd, "generateCrossSectionSignals", "DefaultPixelSize", 0.0191)
        factor = _config_float(self.hd, "Preprocess", "InterpolationFactor", 1.0)
        return default / (2.0 ** factor)


def load_waveform_shape_source_data(ctx: PipelineContext) -> WaveformShapeSourceData:
    return WaveformShapeSources.from_context(ctx).load()


def _source_provenance(
    hd,
    dv,
    labeled_vessels,
    optic_disc_mask,
    optic_disc_center,
    *,
    spatial_axes_swapped: bool,
    dopplerview_analysis_available: bool,
) -> dict[str, object]:
    return {
        "hd_source_file": str(hd.filename or ""),
        "dv_source_file": str(dv.filename or ""),
        "has_retinal_labeled_vessels": labeled_vessels is not None,
        "has_optic_disc_mask": optic_disc_mask is not None,
        "has_optic_disc_center": optic_disc_center is not None,
        "dopplerview_analysis_available": dopplerview_analysis_available,
        "dv_spatial_axes_swapped_to_match_hd": spatial_axes_swapped,
        "beat_index_base": DOPPLERVIEW_BEAT_INDEX_BASE,
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


def _align_analysis_spatial_arrays(
    analysis: dict[str, object] | None,
    spatial_shape: tuple[int, int],
) -> dict[str, object] | None:
    if analysis is None:
        return None
    aligned = dict(analysis)
    for key in (
        "retinal_vessel_velocity",
        "velocity_map_avg",
        "fRMS_avg",
        "fRMS_bkg_avg",
    ):
        aligned[key], _ = _align_spatial_array(aligned[key], spatial_shape, key)
    return aligned


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


def _config_float(source, section: str, key: str, default: float) -> float:
    value = source.config_value(section, key, default)
    array = np.asarray(value).reshape(-1)
    return float(default) if array.size == 0 else float(array[0])


def _config_int(source, section: str, key: str, default: int) -> int:
    return int(round(_config_float(source, section, key, float(default))))


def _config_bool(source, section: str, key: str, default: bool) -> bool:
    value = source.config_value(section, key, default)
    if isinstance(value, str):
        return value.lower() in {"1", "true", "yes"}
    return bool(np.asarray(value).reshape(-1)[0])


def _mean_pair(first, second) -> float | None:
    if first is None or second is None:
        return None
    values = np.asarray([first, second], dtype=np.float32).reshape(-1)
    if values.size < 2 or not np.all(np.isfinite(values[:2])):
        return None
    return float(np.mean(values[:2], dtype=np.float32))
