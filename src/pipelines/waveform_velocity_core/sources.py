"""Source assembly for shared waveform velocity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from calculations.segment_profiles import SegmentProfileSettings
from input_output.schema import (
    DopplerViewSource,
    HD_BAND_HF_PATH,
    HD_BAND_LF_PATH,
    HolodopplerSource,
    PixelPitch,
    RetinalSourceData,
)
from pipelines.vessel_inputs import load_retinal_source_data
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    physical_velocity_provenance,
)

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

    source: RetinalSourceData
    profile_settings: SegmentProfileSettings
    provenance: dict[str, object]


@dataclass(frozen=True)
class WaveformVelocitySources:
    """Typed input adapters needed by the waveform velocity core."""

    hd: HolodopplerSource
    dv: DopplerViewSource
    velocity_estimation_method: str = "doppler_moments"
    band_ratio_frequency_scale_hz: float = (
        DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    )

    @classmethod
    def from_context(cls, ctx: PipelineContext) -> WaveformVelocitySources:
        ctx.require_inputs("hd", "dv")
        return cls(
            hd=ctx.inputs.hd.as_holodoppler(),
            dv=ctx.inputs.dv.as_dopplerview(),
            velocity_estimation_method=ctx.velocity_estimation_method,
            band_ratio_frequency_scale_hz=ctx.band_ratio_frequency_scale_hz,
        )

    def load(self) -> WaveformVelocitySourceData:
        source = load_retinal_source_data(
            self.hd,
            self.dv,
            velocity_estimation_method=self.velocity_estimation_method,
            band_ratio_frequency_scale_hz=(
                self.band_ratio_frequency_scale_hz
            ),
        )
        pixel_pitch = source.holodoppler.pixel_pitch
        return WaveformVelocitySourceData(
            source=source,
            profile_settings=self._profile_settings(pixel_pitch),
            provenance=_source_provenance(self.hd, self.dv, source),
        )

    def _profile_settings(self, pixel_pitch: PixelPitch) -> SegmentProfileSettings:
        return SegmentProfileSettings(
            pixel_size_mm=self._pixel_size(pixel_pitch),
            submask_size_percentile_kept=(CROSS_SECTION_SUBMASK_SIZE_PERCENTILE_KEPT),
        )

    def _pixel_size(self, pixel_pitch: PixelPitch) -> float:
        return pixel_pitch.isotropic_mm


def load_waveform_velocity_source_data(ctx: PipelineContext) -> WaveformVelocitySourceData:
    return WaveformVelocitySources.from_context(ctx).load()


def _load_moment_pair(hd: HolodopplerSource):
    return hd.moment0_dataset(), hd.moment2_dataset()


def _source_provenance(
    hd,
    dv,
    source: RetinalSourceData,
) -> dict[str, object]:
    segmentation = source.segmentation
    optic_disc = segmentation.optic_disc
    return {
        "hd_source_file": str(hd.filename or ""),
        "dv_source_file": str(dv.filename or ""),
        "has_retinal_labeled_vessels": segmentation.vessels.labeled is not None,
        "has_optic_disc_mask": (optic_disc.mask is not None and not optic_disc.is_fallback),
        "has_optic_disc_center": True,
        "optic_disc_geometry_fallback": optic_disc.is_fallback,
        "vein_analysis_enabled": not optic_disc.is_fallback,
        "pixel_pitch_xy_m": list(source.holodoppler.pixel_pitch.xy_m),
        "dv_spatial_axes_swapped_to_match_hd": (
            source.doppler_view.spatial_axes_swapped_to_match_hd
        ),
        "beat_index_base": BEAT_INDEX_BASE,
        "moment_axes": list(MOMENT_AXES),
        "mask_axes": list(MASK_AXES),
        **physical_velocity_provenance(
            velocity_estimation_method=source.velocity_estimation_method,
            band_ratio_frequency_scale_hz=(
                source.band_ratio_frequency_scale_hz
            ),
        ),
        "band_lf_source_path": (
            f"/{HD_BAND_LF_PATH}"
            if source.velocity_estimation_method == "frequency_bands"
            else None
        ),
        "band_hf_source_path": (
            f"/{HD_BAND_HF_PATH}"
            if source.velocity_estimation_method == "frequency_bands"
            else None
        ),
        "frequency_band_axes": list(MOMENT_AXES),
    }
