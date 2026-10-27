"""Source assembly for shared velocity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from calculations.vessel_segments.measurement import SegmentMeasurementSettings
from input_output.schema import (
    DopplerViewSource,
    HolodopplerSource,
    PixelPitch,
    RetinalSourceData,
)
from input_output.canonical_sources import load_retinal_source_data
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
class VelocityAnalysisSourceData:
    """Resolved source data with explicit analysis-axis contracts."""

    source: RetinalSourceData
    profile_settings: SegmentMeasurementSettings
    provenance: dict[str, object]


@dataclass(frozen=True)
class VelocityAnalysisSources:
    """Typed input adapters needed by velocity analysis."""

    hd: HolodopplerSource
    dv: DopplerViewSource
    velocity_estimation_method: str = "doppler_moments"
    band_ratio_frequency_scale_hz: float = (
        DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    )

    @classmethod
    def from_context(cls, ctx: PipelineContext) -> VelocityAnalysisSources:
        ctx.require_inputs("hd", "dv")
        return cls(
            hd=ctx.inputs.hd.as_holodoppler(),
            dv=ctx.inputs.dv.as_dopplerview(),
            velocity_estimation_method=ctx.velocity_estimation_method,
            band_ratio_frequency_scale_hz=ctx.band_ratio_frequency_scale_hz,
        )

    def load(self, source: RetinalSourceData | None = None) -> VelocityAnalysisSourceData:
        if source is None:
            source = load_retinal_source_data(
                self.hd,
                self.dv,
                velocity_estimation_method=self.velocity_estimation_method,
                include_spectrograms=False,
                band_ratio_frequency_scale_hz=self.band_ratio_frequency_scale_hz,
            )
        elif source.velocity_estimation_method != self.velocity_estimation_method:
            raise ValueError("Cached retinal source uses a different velocity method.")
        pixel_pitch = source.holodoppler.pixel_pitch
        return VelocityAnalysisSourceData(
            source=source,
            profile_settings=self._profile_settings(pixel_pitch),
            provenance=_source_provenance(self.hd, self.dv, source),
        )

    def _profile_settings(self, pixel_pitch: PixelPitch) -> SegmentMeasurementSettings:
        return SegmentMeasurementSettings(
            pixel_size_mm=self._pixel_size(pixel_pitch),
            submask_size_percentile_kept=(CROSS_SECTION_SUBMASK_SIZE_PERCENTILE_KEPT),
        )

    def _pixel_size(self, pixel_pitch: PixelPitch) -> float:
        return pixel_pitch.isotropic_mm


def load_velocity_analysis_source_data(ctx: PipelineContext) -> VelocityAnalysisSourceData:
    cached = ctx.state.get("velocity_source")
    return VelocityAnalysisSources.from_context(ctx).load(
        cached if isinstance(cached, RetinalSourceData) else None
    )


def _source_provenance(
    hd,
    dv,
    source: RetinalSourceData,
) -> dict[str, object]:
    segmentation = source.segmentation
    optic_disc = segmentation.optic_disc
    image_maps = source.image_maps
    if source.velocity_estimation_method == "frequency_bands":
        source_attrs = {
            "band_lf_source_path": image_maps.band_lf.name,
            "band_hf_source_path": image_maps.band_hf.name,
            "frequency_band_axes": list(MOMENT_AXES),
        }
    else:
        source_attrs = {
            "moment0_source_path": image_maps.moment0.name,
            "moment2_source_path": image_maps.moment2.name,
        }
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
        **source_attrs,
    }
