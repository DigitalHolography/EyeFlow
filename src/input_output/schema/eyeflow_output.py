"""EyeFlow output HDF5 paths."""

from __future__ import annotations

from dataclasses import dataclass

EYEFLOW_V2_OUTPUT_SCHEMA = "eyeflow_v2"

VELOCITY_WORKFLOW_ROOTS = {
    "doppler_moments": "Processing",
    "frequency_bands": "ProcessingAlt",
}
VELOCITY_WORKFLOW_FOLDERS = {
    "doppler_moments": "moments",
    "frequency_bands": "bandratio",
}


def processing_path(path: str, root: str) -> str:
    """Resolve a canonical processing path for one velocity workflow."""
    normalized = str(path).replace("\\", "/").strip("/")
    if normalized == "Processing" or normalized.startswith("Processing/"):
        return root + normalized[len("Processing"):]
    return normalized


@dataclass(frozen=True)
class DopplerViewAnalysisOutputPaths:
    retinal_artery_velocity_signal: str
    retinal_vein_velocity_signal: str
    retinal_artery_velocity_signal_band_limited: str
    retinal_vein_velocity_signal_band_limited: str
    velocity_map_avg_masked: str
    fRMS_avg: str
    fRMS_bkg_avg: str
    beat_indices: str
    time_per_beat: str


@dataclass(frozen=True)
class SegmentVelocityOutputPaths:
    velocity_signal: str | None
    velocity_signal_band_limited: str | None = None
    velocity_map_per_segment: str | None = None
    segments: str | None = None


@dataclass(frozen=True)
class VelocityPerBeatOutputPaths:
    velocity_signal: str
    velocity_signal_fft_abs: str
    velocity_signal_fft_arg: str
    velocity_signal_band_limited: str
    segment_velocity_signal: str | None = None
    segment_velocity_signal_band_limited: str | None = None


@dataclass(frozen=True)
class OpticDiscSegmentationOutputPaths:
    mask: str
    height: str
    width: str
    center: str


@dataclass(frozen=True)
class VesselSegmentationOutputPaths:
    mask: str
    branch_label_map: str
    segment_map: str
    segment_mask_area: str
    lumen_diameter: str
    delta_radius: str


@dataclass(frozen=True)
class SegmentationOutputPaths:
    optic_disc: OpticDiscSegmentationOutputPaths
    artery: VesselSegmentationOutputPaths
    vein: VesselSegmentationOutputPaths
    pixel_pitch_m: str


@dataclass(frozen=True)
class VesselBloodVolumeRateOutputPaths:
    dynamic_edges: str
    static_edges: str
    masked_edges: str
    total_masked_edges: str


@dataclass(frozen=True)
class BloodVolumeRateOutputPaths:
    artery: VesselBloodVolumeRateOutputPaths
    vein: VesselBloodVolumeRateOutputPaths


@dataclass(frozen=True)
class VelocityProfileOutputPaths:
    transverse_velocity_profile_unmasked: str
    transverse_velocity_profile_masked: str
    longitudinal_velocity_profile_unmasked: str
    longitudinal_velocity_profile_masked: str
    transverse_velocity_profile_unmasked_meaned: str
    transverse_velocity_profile_masked_meaned: str
    longitudinal_velocity_profile_unmasked_meaned: str
    longitudinal_velocity_profile_masked_meaned: str
    transverse_velocity_profile_fft_unmasked: str | None = None
    transverse_velocity_profile_fft_masked: str | None = None

@dataclass(frozen=True)
class CardiacCycleOutputPaths:
    systolic_peak_frame_indices: str
    systolic_cycle_duration_seconds: str
    spectral_fundamental_frequency_hz: str
    spectral_heart_rate_bpm: str
    spectral_heart_rate_standard_error_bpm: str
    spectral_period_seconds: str


@dataclass(frozen=True)
class EyeFlowOutputPaths:
    name: str
    analysis: DopplerViewAnalysisOutputPaths
    artery_segments: SegmentVelocityOutputPaths
    vein_segments: SegmentVelocityOutputPaths
    artery_per_beat: VelocityPerBeatOutputPaths
    vein_per_beat: VelocityPerBeatOutputPaths
    artery_per_beat_safe: SegmentVelocityOutputPaths
    vein_per_beat_safe: SegmentVelocityOutputPaths
    segmentation: SegmentationOutputPaths
    artery_velocity_profiles: VelocityProfileOutputPaths
    vein_velocity_profiles: VelocityProfileOutputPaths
    cardiac_cycle: CardiacCycleOutputPaths
    displacement_map: str
    spatial_gradient_profiles_root: str
    spatial_gradient_metrics_root: str
    velocity_profile_analysis_root: str
    waveform_shape_metrics_root: str
    absolute_waveform_metrics_root: str
    lowrank_waveform_decomposition_root: str
    meta_root: str
    blood_volume_rate: BloodVolumeRateOutputPaths

    @staticmethod
    def active(name: str | None = None) -> "EyeFlowOutputPaths":
        if name is not None and name != EYEFLOW_V2_OUTPUT_SCHEMA:
            raise ValueError(
                f"Unknown EyeFlow output schema '{name}'. "
                f"Known: {EYEFLOW_V2_OUTPUT_SCHEMA}."
            )
        return EYEFLOW_V2_OUTPUT


def _segmentation_paths(root: str) -> SegmentationOutputPaths:
    return SegmentationOutputPaths(
        optic_disc=OpticDiscSegmentationOutputPaths(
            mask=f"{root}/OpticDisc/Mask/value",
            height=f"{root}/OpticDisc/Height/value",
            width=f"{root}/OpticDisc/Width/value",
            center=f"{root}/OpticDisc/Center/value",
        ),
        artery=VesselSegmentationOutputPaths(
            mask=f"{root}/Artery/Mask/value",
            branch_label_map=f"{root}/Artery/BranchLabelMap/value",
            segment_map=f"{root}/Artery/SegmentMap/value",
            segment_mask_area=f"{root}/Artery/SegmentMaskArea/value",
            lumen_diameter=f"{root}/Artery/LumenDiameter/value",
            delta_radius=f"{root}/Artery/DeltaRadius/value",
        ),
        vein=VesselSegmentationOutputPaths(
            mask=f"{root}/Vein/Mask/value",
            branch_label_map=f"{root}/Vein/BranchLabelMap/value",
            segment_map=f"{root}/Vein/SegmentMap/value",
            segment_mask_area=f"{root}/Vein/SegmentMaskArea/value",
            lumen_diameter=f"{root}/Vein/LumenDiameter/value",
            delta_radius=f"{root}/Vein/DeltaRadius/value",
        ),
        pixel_pitch_m=f"{root}/PixelPitch_m/value",
    )


def _blood_volume_rate_paths(root: str) -> BloodVolumeRateOutputPaths:
    def vessel(name: str) -> VesselBloodVolumeRateOutputPaths:
        vessel_root = f"{root}/{name}"
        return VesselBloodVolumeRateOutputPaths(
            dynamic_edges=f"{vessel_root}/dynamicEdges/value",
            static_edges=f"{vessel_root}/staticEdges/value",
            masked_edges=f"{vessel_root}/maskedEdges/value",
            total_masked_edges=f"{vessel_root}/totalMaskedEdges/value",
        )

    return BloodVolumeRateOutputPaths(
        artery=vessel("Artery"),
        vein=vessel("Vein"),
    )


def _velocity_profile_paths(
    root: str,
    *,
    fft_root: str,
) -> VelocityProfileOutputPaths:
    return VelocityProfileOutputPaths(
        transverse_velocity_profile_fft_unmasked=(
            f"{fft_root}/TransverseVelocityProfileUnmasked"
        ),
        transverse_velocity_profile_fft_masked=(
            f"{fft_root}/TransverseVelocityProfileMasked"
        ),
        transverse_velocity_profile_unmasked=(
            f"{root}/Transversal/Unmasked/VelocityProfile/value"
        ),
        transverse_velocity_profile_masked=(
            f"{root}/Transversal/Masked/VelocityProfile/value"
        ),
        longitudinal_velocity_profile_unmasked=(
            f"{root}/Longitudinal/Unmasked/VelocityProfile/value"
        ),
        longitudinal_velocity_profile_masked=(
            f"{root}/Longitudinal/Masked/VelocityProfile/value"
        ),
        transverse_velocity_profile_unmasked_meaned=(
            f"{root}/Transversal/Unmasked/VelocityProfileMeaned/value"
        ),
        transverse_velocity_profile_masked_meaned=(
            f"{root}/Transversal/Masked/VelocityProfileMeaned/value"
        ),
        longitudinal_velocity_profile_unmasked_meaned=(
            f"{root}/Longitudinal/Unmasked/VelocityProfileMeaned/value"
        ),
        longitudinal_velocity_profile_masked_meaned=(
            f"{root}/Longitudinal/Masked/VelocityProfileMeaned/value"
        ),
    )


CARDIAC_CYCLE_OUTPUT = CardiacCycleOutputPaths(
    systolic_peak_frame_indices="Processing/CardiacCycle/Systole/PeakFrameIndices/value",
    systolic_cycle_duration_seconds=(
        "Processing/CardiacCycle/Systole/CycleDurationSeconds/value"
    ),
    spectral_fundamental_frequency_hz=(
        "Processing/CardiacCycle/Spectral/FundamentalFrequencyHz/value"
    ),
    spectral_heart_rate_bpm="Processing/CardiacCycle/Spectral/HeartRateBpm/value",
    spectral_heart_rate_standard_error_bpm=(
        "Processing/CardiacCycle/Spectral/HeartRateStandardErrorBpm/value"
    ),
    spectral_period_seconds="Processing/CardiacCycle/Spectral/PeriodSeconds/value",
)


EYEFLOW_V2_OUTPUT = EyeFlowOutputPaths(
    name=EYEFLOW_V2_OUTPUT_SCHEMA,
    analysis=DopplerViewAnalysisOutputPaths(
        retinal_artery_velocity_signal="Processing/Velocity/global/Artery/Raw/value",
        retinal_vein_velocity_signal="Processing/Velocity/global/Vein/Raw/value",
        retinal_artery_velocity_signal_band_limited=(
            "Processing/Velocity/global/Artery/BandLimited/value"
        ),
        retinal_vein_velocity_signal_band_limited=(
            "Processing/Velocity/global/Vein/BandLimited/value"
        ),
        velocity_map_avg_masked="Processing/Maps/VelocityAverageMasked/value",
        fRMS_avg="Processing/Maps/FRMSAverage/value",
        fRMS_bkg_avg="Processing/Maps/FRMSBackgroundAverage/value",
        beat_indices="Processing/CardiacCycle/Systole/PeakFrameIndices/value",
        time_per_beat="Processing/CardiacCycle/Systole/CycleDurationSeconds/value",
    ),
    artery_segments=SegmentVelocityOutputPaths(
        velocity_signal="Processing/Velocity/segments/Artery/Raw/value",
        velocity_signal_band_limited=(
            "Processing/Velocity/segments/Artery/BandLimited/value"
        ),
        velocity_map_per_segment="Processing/VelocityMapPerSegment/Artery",
        segments="Segmentation/Artery/Segments",
    ),
    vein_segments=SegmentVelocityOutputPaths(
        velocity_signal="Processing/Velocity/segments/Vein/Raw/value",
        velocity_signal_band_limited=(
            "Processing/Velocity/segments/Vein/BandLimited/value"
        ),
        velocity_map_per_segment="Processing/VelocityMapPerSegment/Vein",
        segments="Segmentation/Vein/Segments",
    ),
    artery_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="Processing/VelocityPerBeat/Artery/Raw/value",
        velocity_signal_fft_abs="Processing/VelocityPerBeat/Artery/FFTAbs/value",
        velocity_signal_fft_arg="Processing/VelocityPerBeat/Artery/FFTPhase/value",
        velocity_signal_band_limited=(
            "Processing/VelocityPerBeat/Artery/BandLimited/value"
        ),
        segment_velocity_signal=(
            "Processing/VelocityPerBeat/Artery/Segments/Raw/value"
        ),
        segment_velocity_signal_band_limited=(
            "Processing/VelocityPerBeat/Artery/Segments/BandLimited/value"
        ),
    ),
    vein_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="Processing/VelocityPerBeat/Vein/Raw/value",
        velocity_signal_fft_abs="Processing/VelocityPerBeat/Vein/FFTAbs/value",
        velocity_signal_fft_arg="Processing/VelocityPerBeat/Vein/FFTPhase/value",
        velocity_signal_band_limited=(
            "Processing/VelocityPerBeat/Vein/BandLimited/value"
        ),
        segment_velocity_signal=(
            "Processing/VelocityPerBeat/Vein/Segments/Raw/value"
        ),
        segment_velocity_signal_band_limited=(
            "Processing/VelocityPerBeat/Vein/Segments/BandLimited/value"
        ),
    ),
    artery_per_beat_safe=SegmentVelocityOutputPaths(
        velocity_signal=(
            "Processing/VelocityPerBeatSafe/Artery/Segments/Raw/value"
        ),
        velocity_signal_band_limited=(
            "Processing/VelocityPerBeatSafe/Artery/Segments/BandLimited/value"
        ),
    ),
    vein_per_beat_safe=SegmentVelocityOutputPaths(
        velocity_signal=(
            "Processing/VelocityPerBeatSafe/Vein/Segments/Raw/value"
        ),
        velocity_signal_band_limited=(
            "Processing/VelocityPerBeatSafe/Vein/Segments/BandLimited/value"
        ),
    ),
    segmentation=_segmentation_paths("Segmentation"),
    artery_velocity_profiles=_velocity_profile_paths(
        "Processing/VelocityProfiles/Artery",
        fft_root="Processing/VelocityProfilesFFT/Artery",
    ),
    vein_velocity_profiles=_velocity_profile_paths(
        "Processing/VelocityProfiles/Vein",
        fft_root="Processing/VelocityProfilesFFT/Vein",
    ),
    cardiac_cycle=CARDIAC_CYCLE_OUTPUT,
    displacement_map="Processing/DisplacementMap",
    spatial_gradient_profiles_root="Processing/SpatialGradientProfiles",
    spatial_gradient_metrics_root="Processing/SpatialGradientMetrics",
    velocity_profile_analysis_root="Processing/VelocityProfileAnalysis",
    waveform_shape_metrics_root="Processing/Metrics/waveform_shape_metrics",
    absolute_waveform_metrics_root="Processing/Metrics/absolute_waveform_metrics",
    lowrank_waveform_decomposition_root=(
        "Processing/Metrics/lowrank_waveform_decomposition"
    ),
    meta_root="Meta",
    blood_volume_rate=_blood_volume_rate_paths("Processing/BloodVolumeRate"),
)
