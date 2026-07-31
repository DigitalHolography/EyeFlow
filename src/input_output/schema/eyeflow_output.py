"""Versioned EyeFlow output HDF5 paths."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import h5py

ANGIOEYE_FULL_OUTPUT_SCHEMA = "angioeye_full"
SLIM_TEMP_OUTPUT_SCHEMA = "slim_temp"
EYEFLOW_V2_OUTPUT_SCHEMA = "eyeflow_v2"
ACTIVE_OUTPUT_SCHEMA_VARIANT = EYEFLOW_V2_OUTPUT_SCHEMA


@dataclass(frozen=True)
class DopplerViewAnalysisOutputPaths:
    retinal_velocity_array: str | None
    retinal_artery_velocity_signal: str
    retinal_vein_velocity_signal: str
    retinal_artery_velocity_signal_filtered: str
    retinal_vein_velocity_signal_filtered: str
    velocity_map_avg: str | None
    fRMS_avg: str
    fRMS_bkg_avg: str
    velocitysignal_per_beat: str | None
    velocitysignal_filtered: str | None
    beat_indices: str
    time_per_beat: str


@dataclass(frozen=True)
class VelocityPerBeatOutputPaths:
    velocity_signal: str
    velocity_signal_fft_abs: str
    velocity_signal_fft_arg: str
    velocity_signal_band_limited: str
    vmax_band_limited: str
    vmin_band_limited: str
    vti_per_beat: str
    segment_velocity_signal: str | None = None
    segment_velocity_signal_band_limited: str | None = None


@dataclass(frozen=True)
class OpticDiscTopologyOutputPaths:
    mask: str
    center_xy: str


@dataclass(frozen=True)
class VesselTopologyOutputPaths:
    mask: str | None
    branch_label_map: str
    branch_ids: str
    segment_center_xy: str


@dataclass(frozen=True)
class TopologyOutputPaths:
    optic_disc: OpticDiscTopologyOutputPaths
    artery: VesselTopologyOutputPaths
    vein: VesselTopologyOutputPaths


@dataclass(frozen=True)
class CrossSectionProfileVariantOutputPaths:
    velocity_profile: str
    transverse_coordinate_micrometers: str


@dataclass(frozen=True)
class CrossSectionProfileOutputPaths:
    raw_profile: CrossSectionProfileVariantOutputPaths
    interpolated_profile: CrossSectionProfileVariantOutputPaths

    @property
    def velocity_profile(self) -> str:
        """Compatibility path for the previously exported profile."""
        return self.interpolated_profile.velocity_profile

    @property
    def transverse_coordinate_micrometers(self) -> str:
        """Compatibility path for the previously exported coordinates."""
        return self.interpolated_profile.transverse_coordinate_micrometers


@dataclass(frozen=True)
class HeartbeatOutputPaths:
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
    artery_per_beat: VelocityPerBeatOutputPaths
    vein_per_beat: VelocityPerBeatOutputPaths
    topology: TopologyOutputPaths
    artery_cross_section_profiles: CrossSectionProfileOutputPaths
    vein_cross_section_profiles: CrossSectionProfileOutputPaths
    heartbeat: HeartbeatOutputPaths
    beat_period_seconds: str
    waveform_shape_metrics_root: str
    meta_root: str

    @classmethod
    def active(cls, name: str | None = None) -> "EyeFlowOutputPaths":
        name = ACTIVE_OUTPUT_SCHEMA_VARIANT if name is None else name
        try:
            return OUTPUT_PATH_VARIANTS[name]
        except KeyError as exc:
            known = ", ".join(sorted(OUTPUT_PATH_VARIANTS))
            raise ValueError(
                f"Unknown EyeFlow output schema '{name}'. Known: {known}."
            ) from exc


def _topology_paths(root: str) -> TopologyOutputPaths:
    return TopologyOutputPaths(
        optic_disc=OpticDiscTopologyOutputPaths(
            mask=f"{root}/OpticDisc/Mask/value",
            center_xy=f"{root}/OpticDisc/CenterXY/value",
        ),
        artery=VesselTopologyOutputPaths(
            mask=f"{root}/Artery/Mask/value" if root == "Segmentation" else None,
            branch_label_map=f"{root}/Artery/BranchLabelMap/value",
            branch_ids=f"{root}/Artery/BranchIds/value",
            segment_center_xy=f"{root}/Artery/SegmentCenterXY/value",
        ),
        vein=VesselTopologyOutputPaths(
            mask=f"{root}/Vein/Mask/value" if root == "Segmentation" else None,
            branch_label_map=f"{root}/Vein/BranchLabelMap/value",
            branch_ids=f"{root}/Vein/BranchIds/value",
            segment_center_xy=f"{root}/Vein/SegmentCenterXY/value",
        ),
    )


def _cross_section_profile_paths(
    root: str,
    *,
    velocity_profile_name: str = "VelocityProfile",
) -> CrossSectionProfileOutputPaths:
    def variant_paths(prefix: str) -> CrossSectionProfileVariantOutputPaths:
        variant_root = f"{root}/{prefix}"
        return CrossSectionProfileVariantOutputPaths(
            velocity_profile=f"{variant_root}/{velocity_profile_name}/value",
            transverse_coordinate_micrometers=(
                f"{variant_root}/TransverseCoordinateMicrometers/value"
            ),
        )

    return CrossSectionProfileOutputPaths(
        raw_profile=variant_paths("RawProfile"),
        interpolated_profile=variant_paths("InterpolatedProfile"),
    )


LEGACY_HEARTBEAT_OUTPUT = HeartbeatOutputPaths(
    systolic_peak_frame_indices="analysis/heartbeat/systolic_peak_frame_indices",
    systolic_cycle_duration_seconds="analysis/heartbeat/systolic_cycle_duration_seconds",
    spectral_fundamental_frequency_hz="analysis/heartbeat/spectral_fundamental_frequency_hz",
    spectral_heart_rate_bpm="analysis/heartbeat/spectral_heart_rate_bpm",
    spectral_heart_rate_standard_error_bpm=(
        "analysis/heartbeat/spectral_heart_rate_standard_error_bpm"
    ),
    spectral_period_seconds="analysis/heartbeat/spectral_period_seconds",
)

HEARTBEAT_OUTPUT = HeartbeatOutputPaths(
    systolic_peak_frame_indices="Processing/Heartbeat/Systole/PeakFrameIndices/value",
    systolic_cycle_duration_seconds=(
        "Processing/Heartbeat/Systole/CycleDurationSeconds/value"
    ),
    spectral_fundamental_frequency_hz=(
        "Processing/Heartbeat/Spectral/FundamentalFrequencyHz/value"
    ),
    spectral_heart_rate_bpm="Processing/Heartbeat/Spectral/HeartRateBpm/value",
    spectral_heart_rate_standard_error_bpm=(
        "Processing/Heartbeat/Spectral/HeartRateStandardErrorBpm/value"
    ),
    spectral_period_seconds="Processing/Heartbeat/Spectral/PeriodSeconds/value",
)


ANGIOEYE_FULL_OUTPUT = EyeFlowOutputPaths(
    name=ANGIOEYE_FULL_OUTPUT_SCHEMA,
    analysis=DopplerViewAnalysisOutputPaths(
        retinal_velocity_array="analysis/retinal_velocity_array",
        retinal_artery_velocity_signal="analysis/retinal_artery_velocity_signal",
        retinal_vein_velocity_signal="analysis/retinal_vein_velocity_signal",
        retinal_artery_velocity_signal_filtered="analysis/velocitysignal_filtered",
        retinal_vein_velocity_signal_filtered="analysis/vein_velocitysignal_filtered",
        velocity_map_avg="analysis/velocity_map_avg",
        fRMS_avg="analysis/fRMS_avg",
        fRMS_bkg_avg="analysis/fRMS_bkg_avg",
        velocitysignal_per_beat="analysis/velocitysignal_per_beat",
        velocitysignal_filtered="analysis/velocitysignal_filtered",
        beat_indices="analysis/beat_indices",
        time_per_beat="analysis/time_per_beat",
    ),
    artery_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="Artery/VelocityPerBeat/VelocitySignalPerBeat/value",
        velocity_signal_fft_abs="Artery/VelocityPerBeat/VelocitySignalPerBeatFFT_abs/value",
        velocity_signal_fft_arg="Artery/VelocityPerBeat/VelocitySignalPerBeatFFT_arg/value",
        velocity_signal_band_limited=(
            "Artery/VelocityPerBeat/VelocitySignalPerBeatBandLimited/value"
        ),
        vmax_band_limited="Artery/VelocityPerBeat/VmaxPerBeatBandLimited/value",
        vmin_band_limited="Artery/VelocityPerBeat/VminPerBeatBandLimited/value",
        vti_per_beat="Artery/VelocityPerBeat/VTIPerBeat/value",
        segment_velocity_signal=(
            "Artery/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegment/value"
        ),
        segment_velocity_signal_band_limited=(
            "Artery/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegmentBandLimited/value"
        ),
    ),
    vein_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="Vein/VelocityPerBeat/VelocitySignalPerBeat/value",
        velocity_signal_fft_abs="Vein/VelocityPerBeat/VelocitySignalPerBeatFFT_abs/value",
        velocity_signal_fft_arg="Vein/VelocityPerBeat/VelocitySignalPerBeatFFT_arg/value",
        velocity_signal_band_limited=(
            "Vein/VelocityPerBeat/VelocitySignalPerBeatBandLimited/value"
        ),
        vmax_band_limited="Vein/VelocityPerBeat/VmaxPerBeatBandLimited/value",
        vmin_band_limited="Vein/VelocityPerBeat/VminPerBeatBandLimited/value",
        vti_per_beat="Vein/VelocityPerBeat/VTIPerBeat/value",
        segment_velocity_signal=(
            "Vein/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegment/value"
        ),
        segment_velocity_signal_band_limited=(
            "Vein/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegmentBandLimited/value"
        ),
    ),
    topology=_topology_paths("Topology"),
    artery_cross_section_profiles=_cross_section_profile_paths(
        "Artery/CrossSections",
        velocity_profile_name="VelocityProfileSeg",
    ),
    vein_cross_section_profiles=_cross_section_profile_paths(
        "Vein/CrossSections",
        velocity_profile_name="VelocityProfileSeg",
    ),
    heartbeat=LEGACY_HEARTBEAT_OUTPUT,
    beat_period_seconds="Artery/VelocityPerBeat/beatPeriodSeconds/value",
    waveform_shape_metrics_root="Metrics/waveform_shape_metrics",
    meta_root="Meta",
)


SLIM_TEMP_OUTPUT = EyeFlowOutputPaths(
    name=SLIM_TEMP_OUTPUT_SCHEMA,
    analysis=DopplerViewAnalysisOutputPaths(
        retinal_velocity_array=None,
        retinal_artery_velocity_signal="artery/velocity/signal/value",
        retinal_vein_velocity_signal="vein/velocity/signal/value",
        retinal_artery_velocity_signal_filtered="artery/velocity/filtered_signal/value",
        retinal_vein_velocity_signal_filtered="vein/velocity/filtered_signal/value",
        velocity_map_avg="topo/velocity_map_avg/value",
        fRMS_avg="topo/fRMS_avg/value",
        fRMS_bkg_avg="topo/fRMS_bkg_avg/value",
        velocitysignal_per_beat="artery/velocity/perbeat/filtered_signal/value",
        velocitysignal_filtered="artery/velocity/filtered_signal/value",
        beat_indices="perbeat/beat_indices/value",
        time_per_beat="perbeat/time_per_beat/value",
    ),
    artery_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="artery/velocity/perbeat/signal/value",
        velocity_signal_fft_abs="artery/velocity/perbeat/fft_abs/value",
        velocity_signal_fft_arg="artery/velocity/perbeat/fft_arg/value",
        velocity_signal_band_limited="artery/velocity/perbeat/band_limited/value",
        vmax_band_limited="artery/velocity/perbeat/vmax_band_limited/value",
        vmin_band_limited="artery/velocity/perbeat/vmin_band_limited/value",
        vti_per_beat="artery/velocity/perbeat/vti/value",
        segment_velocity_signal="artery/velocity/perbeat/segments/signal/value",
        segment_velocity_signal_band_limited=(
            "artery/velocity/perbeat/segments/band_limited/value"
        ),
    ),
    vein_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="vein/velocity/perbeat/signal/value",
        velocity_signal_fft_abs="vein/velocity/perbeat/fft_abs/value",
        velocity_signal_fft_arg="vein/velocity/perbeat/fft_arg/value",
        velocity_signal_band_limited="vein/velocity/perbeat/band_limited/value",
        vmax_band_limited="vein/velocity/perbeat/vmax_band_limited/value",
        vmin_band_limited="vein/velocity/perbeat/vmin_band_limited/value",
        vti_per_beat="vein/velocity/perbeat/vti/value",
        segment_velocity_signal="vein/velocity/perbeat/segments/signal/value",
        segment_velocity_signal_band_limited=(
            "vein/velocity/perbeat/segments/band_limited/value"
        ),
    ),
    topology=_topology_paths("Topology"),
    artery_cross_section_profiles=_cross_section_profile_paths(
        "artery/cross_sections"
    ),
    vein_cross_section_profiles=_cross_section_profile_paths(
        "vein/cross_sections"
    ),
    heartbeat=LEGACY_HEARTBEAT_OUTPUT,
    beat_period_seconds="perbeat/beat_period_seconds/value",
    waveform_shape_metrics_root="Metrics/waveform_shape_metrics",
    meta_root="Meta",
)


EYEFLOW_V2_OUTPUT = EyeFlowOutputPaths(
    name=EYEFLOW_V2_OUTPUT_SCHEMA,
    analysis=DopplerViewAnalysisOutputPaths(
        retinal_velocity_array=None,
        retinal_artery_velocity_signal="Processing/Velocity/Artery/Raw/value",
        retinal_vein_velocity_signal="Processing/Velocity/Vein/Raw/value",
        retinal_artery_velocity_signal_filtered=(
            "Processing/Velocity/Artery/Filtered/value"
        ),
        retinal_vein_velocity_signal_filtered=(
            "Processing/Velocity/Vein/Filtered/value"
        ),
        velocity_map_avg=None,
        fRMS_avg="Processing/FrequencyMaps/fRMS_avg/value",
        fRMS_bkg_avg="Processing/FrequencyMaps/fRMS_bkg_avg/value",
        velocitysignal_per_beat=None,
        velocitysignal_filtered=None,
        beat_indices="Processing/Heartbeat/Systole/PeakFrameIndices/value",
        time_per_beat="Processing/Heartbeat/Systole/CycleDurationSeconds/value",
    ),
    artery_per_beat=VelocityPerBeatOutputPaths(
        velocity_signal="Processing/VelocityPerBeat/Artery/Raw/value",
        velocity_signal_fft_abs="Processing/VelocityPerBeat/Artery/FFTAbs/value",
        velocity_signal_fft_arg="Processing/VelocityPerBeat/Artery/FFTPhase/value",
        velocity_signal_band_limited=(
            "Processing/VelocityPerBeat/Artery/BandLimited/value"
        ),
        vmax_band_limited="Processing/VelocityPerBeat/Artery/VmaxBandLimited/value",
        vmin_band_limited="Processing/VelocityPerBeat/Artery/VminBandLimited/value",
        vti_per_beat="Processing/VelocityPerBeat/Artery/VTI/value",
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
        vmax_band_limited="Processing/VelocityPerBeat/Vein/VmaxBandLimited/value",
        vmin_band_limited="Processing/VelocityPerBeat/Vein/VminBandLimited/value",
        vti_per_beat="Processing/VelocityPerBeat/Vein/VTI/value",
        segment_velocity_signal=(
            "Processing/VelocityPerBeat/Vein/Segments/Raw/value"
        ),
        segment_velocity_signal_band_limited=(
            "Processing/VelocityPerBeat/Vein/Segments/BandLimited/value"
        ),
    ),
    topology=_topology_paths("Segmentation"),
    artery_cross_section_profiles=_cross_section_profile_paths(
        "Processing/CrossSections/Artery"
    ),
    vein_cross_section_profiles=_cross_section_profile_paths(
        "Processing/CrossSections/Vein"
    ),
    heartbeat=HEARTBEAT_OUTPUT,
    beat_period_seconds="Processing/VelocityPerBeat/BeatPeriodSeconds/value",
    waveform_shape_metrics_root="Processing/Metrics/waveform_shape_metrics",
    meta_root="Meta",
)


OUTPUT_PATH_VARIANTS = {
    ANGIOEYE_FULL_OUTPUT_SCHEMA: ANGIOEYE_FULL_OUTPUT,
    SLIM_TEMP_OUTPUT_SCHEMA: SLIM_TEMP_OUTPUT,
    EYEFLOW_V2_OUTPUT_SCHEMA: EYEFLOW_V2_OUTPUT,
}

ZERO_BASED_INDEX_PATHS = frozenset(
    paths.analysis.beat_indices for paths in OUTPUT_PATH_VARIANTS.values()
)


def systolic_index_base_for_path(path: str) -> int | None:
    from input_output.writers.h5 import normalize_h5_path

    normalized = normalize_h5_path(path)
    return 0 if normalized in ZERO_BASED_INDEX_PATHS else None


def iter_metric_datasets(group: h5py.Group) -> Iterator[tuple[str, h5py.Dataset]]:
    def visitor(name: str, obj: h5py.Group | h5py.Dataset) -> None:
        if isinstance(obj, h5py.Dataset):
            datasets.append((name, obj))

    datasets: list[tuple[str, h5py.Dataset]] = []
    group.visititems(visitor)
    yield from datasets
