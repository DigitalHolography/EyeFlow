"""EyeFlow output HDF5 paths for AngioEye and temporary schemas."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import h5py

ANGIOEYE_FULL_OUTPUT_SCHEMA = "angioeye_full"
SLIM_TEMP_OUTPUT_SCHEMA = "slim_temp"
ACTIVE_OUTPUT_SCHEMA_VARIANT = ANGIOEYE_FULL_OUTPUT_SCHEMA


@dataclass(frozen=True)
class DopplerViewAnalysisOutputPaths:
    retinal_velocity_array: str | None
    retinal_artery_velocity_signal: str
    retinal_vein_velocity_signal: str
    velocity_map_avg: str
    fRMS_avg: str
    fRMS_bkg_avg: str
    velocitysignal_per_beat: str
    velocitysignal_filtered: str
    beat_indices: str
    time_per_beat: str


@dataclass(frozen=True)
class VelocityPerBeatOutputPaths:
    velocity_signal: str
    velocity_signal_fft_abs: str
    velocity_signal_fft_arg: str
    velocity_signal_band_limited: str
    segment_velocity_signal: str | None = None
    segment_velocity_signal_band_limited: str | None = None


@dataclass(frozen=True)
class BeatQualityOutputPaths:
    candidate_start_index: str
    candidate_stop_index: str
    accepted_mask: str
    accepted_candidate_index: str
    rejection_flags: str
    quality_score: str
    period_ratio: str
    boundary_recovered_mask: str
    nominal_period_samples: str
    artery_raw_bandlimited_residual: str
    vein_raw_bandlimited_residual: str
    artery_template_correlation: str
    vein_template_correlation: str
    candidate_count: str
    accepted_count: str
    rejected_count: str
    recovered_boundary_count: str


@dataclass(frozen=True)
class OpticDiscTopologyOutputPaths:
    mask: str
    center_xy: str


@dataclass(frozen=True)
class VesselTopologyOutputPaths:
    branch_label_map: str
    branch_ids: str
    segment_center_xy: str


@dataclass(frozen=True)
class TopologyOutputPaths:
    optic_disc: OpticDiscTopologyOutputPaths
    artery: VesselTopologyOutputPaths
    vein: VesselTopologyOutputPaths


@dataclass(frozen=True)
class EyeFlowOutputPaths:
    name: str
    analysis: DopplerViewAnalysisOutputPaths
    artery_per_beat: VelocityPerBeatOutputPaths
    vein_per_beat: VelocityPerBeatOutputPaths
    topology: TopologyOutputPaths
    beat_period_idx: str
    beat_period_seconds: str
    waveform_shape_metrics_root: str
    beat_quality: BeatQualityOutputPaths | None = None

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


TOPOLOGY_OUTPUT = TopologyOutputPaths(
    optic_disc=OpticDiscTopologyOutputPaths(
        mask="Topology/OpticDisc/Mask/value",
        center_xy="Topology/OpticDisc/CenterXY/value",
    ),
    artery=VesselTopologyOutputPaths(
        branch_label_map="Topology/Artery/BranchLabelMap/value",
        branch_ids="Topology/Artery/BranchIds/value",
        segment_center_xy="Topology/Artery/SegmentCenterXY/value",
    ),
    vein=VesselTopologyOutputPaths(
        branch_label_map="Topology/Vein/BranchLabelMap/value",
        branch_ids="Topology/Vein/BranchIds/value",
        segment_center_xy="Topology/Vein/SegmentCenterXY/value",
    ),
)


ANGIOEYE_FULL_OUTPUT = EyeFlowOutputPaths(
    name=ANGIOEYE_FULL_OUTPUT_SCHEMA,
    analysis=DopplerViewAnalysisOutputPaths(
        retinal_velocity_array="analysis/retinal_velocity_array",
        retinal_artery_velocity_signal="analysis/retinal_artery_velocity_signal",
        retinal_vein_velocity_signal="analysis/retinal_vein_velocity_signal",
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
        segment_velocity_signal=(
            "Artery/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegment/value"
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
        segment_velocity_signal=(
            "Vein/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegment/value"
        ),
        segment_velocity_signal_band_limited=(
            "Vein/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegmentBandLimited/value"
        ),
    ),
    topology=TOPOLOGY_OUTPUT,
    beat_period_idx="Artery/VelocityPerBeat/beatPeriodIdx/value",
    beat_period_seconds="Artery/VelocityPerBeat/beatPeriodSeconds/value",
    waveform_shape_metrics_root="Metrics/waveform_shape_metrics",
    beat_quality=BeatQualityOutputPaths(
        candidate_start_index="QualityControl/Beat/CandidateStartIndex/value",
        candidate_stop_index="QualityControl/Beat/CandidateStopIndex/value",
        accepted_mask="QualityControl/Beat/AcceptedMask/value",
        accepted_candidate_index="QualityControl/Beat/AcceptedCandidateIndex/value",
        rejection_flags="QualityControl/Beat/RejectionFlags/value",
        quality_score="QualityControl/Beat/QualityScore/value",
        period_ratio="QualityControl/Beat/PeriodRatio/value",
        boundary_recovered_mask="QualityControl/Beat/BoundaryRecoveredMask/value",
        nominal_period_samples="QualityControl/Beat/NominalPeriodSamples/value",
        artery_raw_bandlimited_residual=(
            "QualityControl/Beat/ArteryRawBandlimitedResidual/value"
        ),
        vein_raw_bandlimited_residual=(
            "QualityControl/Beat/VeinRawBandlimitedResidual/value"
        ),
        artery_template_correlation=(
            "QualityControl/Beat/ArteryTemplateCorrelation/value"
        ),
        vein_template_correlation=(
            "QualityControl/Beat/VeinTemplateCorrelation/value"
        ),
        candidate_count="QualityControl/Beat/CandidateCount/value",
        accepted_count="QualityControl/Beat/AcceptedCount/value",
        rejected_count="QualityControl/Beat/RejectedCount/value",
        recovered_boundary_count=(
            "QualityControl/Beat/RecoveredBoundaryCount/value"
        ),
    ),
)

SLIM_TEMP_OUTPUT = EyeFlowOutputPaths(
    name=SLIM_TEMP_OUTPUT_SCHEMA,
    analysis=DopplerViewAnalysisOutputPaths(
        retinal_velocity_array=None,
        retinal_artery_velocity_signal="artery/velocity/signal/value",
        retinal_vein_velocity_signal="vein/velocity/signal/value",
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
        segment_velocity_signal="vein/velocity/perbeat/segments/signal/value",
        segment_velocity_signal_band_limited=(
            "vein/velocity/perbeat/segments/band_limited/value"
        ),
    ),
    topology=TOPOLOGY_OUTPUT,
    beat_period_idx="perbeat/beat_period_idx/value",
    beat_period_seconds="perbeat/beat_period_seconds/value",
    waveform_shape_metrics_root="Metrics/waveform_shape_metrics",
    beat_quality=BeatQualityOutputPaths(
        candidate_start_index="quality/beat/candidate_start_index/value",
        candidate_stop_index="quality/beat/candidate_stop_index/value",
        accepted_mask="quality/beat/accepted_mask/value",
        accepted_candidate_index="quality/beat/accepted_candidate_index/value",
        rejection_flags="quality/beat/rejection_flags/value",
        quality_score="quality/beat/quality_score/value",
        period_ratio="quality/beat/period_ratio/value",
        boundary_recovered_mask="quality/beat/boundary_recovered_mask/value",
        nominal_period_samples="quality/beat/nominal_period_samples/value",
        artery_raw_bandlimited_residual=(
            "quality/beat/artery_raw_bandlimited_residual/value"
        ),
        vein_raw_bandlimited_residual=(
            "quality/beat/vein_raw_bandlimited_residual/value"
        ),
        artery_template_correlation="quality/beat/artery_template_correlation/value",
        vein_template_correlation="quality/beat/vein_template_correlation/value",
        candidate_count="quality/beat/candidate_count/value",
        accepted_count="quality/beat/accepted_count/value",
        rejected_count="quality/beat/rejected_count/value",
        recovered_boundary_count="quality/beat/recovered_boundary_count/value",
    ),
)

OUTPUT_PATH_VARIANTS = {
    ANGIOEYE_FULL_OUTPUT_SCHEMA: ANGIOEYE_FULL_OUTPUT,
    SLIM_TEMP_OUTPUT_SCHEMA: SLIM_TEMP_OUTPUT,
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
    for item in datasets:
        yield item
