"""Prepare and analyze velocity segments for configured vessel classes."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping

from calculations.topology import AnnulusGeometry, OpticDisc
from calculations.vessel_segments.sampling.cache import TopologyCacheKey
from calculations.vessel_segments.sampling.models import SegmentSamplingPlan
from calculations.vessel_segments.profiles.fft import (
    DEFAULT_PROFILE_MASK_DILATION_PIXELS,
    SegmentFftAccumulator,
)
from calculations.vessel_segments.measurement import (
    SegmentMeasurements,
    SegmentMeasurementSettings,
    analyze_segment_profiles,
)
from utils.logger import Logger

from ..models import VelocitySegmentResult


def analyze_velocity_segment_profiles(
    velocity_map,
    vessel_masks: Mapping[str, object],
    optic_disc: OpticDisc,
    ring_settings: AnnulusGeometry,
    profile_settings: SegmentMeasurementSettings,
    *,
    source_id: str = "",
    topology_cache: MutableMapping[TopologyCacheKey, object] | None = None,
    retain_velocity_maps: bool = False,
    cycle_boundary_indexes=None,
    velocity_profile_fft: bool = False,
    index_base: int = 0,
    transverse_mask_dilation_pixels: int | None = None,
    prepared_topologies: Mapping[str, SegmentSamplingPlan] | None = None,
    transform_mode: str = "fused",
    post_interpolation=None,
    temporal_halo: int = 0,
    scratch_array_count: int | None = None,
) -> dict[str, VelocitySegmentResult]:
    """Measure velocity profiles and attach velocity-only optional products."""

    if not vessel_masks:
        return {}
    if velocity_profile_fft and cycle_boundary_indexes is None:
        raise ValueError(
            "cycle_boundary_indexes are required for velocity FFT profiles."
    )
    fft_profiles: dict[str, SegmentFftAccumulator] = {}

    def segment_observer_factory(name, topology):
        geometry = topology.native
        if not velocity_profile_fft:
            return None
        accumulator = SegmentFftAccumulator(
            frame_count=int(velocity_map.shape[0]),
            ring_count=int(geometry.annulus_masks.shape[0]),
            branch_count=int(geometry.branch_ids.size),
            canvas_side=int(topology.rotated_masks.shape[-1]),
            cycle_boundary_indexes=cycle_boundary_indexes,
            index_base=index_base,
        )
        fft_profiles[name] = accumulator
        return accumulator.observe

    dilation_pixels = (
        {
            str(name): _legacy_profile_dilation_pixels(str(name))
            for name in vessel_masks
        }
        if transverse_mask_dilation_pixels is None
        else int(transverse_mask_dilation_pixels)
    )
    profile_results = analyze_segment_profiles(
        velocity_map,
        vessel_masks,
        optic_disc,
        ring_settings,
        profile_settings,
        source_id=source_id,
        topology_cache=topology_cache,
        retain_segment_maps=retain_velocity_maps,
        transverse_mask_dilation_pixels=dilation_pixels,
        prepared_topologies=prepared_topologies,
        transform_mode=transform_mode,
        post_interpolation=post_interpolation,
        temporal_halo=temporal_halo,
        scratch_array_count=scratch_array_count,
        segment_observer_factory=segment_observer_factory,
    )

    results: dict[str, VelocitySegmentResult] = {}
    for name, profile_result in profile_results.items():
        accumulator = fft_profiles.get(name)
        results[name] = _velocity_result(
            profile_result,
            fft_profiles=accumulator,
        )
        if accumulator is not None:
            Logger.log(
                f"Completed {name} optional velocity-profile FFT in "
                f"{accumulator.elapsed_seconds:.2f}s."
            )
    return results


def _velocity_result(
    profiles: SegmentMeasurements,
    *,
    fft_profiles: SegmentFftAccumulator | None,
) -> VelocitySegmentResult:
    """Attach velocity-only optional products to neutral segment profiles."""

    return VelocitySegmentResult.from_profile_result(
        profiles,
        transverse_fft_profiles_unmasked=(
            fft_profiles.unmasked if fft_profiles else None
        ),
        transverse_fft_profiles_masked=(
            fft_profiles.masked if fft_profiles else None
        ),
    )


def _legacy_profile_dilation_pixels(vessel_name: str) -> int:
    """Retain the historical artery-only transverse profile expansion."""

    return (
        DEFAULT_PROFILE_MASK_DILATION_PIXELS
        if str(vessel_name).lower() == "artery"
        else 0
    )
