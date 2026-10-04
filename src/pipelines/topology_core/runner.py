"""Prepare map-independent vessel topology once per pipeline run."""

from __future__ import annotations

from collections.abc import Mapping

from calculations.topology import (
    PreparedTopology,
    prepare_topologies,
    resolve_segment_rotations,
    run_topology_cache,
    topology_source_id,
)
from pipeline_engine.imports import read_int_setting
from pipelines.vessel_inputs import load_vessel_topology_inputs

TOPOLOGY_CORE_STATE = "topology_core.prepared"
_DEFAULT_NUMBER_OF_RADII_IN_FOV = 25
_TOPOLOGY_WINDOW_SIZE_PERCENTILE_KEPT = 0.95


def run_topology_core(ctx) -> dict[str, object]:
    """Prepare and publish the canonical artery and vein topology."""

    inputs = load_vessel_topology_inputs(ctx)
    images = inputs.image_maps
    segmentation = inputs.segmentation
    vessels = segmentation.vessels
    number_of_radii = read_int_setting(
        ctx,
        default=_DEFAULT_NUMBER_OF_RADII_IN_FOV,
        keys=(
            "number_of_radii_in_FOV",
            "number_of_radii_in_fov",
            "NumberOfRadiiInFOV",
            "number_of_radii_over_FOV",
            "number_of_radii_over_fov",
            "NumberOfRadiiOverFOV",
        ),
    )
    if number_of_radii < 1:
        raise ValueError("number_of_radii_in_FOV must be positive.")

    topology_reference = (
        images.band_lf
        if inputs.velocity_estimation_method == "frequency_bands"
        else images.moment0
    )
    ring_settings = segmentation.optic_disc.annulus_geometry(
        tuple(int(size) for size in topology_reference.shape[-2:]),
        number_of_radii_in_fov=number_of_radii,
    )
    prepared = prepare_topologies(
        {
            "artery": vessels.artery,
            "vein": vessels.vein,
        },
        segmentation.optic_disc,
        ring_settings,
        source_id=topology_source_id(
            ctx.inputs.hd.filename,
            ctx.inputs.dv.filename,
        ),
        cache=run_topology_cache(ctx.state.raw),
        window_size_percentile_kept=(_TOPOLOGY_WINDOW_SIZE_PERCENTILE_KEPT),
    )
    prepared = {
        name: resolve_segment_rotations(topology, topology_reference)
        for name, topology in prepared.items()
    }
    ctx.state.set(TOPOLOGY_CORE_STATE, prepared)
    # Import lazily to retain the legacy module as a compatibility surface
    # without coupling topology registration to the waveform package import.
    from pipelines.waveform_velocity_core.segmentation import pack_topology_outputs

    return pack_topology_outputs(
        vessels.artery,
        vessels.vein,
        segmentation.optic_disc,
        prepared,
        inputs.holodoppler.pixel_pitch,
    )


def prepared_topologies(ctx) -> Mapping[str, PreparedTopology]:
    """Return the canonical topology produced by the declared DAG dependency."""

    value = ctx.state.get(TOPOLOGY_CORE_STATE)
    if not isinstance(value, Mapping) or set(value) != {"artery", "vein"}:
        raise RuntimeError("Prepared topology state is unavailable; check the pipeline DAG.")
    if not all(isinstance(item, PreparedTopology) for item in value.values()):
        raise TypeError("Prepared topology state contains an invalid value.")
    return value


__all__ = [
    "TOPOLOGY_CORE_STATE",
    "prepared_topologies",
    "run_topology_core",
]
