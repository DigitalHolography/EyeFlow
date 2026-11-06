"""Prepare map-independent vessel topology once per pipeline run."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping

import numpy as np

from calculations.topology import OpticDisc, capillary_analysis_zones
from calculations.topology.capillary import (
    GAPS_CUTOFF_DIAMETER_WIDTH_FRACTION,
    capillary_disc_support,
)
from calculations.topology.optic_disc import DEFAULT_OPTIC_DISC_DIAMETER_MM
from calculations.vessel_segments.sampling import (
    SegmentSamplingPlan,
    prepare_sampling_plans,
    resolve_segment_rotations,
)
from calculations.vessel_segments.sampling.cache import topology_source_id, TopologyCacheKey
from data_manager import InvocationResource
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.doppler_view import (
    DV_CHOROID_VESSEL_MASK_PATH,
    DV_RETINAL_VESSEL_LOGITS_PATH,
)
from pipeline_engine import read_int_setting, ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime
from pipelines.vessel_inputs import load_vessel_topology_inputs
from utils.logger import Logger

from .gap_masks import gaps_from_m0ff
from .outputs import (
    pack_capillary_analysis_zones,
    pack_gap_mask,
    pack_topology_outputs,
    pack_rotated_segment_masks,
)

TOPOLOGY_CORE_STATE = "topology_core.prepared"
CAPILLARY_ANALYSIS_ZONES_STATE = "topology_core.capillary_analysis_zones"
_DEFAULT_NUMBER_OF_RADII_IN_FOV = 25
_TOPOLOGY_WINDOW_SIZE_PERCENTILE_KEPT = 0.95


def _inputs(ctx):
    ctx.dag.handoff("inputs", InvocationResource(load_vessel_topology_inputs(ctx)))


def _gaps(ctx):
    inputs = ctx.dag.require("inputs")
    legacy, disc = _prepare_capillary_zones(inputs)
    zones, outputs = _prepare_gap_masks(ctx, inputs, legacy, disc)
    if legacy is not None:
        outputs.update(pack_capillary_analysis_zones(legacy, disc))
    accepted_zones, _ = ctx.dag.handoff("gaps", (zones, outputs))
    ctx.state.set(CAPILLARY_ANALYSIS_ZONES_STATE, accepted_zones)


def _sampling(ctx):
    inputs = ctx.dag.require("inputs")
    images, segmentation = inputs.image_maps, inputs.segmentation
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
    prepared = prepare_sampling_plans(
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
    accepted = ctx.dag.handoff("topology", prepared)
    ctx.state.set(TOPOLOGY_CORE_STATE, accepted)


def _pack_topology(prepared, ctx):
    inputs = ctx.dag.require("inputs")
    vessels = inputs.segmentation.vessels
    outputs = pack_topology_outputs(vessels.artery, vessels.vein, inputs.segmentation.optic_disc,
                                    prepared, inputs.holodoppler.pixel_pitch)
    if "segment_velocity_maps" in ctx.options_for("velocity_analysis"):
        outputs.update(pack_rotated_segment_masks(prepared))
    return outputs


def topology_graph(*, options=()):
    return StageGraph((
        StageSpec("inputs", _inputs, products=(ProductSpec("inputs"),)),
        StageSpec("gaps", _gaps, requires=("inputs",),
                  products=(ProductSpec("gaps", pack=lambda value, ctx: value[1]),)),
        StageSpec("sampling", _sampling, requires=("inputs",),
                  products=(ProductSpec("topology", pack=_pack_topology),)),
    ))


def run_topology_core(ctx):
    require_managed_runtime(ctx, "topology_core")
    return topology_graph(options=ctx.options_for("topology_core"))


def _prepare_gap_masks(ctx, inputs, legacy, optic_disc):
    schema = EyeFlowOutputPaths.active().segmentation
    volume = ctx.inputs.hd.as_holodoppler().moment0_flat_field_dataset()
    outputs = {}
    sources = []
    if legacy is not None:
        sources.append(f"/{schema.gaps_legacy}")
    if volume is None:
        Logger.log_warning(
            "HD moment0ff/M0FF is missing; gaps Mask uses Legacy when available."
        )
        combined = None if legacy is None else legacy.copy()
    else:
        expected_shape = inputs.segmentation.vessels.artery.shape
        if tuple(volume.shape[-2:]) != expected_shape:
            raise ValueError(
                f"moment0ff spatial shape {volume.shape[-2:]} does not match "
                f"the native HD frame {expected_shape}."
            )
        mask, attrs = gaps_from_m0ff(volume)
        attrs["source_path"] = volume.name
        outputs.update(pack_gap_mask(mask, schema.gaps_from_m0ff, attrs))
        sources.append(f"/{schema.gaps_from_m0ff}")
        # Persisted inputs are white-gap masks. OR their vessel/exclusion
        # complements first, then invert: a gap must be present in both inputs.
        combined = mask if legacy is None else ~(~legacy | ~mask)
    if combined is not None:
        combined = combined & capillary_disc_support(optic_disc.mask_for(combined.shape))
        disc_attrs = {
            "cutoff_disc_center": "image_center",
            "cutoff_disc_diameter_width_fraction": GAPS_CUTOFF_DIAMETER_WIDTH_FRACTION,
            "cutoff_disc_diameter_pixels": GAPS_CUTOFF_DIAMETER_WIDTH_FRACTION * combined.shape[1],
            "optic_disc_mask_source": (
                "frame_center_1.91_mm_diameter_fallback"
                if optic_disc.is_fallback else "dopplerview_segmentation"
            ),
        }
        if optic_disc.is_fallback:
            disc_attrs["optic_disc_fallback_diameter_mm"] = DEFAULT_OPTIC_DISC_DIAMETER_MM
        outputs.update(pack_gap_mask(combined, schema.gaps, {
            "definition": (
                "NOT ((NOT Legacy) OR (NOT GapsFromM0ff)) "
                "AND cutoff_disc AND NOT optic_disc"
            ),
            "source_paths": sources,
            "combination_operation": "OR",
            "combination_domain": "vessel_exclusion_masks_before_gap_inversion",
            "equivalent_gap_operation": "AND",
            "legacy_fallback": volume is None,
            **disc_attrs,
        }))
    return combined, outputs


def _prepare_capillary_zones(inputs) -> tuple[np.ndarray | None, OpticDisc]:
    segmentation = inputs.segmentation
    vessels = segmentation.vessels
    missing = [
        f"/{path}"
        for path, mask in (
            (DV_RETINAL_VESSEL_LOGITS_PATH, vessels.retinal_vessel_mask),
            (DV_CHOROID_VESSEL_MASK_PATH, vessels.choroid_vessel_mask),
        )
        if mask is None
    ]
    optic_disc = segmentation.optic_disc
    image_shape = vessels.artery.shape
    # Both legacy support and the final OR use exactly the same disc geometry,
    # including when the DV vessel masks required by Legacy are absent.
    if optic_disc.mask is None:
        optic_disc = OpticDisc.centered_fallback(
            image_shape,
            pixel_pitch_mm=inputs.holodoppler.pixel_pitch.isotropic_mm,
        )
    if missing:
        Logger.log_warning(
            "Omitting legacy capillary analysis zones: missing DopplerView masks "
            + ", ".join(missing)
            + ". Existing vessel analyses remain available."
        )
        return None, optic_disc
    return capillary_analysis_zones(
        vessels.retinal_vessel_mask,
        optic_disc.mask_for(image_shape),
        vessels.choroid_vessel_mask,
    ), optic_disc


def optional_capillary_zones(ctx) -> np.ndarray | None:
    """Return the combined native gap mask, or None when both sources are absent."""
    if ctx.state.get(CAPILLARY_ANALYSIS_ZONES_STATE) is None:
        return None
    return capillary_zones(ctx)


def capillary_zones(ctx) -> np.ndarray:
    """Return native, read-only support produced by ``capillary_analysis_zones``."""
    value = ctx.state.get(CAPILLARY_ANALYSIS_ZONES_STATE)
    if value is None:
        raise RuntimeError(
            "Capillary analysis zones are unavailable; check topology_core "
            "and the HD moment0ff or retinal/choroid DopplerView mask inputs."
        )
    if not isinstance(value, np.ndarray) or value.ndim != 2 or value.dtype != np.bool_:
        raise TypeError("Capillary analysis zones state must be a 2-D boolean array.")
    return value


def prepared_topologies(ctx) -> Mapping[str, SegmentSamplingPlan]:
    """Return the canonical topology produced by the declared DAG dependency."""

    value = ctx.state.get(TOPOLOGY_CORE_STATE)
    if not isinstance(value, Mapping) or set(value) != {"artery", "vein"}:
        raise RuntimeError("Prepared topology state is unavailable; check the pipeline DAG.")
    if not all(isinstance(item, SegmentSamplingPlan) for item in value.values()):
        raise TypeError("Prepared topology state contains an invalid value.")
    return value


__all__ = [
    "CAPILLARY_ANALYSIS_ZONES_STATE",
    "TOPOLOGY_CORE_STATE",
    "capillary_zones",
    "optional_capillary_zones",
    "prepared_topologies",
    "run_topology_core",
]


TOPOLOGY_CACHE_STATE = "topology.prepared"


def run_topology_cache(
    run_state: MutableMapping[str, object],
) -> dict[TopologyCacheKey, object]:
    """Return the sampling-plan dictionary owned by one pipeline run."""

    found = run_state.get(TOPOLOGY_CACHE_STATE)
    if found is None:
        cache: dict[TopologyCacheKey, object] = {}
        run_state[TOPOLOGY_CACHE_STATE] = cache
        return cache
    if not isinstance(found, dict):
        raise TypeError(
            f"Run state '{TOPOLOGY_CACHE_STATE}' must contain a dictionary."
        )
    return found
