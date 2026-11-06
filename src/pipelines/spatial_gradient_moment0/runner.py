"""Run optional segment-level spatial-gradient and lumen measurements."""

from __future__ import annotations

from dataclasses import dataclass

from calculations.blood_flow_velocity.signal_analysis.waveform import (
    mean_period_seconds,
)
from calculations.vessel_segments.measurement import SegmentMeasurementSettings
from data_manager import InvocationResource
from pipeline_engine import ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime
from pipelines.topology_core.runner import prepared_topologies
from pipelines.velocity.runner import cardiac_cycle_indexes
from pipelines.vessel_inputs import load_vessel_topology_inputs

from .lumen_size import export_lumen_size_pngs
from .profiles import (
    SPATIAL_GRADIENT_METRICS_ROOT,
    extract_spatial_gradient_segments,
    pack_spatial_gradient_profile_outputs,
)

SPATIAL_GRADIENT_PRODUCTS_STATE = "spatial_gradient_moment0.products"


@dataclass(frozen=True)
class SpatialGradientProducts:
    """In-run gradient products shared without coupling their calculation."""

    artery_segments: object
    vein_segments: object
    outputs: dict[str, object]


def _inputs(ctx):
    inputs = load_vessel_topology_inputs(ctx)
    ctx.dag.handoff("inputs", (InvocationResource(inputs), cardiac_cycle_indexes(ctx)))


def _profiles(ctx):
    inputs, boundaries = ctx.dag.require("inputs")
    artery, vein = extract_spatial_gradient_segments(ctx, inputs, prepared_topologies(ctx),
        profile_settings=SegmentMeasurementSettings(pixel_size_mm=inputs.holodoppler.pixel_pitch.isotropic_mm))
    outputs = pack_spatial_gradient_profile_outputs(artery, vein, boundaries, index_base=0)
    products = SpatialGradientProducts(artery, vein, outputs)
    accepted = ctx.dag.handoff("profiles", products)
    ctx.state.set(SPATIAL_GRADIENT_PRODUCTS_STATE, accepted)


def _pack_profiles(products, ctx):
    return products.outputs if ctx.pipeline_targeted("spatial_gradient_moment0") else {}


def _lumen_figures(ctx):
    if not ctx.output.available:
        return
    inputs, boundaries = ctx.dag.require("inputs")
    products = ctx.dag.require("profiles")
    period = mean_period_seconds(boundaries, float(inputs.holodoppler.timing.dt_seconds))
    for vessel, segments in (("Artery", products.artery_segments), ("Vein", products.vein_segments)):
        path = f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel}/Transverse/Masked/tk/lumen_size"
        payload = products.outputs.get(path)
        if payload is not None:
            export_lumen_size_pngs(ctx.output, payload.data, segments.topology.native.branch_ids,
                                  vessel_name=vessel, period_seconds=period)


def spatial_gradient_graph(*, options=()):
    return StageGraph((
        StageSpec("inputs", _inputs, products=(ProductSpec("inputs"),)),
        StageSpec("profiles", _profiles, requires=("inputs",),
                  products=(ProductSpec("profiles", pack=_pack_profiles),)),
        StageSpec("lumen_figures", _lumen_figures, requires=("inputs", "profiles"), artifact=True),
    ))


def run_spatial_gradient_moment0(ctx):
    require_managed_runtime(ctx, "spatial_gradient_moment0")
    return spatial_gradient_graph(options=ctx.options_for("spatial_gradient_moment0"))


__all__ = [
    "SPATIAL_GRADIENT_PRODUCTS_STATE",
    "SpatialGradientProducts",
    "run_spatial_gradient_moment0",
]
