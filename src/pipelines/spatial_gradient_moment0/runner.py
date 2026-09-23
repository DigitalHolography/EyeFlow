"""Run optional segment-level spatial-gradient and lumen measurements."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.waveform import (
    mean_period_seconds,
)
from pipelines.waveform_velocity_core.runner import WAVEFORM_CONTEXT_STATE

from .lumen_size import export_lumen_size_pngs
from .profiles import (
    SPATIAL_GRADIENT_METRICS_ROOT,
    extract_spatial_gradient_segments,
    pack_blood_volume_rate_outputs,
    pack_bvr_velocity_profile_outputs,
    pack_spatial_gradient_profile_outputs,
)


SPATIAL_GRADIENT_PRODUCTS_STATE = "spatial_gradient_moment0.products"


@dataclass(frozen=True)
class SpatialGradientProducts:
    """In-run gradient products shared without coupling their calculation."""

    artery_segments: object
    vein_segments: object
    outputs: dict[str, object]


def run_spatial_gradient_moment0(ctx) -> dict[str, object]:
    """Calculate staged gradient profiles only when this pipeline is selected."""

    context = ctx.state.get(WAVEFORM_CONTEXT_STATE)
    if context is None:
        raise RuntimeError(
            "Spatial-gradient profiles require the waveform velocity core state."
        )
    if context.artery_segment_result is None or context.vein_segment_result is None:
        raise RuntimeError(
            "Spatial-gradient profiles require prepared artery and vein topology."
        )

    cycle_boundaries = context.per_beat_analysis.cycle_boundary_indexes
    index_base = int(context.source_data.provenance["beat_index_base"])
    artery_segments, vein_segments = extract_spatial_gradient_segments(ctx, context)
    gradient_outputs = pack_spatial_gradient_profile_outputs(
        artery_segments,
        vein_segments,
        cycle_boundaries,
        index_base=index_base,
    )

    # BVR is downstream of gradient-derived lumen edges. Build the matching
    # velocity profiles from the already-retained core segment results.
    velocity_profile_outputs = pack_bvr_velocity_profile_outputs(
        context.artery_segment_result,
        context.vein_segment_result,
        cycle_boundaries,
        index_base=index_base,
    )
    _validate_profile_segment_alignment(
        "Artery",
        context.artery_segment_result,
        artery_segments,
    )
    _validate_profile_segment_alignment(
        "Vein",
        context.vein_segment_result,
        vein_segments,
    )
    outputs = dict(gradient_outputs)
    outputs.update(
        pack_blood_volume_rate_outputs(
            velocity_profile_outputs,
            gradient_outputs,
        )
    )
    # Keep BVR source paths self-contained when waveform_velocity was not
    # selected. If it was selected, that pipeline publishes the same profiles.
    if not ctx.pipeline_scheduled("waveform_velocity"):
        outputs.update(velocity_profile_outputs)

    output = getattr(ctx, "output", None)
    if getattr(output, "available", False):
        period_seconds = mean_period_seconds(
            cycle_boundaries,
            float(context.source_data.timing.dt_seconds),
        )
        for vessel_name, segments in (
            ("Artery", artery_segments),
            ("Vein", vein_segments),
        ):
            lumen_path = (
                f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel_name}/"
                "Transverse/Masked/tk/lumen_size"
            )
            lumen_size = gradient_outputs.get(lumen_path)
            if lumen_size is not None:
                export_lumen_size_pngs(
                    output,
                    lumen_size.data,
                    segments.branch_ids,
                    vessel_name=vessel_name,
                    period_seconds=period_seconds,
                )

    ctx.state.set(
        SPATIAL_GRADIENT_PRODUCTS_STATE,
        SpatialGradientProducts(
            artery_segments=artery_segments,
            vein_segments=vein_segments,
            outputs=outputs,
        ),
    )
    return outputs


def _validate_profile_segment_alignment(
    vessel_name: str,
    velocity_segments,
    gradient_segments,
) -> None:
    """Ensure velocity and gradient profiles share branch/radius identities."""

    if velocity_segments is None or gradient_segments is None:
        raise RuntimeError(
            f"{vessel_name} velocity and gradient segments are required for alignment."
        )
    for field in ("labels", "branch_ids"):
        if not np.array_equal(
            np.asarray(getattr(velocity_segments, field)),
            np.asarray(getattr(gradient_segments, field)),
        ):
            raise RuntimeError(
                f"{vessel_name} velocity and gradient segment {field} do not match."
            )
    velocity_centers = _ring_branch_segment_centers(velocity_segments)
    gradient_centers = _ring_branch_segment_centers(gradient_segments)
    if velocity_centers.shape != gradient_centers.shape or not np.allclose(
        velocity_centers,
        gradient_centers,
        equal_nan=True,
    ):
        raise RuntimeError(
            f"{vessel_name} velocity and gradient segment centers do not match."
        )
    velocity_shape = tuple(_unmasked_transverse_profiles(velocity_segments).shape[:2])
    gradient_shape = tuple(_unmasked_transverse_profiles(gradient_segments).shape[:2])
    if velocity_shape != gradient_shape:
        raise RuntimeError(
            f"{vessel_name} velocity and gradient (radius, branch) dimensions "
            f"do not match: {velocity_shape} != {gradient_shape}."
        )


def _ring_branch_segment_centers(segments) -> np.ndarray:
    if hasattr(segments, "segment_centers_xy"):
        return np.asarray(segments.segment_centers_xy)
    return np.transpose(np.asarray(segments.segment_center_xy), (1, 0, 2))


def _unmasked_transverse_profiles(segments) -> np.ndarray:
    if hasattr(segments, "transverse_profiles_unmasked"):
        return np.asarray(segments.transverse_profiles_unmasked)
    return np.asarray(segments.velocity_profiles)


__all__ = [
    "SPATIAL_GRADIENT_PRODUCTS_STATE",
    "SpatialGradientProducts",
    "_validate_profile_segment_alignment",
    "run_spatial_gradient_moment0",
]
