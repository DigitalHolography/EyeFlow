"""Resolve the minimal inputs required to detect heartbeat boundaries."""

from __future__ import annotations

from dataclasses import dataclass

from calculations.topology import OpticDisc
from input_output.schema import HolodopplerTiming
from pipelines.vessel_inputs import (
    align_mask as _align_mask,
    apply_mask_alignment as _apply_mask_alignment,
    load_vessel_topology_inputs,
)


@dataclass(frozen=True)
class HeartbeatInputs:
    moment0: object
    moment2: object
    artery_mask: np.ndarray
    vein_mask: np.ndarray
    optic_disc: OpticDisc
    timing: HolodopplerTiming
    local_background_dist: int
    index_base: int = 0


def load_heartbeat_inputs(ctx) -> HeartbeatInputs:
    """Load and spatially align only the arrays used for beat detection."""

    topology = load_vessel_topology_inputs(ctx)
    hd = ctx.inputs.hd.as_holodoppler()
    dv = ctx.inputs.dv.as_dopplerview()
    moment2 = hd.moment2_dataset()
    return HeartbeatInputs(
        moment0=topology.moment0,
        moment2=moment2,
        artery_mask=topology.artery_mask,
        vein_mask=topology.vein_mask,
        optic_disc=topology.optic_disc,
        timing=hd.timing(),
        local_background_dist=dv.local_background_dist(),
    )


__all__ = ["HeartbeatInputs", "load_heartbeat_inputs"]
