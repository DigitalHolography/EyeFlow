"""Select one shared cardiac timing source for run-once displacement outputs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pipelines.velocity import velocity_workflows
from pipelines.velocity.models import RetinalVelocity


@dataclass(frozen=True)
class DisplacementTiming:
    velocity: RetinalVelocity
    workflow: str

    @property
    def boundaries(self) -> np.ndarray:
        return np.asarray(self.velocity.cycle_boundary_indexes, dtype=np.int64)

    @property
    def attrs(self) -> dict[str, object]:
        return {
            "cycle_boundary_workflow": self.workflow,
            "cycle_boundary_vessel_source": self.velocity.cardiac_cycle_source,
            "cycle_boundary_indices": self.boundaries,
            "cycle_boundary_index_base": np.int32(0),
            "cycle_duration_seconds": self.velocity.cycle_durations_seconds,
            "cycle_interpolation": "periodic_fourier_resampling",
        }


def select_displacement_timing(ctx) -> DisplacementTiming | None:
    """Prefer moments timing, then fall back to calibrated band-ratio timing."""

    workflows = velocity_workflows(ctx)
    for workflow, candidate in (
        ("doppler_moments", workflows.moments),
        ("frequency_bands", workflows.bandratio),
    ):
        if (
            candidate is not None
            and candidate.cardiac_cycle_source != "none"
            and candidate.cycle_boundary_indexes.size >= 2
        ):
            return DisplacementTiming(candidate, workflow)
    return None


__all__ = ["DisplacementTiming", "select_displacement_timing"]
