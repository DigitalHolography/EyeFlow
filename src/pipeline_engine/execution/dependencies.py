"""Availability of whole-invocation prerequisites for dispatch and input planning.

Product requirements use accepted handles independently: a partial invocation can
still provide a committed product without satisfying a whole-invocation edge.
"""

from enum import Enum


class DependencyAvailability(Enum):
    WAITING = "waiting"
    SATISFIED = "satisfied"
    UNAVAILABLE = "unavailable"


def dependency_availability(status: str | None) -> DependencyAvailability:
    if status == "succeeded":
        return DependencyAvailability.SATISFIED
    if status in (None, "pending", "running"):
        return DependencyAvailability.WAITING
    # Partial, failed and skipped invocations cannot satisfy whole-invocation
    # dependencies. Unknown terminal outcomes must not authorize speculative I/O.
    return DependencyAvailability.UNAVAILABLE
