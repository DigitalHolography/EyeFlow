"""Shared units and display labels for velocity-estimation outputs."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class VelocitySemantics:
    """Physical meaning attached to values emitted by the velocity pipeline."""

    method: str
    quantity: str
    unit: str
    label: str

    @property
    def axis_label(self) -> str:
        if self.quantity == "relative_velocity_index":
            return self.label
        return f"Velocity ({self.unit})"

    @property
    def value_suffix(self) -> str:
        return "" if self.unit == "1" else f" {self.unit}"

    def dataset_attrs(self) -> dict[str, str]:
        return {
            "unit": self.unit,
            "velocity_estimation_method": self.method,
            "velocity_quantity": self.quantity,
        }


_PHYSICAL = VelocitySemantics(
    method="doppler_moments",
    quantity="physical_velocity",
    unit="mm/s",
    label="Velocity",
)
_RELATIVE = VelocitySemantics(
    method="frequency_bands",
    quantity="relative_velocity_index",
    unit="1",
    label="Relative velocity index",
)


def resolve_velocity_semantics(
    metadata: Mapping[str, object] | None = None,
    *,
    unit: str | None = None,
) -> VelocitySemantics:
    """Resolve velocity meaning, defaulting legacy callers to physical velocity.

    Estimator results expose the three top-level provenance keys used here.  A
    nested ``provenance`` mapping is also accepted for imported/external
    analyses.  The method or quantity takes precedence over a stale unit so a
    relative result can never be presented as calibrated ``mm/s``.
    """

    provenance = _metadata_mapping(metadata, "provenance")
    method = _metadata_value(metadata, provenance, "velocity_estimation_method")
    quantity = _metadata_value(metadata, provenance, "velocity_quantity")
    resolved_unit = unit or _metadata_value(metadata, provenance, "velocity_unit")

    if method == "frequency_bands" or quantity == "relative_velocity_index":
        return _RELATIVE
    if method == "doppler_moments" or quantity == "physical_velocity":
        return _PHYSICAL
    if resolved_unit == "1":
        return _RELATIVE
    return _PHYSICAL


def velocity_unit_from_payload(value: object, *, default: str = "mm/s") -> str:
    """Read an HDF-style unit attribute from a packed metric payload."""

    attrs = getattr(value, "attrs", None)
    if attrs is None and isinstance(value, tuple) and len(value) == 2:
        attrs = value[1] if isinstance(value[1], Mapping) else None
    if isinstance(attrs, Mapping):
        raw_unit = attrs.get("unit")
        if raw_unit is not None:
            return str(raw_unit)
    return default


def _metadata_mapping(
    metadata: Mapping[str, object] | None,
    key: str,
) -> Mapping[str, object] | None:
    if not isinstance(metadata, Mapping):
        return None
    value = metadata.get(key)
    return value if isinstance(value, Mapping) else None


def _metadata_value(
    metadata: Mapping[str, object] | None,
    provenance: Mapping[str, object] | None,
    key: str,
) -> str | None:
    for source in (metadata, provenance):
        if isinstance(source, Mapping):
            value = source.get(key)
            if value is not None:
                return str(value)
    return None


__all__ = [
    "VelocitySemantics",
    "resolve_velocity_semantics",
    "velocity_unit_from_payload",
]
