"""Shared units and display labels for velocity-estimation outputs."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from velocity_calibration import calibration_attrs_from_metadata


@dataclass(frozen=True)
class VelocitySemantics:
    """Physical meaning attached to values emitted by the velocity pipeline."""

    method: str
    quantity: str
    unit: str
    label: str

    @property
    def axis_label(self) -> str:
        return f"Velocity ({self.unit})"

    @property
    def value_suffix(self) -> str:
        return f" {self.unit}"

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
_FREQUENCY_BANDS_PHYSICAL = VelocitySemantics(
    method="frequency_bands",
    quantity="physical_velocity",
    unit="mm/s",
    label="Velocity",
)


def resolve_velocity_semantics(
    metadata: Mapping[str, object] | None = None,
    *,
    unit: str | None = None,
) -> VelocitySemantics:
    """Resolve velocity meaning and reject non-physical velocity metadata.

    Estimator results expose the three top-level provenance keys used here.  A
    nested ``provenance`` mapping is also accepted for imported/external
    analyses. EyeFlow velocity is always physical and expressed in ``mm/s``.
    """

    provenance = _metadata_mapping(metadata, "provenance")
    method = _metadata_value(metadata, provenance, "velocity_estimation_method")
    quantity = _metadata_value(metadata, provenance, "velocity_quantity")
    resolved_unit = unit or _metadata_value(metadata, provenance, "velocity_unit")

    if quantity not in {None, "physical_velocity"}:
        raise ValueError(
            "EyeFlow velocity must have velocity_quantity='physical_velocity'; "
            f"got {quantity!r}."
        )
    if resolved_unit not in {None, "mm/s"}:
        raise ValueError(
            "EyeFlow velocity must use unit 'mm/s'; "
            f"got {resolved_unit!r}."
        )
    if method == "frequency_bands":
        return _FREQUENCY_BANDS_PHYSICAL
    if method in {None, "doppler_moments"}:
        return _PHYSICAL
    raise ValueError(f"Unsupported velocity_estimation_method {method!r}.")


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


def velocity_dataset_attrs(
    metadata: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Return physical velocity and calibration provenance for one dataset."""

    semantics = resolve_velocity_semantics(metadata)
    attrs: dict[str, object] = semantics.dataset_attrs()
    provenance = _metadata_mapping(metadata, "provenance")
    attrs.update(calibration_attrs_from_metadata(provenance))
    attrs.update(calibration_attrs_from_metadata(metadata))
    return attrs


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
    "velocity_dataset_attrs",
    "velocity_unit_from_payload",
]
