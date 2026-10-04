"""Shared physical calibration contract for retinal velocity estimation."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

DEFAULT_LASER_WAVELENGTH_METERS = 8.52e-7
DEFAULT_NUMERICAL_APERTURE = 0.124
DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ = 1.0
BAND_RATIO_CALIBRATION_MODEL = "linear_origin"
BAND_RATIO_CALIBRATION_SOURCE = "eyeflow_setting"
BAND_RATIO_CALIBRATION_VERSION = "1"
BAND_LF_LOW_RELATIVE_THRESHOLD = 1e-6


def validate_band_ratio_frequency_scale_hz(value: object) -> float:
    """Return a finite, positive Hz-per-ratio calibration factor."""

    try:
        scale_hz = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "band_ratio_frequency_scale_hz must be a finite positive value "
            "in Hz per ratio unit."
        ) from exc
    if not np.isfinite(scale_hz) or scale_hz <= 0.0:
        raise ValueError(
            "band_ratio_frequency_scale_hz must be a finite positive value "
            "in Hz per ratio unit."
        )
    return scale_hz


def physical_velocity_provenance(
    *,
    velocity_estimation_method: str,
    band_ratio_frequency_scale_hz: float,
    laser_wavelength_m: float = DEFAULT_LASER_WAVELENGTH_METERS,
    numerical_aperture: float = DEFAULT_NUMERICAL_APERTURE,
) -> dict[str, object]:
    """Return stable physical and method-specific calibration metadata."""

    provenance: dict[str, object] = {
        "velocity_estimation_method": str(velocity_estimation_method),
        "velocity_quantity": "physical_velocity",
        "velocity_unit": "mm/s",
        "laser_wavelength_m": float(laser_wavelength_m),
        "numerical_aperture": float(numerical_aperture),
    }
    if velocity_estimation_method == "frequency_bands":
        provenance.update(
            {
                "band_ratio_frequency_scale_hz": (
                    validate_band_ratio_frequency_scale_hz(
                        band_ratio_frequency_scale_hz
                    )
                ),
                "band_ratio_calibration_model": BAND_RATIO_CALIBRATION_MODEL,
                "band_ratio_calibration_source": BAND_RATIO_CALIBRATION_SOURCE,
                "band_ratio_calibration_version": BAND_RATIO_CALIBRATION_VERSION,
            }
        )
    return provenance


def calibration_attrs_from_metadata(
    metadata: Mapping[str, object] | None,
) -> dict[str, object]:
    """Extract calibration attributes suitable for velocity datasets."""

    if not isinstance(metadata, Mapping):
        return {}
    keys = (
        "band_ratio_frequency_scale_hz",
        "band_ratio_calibration_model",
        "band_ratio_calibration_source",
        "band_ratio_calibration_version",
        "laser_wavelength_m",
        "numerical_aperture",
    )
    return {key: metadata[key] for key in keys if metadata.get(key) is not None}


__all__ = [
    "BAND_LF_LOW_RELATIVE_THRESHOLD",
    "BAND_RATIO_CALIBRATION_MODEL",
    "BAND_RATIO_CALIBRATION_SOURCE",
    "BAND_RATIO_CALIBRATION_VERSION",
    "DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ",
    "DEFAULT_LASER_WAVELENGTH_METERS",
    "DEFAULT_NUMERICAL_APERTURE",
    "calibration_attrs_from_metadata",
    "physical_velocity_provenance",
    "validate_band_ratio_frequency_scale_hz",
]
