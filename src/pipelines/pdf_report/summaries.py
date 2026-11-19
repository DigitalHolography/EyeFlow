"""Presentation-neutral summary rows for moments reports."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from .data import ReportData, VesselGeometryReportData, VolumeRateReportData


@dataclass(frozen=True)
class SummaryValue:
    value: float | int | None
    unit: str = ""
    decimals: int = 2
    regular_unit: bool = False


@dataclass(frozen=True)
class SummaryRow:
    label: str
    value: SummaryValue


def heart_rate_summary(parameters: Mapping[str, Any]) -> SummaryValue:
    return _parameter_summary(parameters, "heart_beat", decimals=1)


def velocity_summary_rows(
    parameters: Mapping[str, Any], vessel: str,
) -> tuple[SummaryRow, ...]:
    arterial = vessel == "artery"
    prefix = "Arterial" if arterial else "Venous"
    metric_prefix = "A" if arterial else "V"
    specs = (
        ("Minimum filtered velocity", f"Min_{prefix}_Velocity"),
        ("Mean filtered velocity", f"Average_{prefix}_Velocity"),
        ("Maximum filtered velocity", f"Max_{prefix}_Velocity"),
        ("Resistivity index", f"{metric_prefix}RI"),
        ("Pulsatility index", f"{metric_prefix}PI"),
    )
    return tuple(
        SummaryRow(label, _parameter_summary(parameters, key))
        for label, key in specs
    )


def analysis_summary_rows(
    report: ReportData, vessel: str,
) -> tuple[SummaryRow, ...]:
    geometry = report.vessel_geometry.get(vessel)
    rate = report.mask_volume_rates.get(vessel)
    rate_values = _volume_rate_values(rate)
    return (
        SummaryRow(
            "Sampled branches",
            _geometry_count(geometry, "sampled_branches"),
        ),
        SummaryRow(
            "Sampling regions",
            _geometry_count(geometry, "valid_sampling_regions"),
        ),
        SummaryRow(
            "Mean lumen diameter",
            _geometry_value(geometry, "mean_lumen_diameter_microns"),
        ),
        SummaryRow(
            "Lumen diameter SD",
            _geometry_value(
                geometry, "lumen_diameter_standard_deviation_microns"
            ),
        ),
        SummaryRow("Median BVR", _rate_value(rate_values["median"], rate)),
        SummaryRow("Minimum BVR", _rate_value(rate_values["minimum"], rate)),
        SummaryRow("Maximum BVR", _rate_value(rate_values["maximum"], rate)),
        SummaryRow("Pulsatile range", _rate_value(rate_values["range"], rate)),
    )


def format_plain_value(value: SummaryValue) -> str:
    """Format one value for text/HTML without Matplotlib math markup."""
    if value.value is None:
        return "Unavailable"
    number = f"{value.value:.{value.decimals}f}"
    unit = _plain_unit(value.unit)
    return f"{number}{f' {unit}' if unit else ''}"


def _parameter_summary(
    parameters: Mapping[str, Any], name: str, *, decimals: int = 2,
) -> SummaryValue:
    parameter = parameters.get(name)
    if not isinstance(parameter, Mapping):
        return SummaryValue(None, decimals=decimals)
    try:
        value = float(parameter["value"])
    except (KeyError, TypeError, ValueError):
        return SummaryValue(None, decimals=decimals)
    if not np.isfinite(value):
        return SummaryValue(None, decimals=decimals)
    return SummaryValue(
        value,
        unit=str(parameter.get("unit", "")).strip(),
        decimals=decimals,
    )


def _volume_rate_values(
    data: VolumeRateReportData | None,
) -> dict[str, float | None]:
    finite = (
        np.asarray([], dtype=float)
        if data is None
        else data.phase_curve[np.isfinite(data.phase_curve)]
    )
    minimum = float(np.min(finite)) if finite.size else None
    maximum = float(np.max(finite)) if finite.size else None
    return {
        "median": None if data is None else data.overall_value,
        "minimum": minimum,
        "maximum": maximum,
        "range": None if not finite.size else maximum - minimum,
    }


def _geometry_count(
    data: VesselGeometryReportData | None, field_name: str,
) -> SummaryValue:
    value = None if data is None else int(getattr(data, field_name))
    return SummaryValue(value, decimals=0)


def _geometry_value(
    data: VesselGeometryReportData | None, field_name: str,
) -> SummaryValue:
    value = None if data is None else getattr(data, field_name)
    return SummaryValue(value, unit="µm", regular_unit=True)


def _rate_value(
    value: float | None, data: VolumeRateReportData | None,
) -> SummaryValue:
    return SummaryValue(value, unit="mm^3/s" if data is None else data.unit)


def _plain_unit(unit: str) -> str:
    normalized = str(unit).strip().replace(" ", "")
    return {
        "mm/s": "mm.s⁻¹",
        "mm.s^-1": "mm.s⁻¹",
        "mm^3/s": "mm³.s⁻¹",
        "mm3/s": "mm³.s⁻¹",
        "µm": "µm",
    }.get(normalized, str(unit).strip())


__all__ = [
    "SummaryRow",
    "SummaryValue",
    "analysis_summary_rows",
    "format_plain_value",
    "heart_rate_summary",
    "velocity_summary_rows",
]
