"""Pack weighted parabolic velocity-profile fits into EyeFlow outputs."""

from __future__ import annotations

from collections.abc import Mapping
from time import perf_counter

from calculations.vessel_segments.profiles.fits import (
    COUNT_OUTPUTS,
    DEFAULT_TIME_BLOCK_SIZE,
    DEFAULT_WEIGHT_POWER,
    FLOAT_OUTPUTS,
    border_weights,
    fit_quadratic_profiles,
)
from input_output import EyeFlowOutputPaths
from pipeline_engine.base import DatasetValue
from utils.logger import Logger

_OUTPUT_PATHS = EyeFlowOutputPaths.active()
SOURCE_PATHS = {
    "Artery": (
        "/" + _OUTPUT_PATHS.artery_velocity_profiles.transverse_velocity_profile_masked
    ),
    "Vein": (
        "/" + _OUTPUT_PATHS.vein_velocity_profiles.transverse_velocity_profile_masked
    ),
}
OUTPUT_ROOT = "/" + _OUTPUT_PATHS.velocity_profile_analysis_root
_MISSING = object()


def pack_velocity_profile_analysis_outputs(
    profile_outputs: Mapping[str, object],
) -> dict[str, object]:
    """Fit both vessel classes from waveform profile output payloads."""

    datasets = {
        vessel: _source_values(profile_outputs, vessel, source_path)
        for vessel, source_path in SOURCE_PATHS.items()
    }

    outputs: dict[str, object] = {}
    for vessel, values in datasets.items():
        started = perf_counter()
        results = fit_quadratic_profiles(values)
        Logger.log(
            f"Completed {vessel.lower()} weighted velocity-profile analysis in "
            f"{perf_counter() - started:.1f}s."
        )
        attrs = {
            "dimDesc": ["time", "beat", "branch", "radius"],
            "source_path": SOURCE_PATHS[vessel],
            "index_base": 0,
            "model": "a*x^2 + b*x + c",
            "fit_method": "weighted_least_squares",
            "weight_definition": "u=x/(Nx-1); d=abs(2*u-1); w=1-d^p",
            "weight_power": DEFAULT_WEIGHT_POWER,
            "integration_method": (
                "unweighted_sum_of_finite_observed_integer_indexes_between_roots"
            ),
            "geometry_policy": "downward_opening_only; fractional_indexes",
            "vessel": vessel.lower(),
        }
        for name, value in results.items():
            outputs[f"{OUTPUT_ROOT}/{vessel}/{name}/value"] = DatasetValue(
                data=value,
                attrs=dict(attrs),
            )
    return outputs


def _source_values(
    profile_outputs: Mapping[str, object],
    vessel: str,
    source_path: str,
):
    payload = _MISSING
    for candidate in (source_path, source_path.lstrip("/")):
        if candidate in profile_outputs:
            payload = profile_outputs[candidate]
            break
    if payload is _MISSING:
        raise KeyError(
            f"Required {vessel.lower()} velocity-profile output is missing: "
            f"{source_path}. velocity_analysis must pack both vessel profiles "
            "before this analysis."
        )
    if isinstance(payload, DatasetValue):
        return payload.data
    if isinstance(payload, tuple) and payload:
        return payload[0]
    return payload


# Preserve the former pipeline-local calculation name for downstream imports.
analyze_velocity_profiles = fit_quadratic_profiles


__all__ = [
    "COUNT_OUTPUTS",
    "DEFAULT_TIME_BLOCK_SIZE",
    "DEFAULT_WEIGHT_POWER",
    "FLOAT_OUTPUTS",
    "OUTPUT_ROOT",
    "SOURCE_PATHS",
    "analyze_velocity_profiles",
    "border_weights",
    "pack_velocity_profile_analysis_outputs",
]
