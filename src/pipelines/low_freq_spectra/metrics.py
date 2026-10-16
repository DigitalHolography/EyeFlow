"""Peak, valley, and inclusive-area metrics for gap spectra."""

from __future__ import annotations

from numbers import Integral

import numpy as np

from calculations.math import nanmedian
from input_output.schema import EyeFlowOutputPaths
from utils.logger import Logger

from .per_beat import early_late_phase_windows

DEFAULT_VALLEY_NEIGHBORHOOD_FRACTION = 1 / 4
PHASE_METRIC_NAMES = ("A0", "A1", "A2", "R1", "R2")

SPECTRUM_METRIC_DEFINITIONS = {
    "index_max1": "zero-based index of the maximum in the first half of the bin domain",
    "max1": "maximum spectrum value in the first half of the bin domain, ignoring NaNs",
    "index_max2": "zero-based index of the maximum in the second half of the bin domain",
    "max2": "maximum spectrum value in the second half of the bin domain, ignoring NaNs",
    "index_min1": "zero-based central valley index floor(bin_count/2) - 1",
    "min1": "spectrum value at zero-based index floor(bin_count/2) - 1",
    "index_min2": "zero-based central valley index floor(bin_count/2)",
    "min2": "spectrum value at zero-based index floor(bin_count/2)",
    "A0": "arithmetic mean of max1 and max2 times the inclusive peak-to-peak bin count",
    "A1": "unweighted spectrum sum from index_max1 through index_max2, inclusive; NaNs propagate",
    "A2": "NaN-aware central-neighborhood median height times the inclusive peak-to-peak bin count",
    "R1": "A1 divided by A0; zero or missing A0 produces NaN",
    "R2": "A2 divided by A0; zero or missing A0 produces NaN",
}

PHASE_METRIC_DEFINITIONS = {
    "A0": "inclusive fixed peak-to-peak bin count times mean spectrum height at those two bins per phase",
    "A1": "inclusive spectrum sum between the fixed peak bins per phase; interval NaNs propagate",
    "A2": "inclusive fixed peak-to-peak bin count times NaN-aware central-neighborhood median per phase",
    "R1": "A1_t divided by A0_t per phase; zero or missing A0_t produces NaN",
    "R2": "A2_t divided by A0_t per phase; zero or missing A0_t produces NaN",
}


def analyze_spectrum(
    spectrum, *, valley_neighborhood_length: int | None = None,
) -> dict[str, float]:
    """Analyze a spectrum in its stored frequency order.

    Select one maximum in each half, including endpoints and ignoring NaNs.
    Ties select the lowest index. An all-NaN half has a NaN index and value.
    Valley indexes are the two central bins, even if their values are NaN.
    The peak-height mean requires both inputs. Areas include both peak bins without wrapping,
    interpolation, or frequency-spacing weights. A1 propagates interval NaNs.
    The even-length valley neighborhood is centered between the central bins;
    its default length is one quarter of the bin count, rounded to the nearest
    even number (ties upward), with at least two bins. Its median ignores NaNs, returning
    NaN if all samples are missing. A2 uses the same inclusive width as A0.
    Any one-dimensional domain with at least two bins is supported. Odd domains
    split at floor(bin_count/2). Unsupported shapes return NaN metrics with a warning.
    """
    values = np.asarray(spectrum)
    result = {name: float("nan") for name in SPECTRUM_METRIC_DEFINITIONS}
    if not _supported_bin_domain(values):
        Logger.log_warning(
            "Gap spectrum metrics require a one-dimensional bin domain with at least "
            f"2 bins; received shape {values.shape}. "
            "Returning NaN for all spectrum metrics."
        )
        return result
    values = np.asarray(values, dtype=np.float64)
    middle = values.size // 2
    valley_bins = _valley_neighborhood_indices(valley_neighborhood_length, values.size)
    result["index_max1"], result["max1"] = _peak(values, 0, middle)
    result["index_max2"], result["max2"] = _peak(values, middle, values.size)
    for name, start, stop in (("1", 0, middle), ("2", middle, values.size)):
        if np.isnan(result[f"index_max{name}"]):
            Logger.log_warning(
                f"Gap spectrum bins [{start}, {stop}) contain only NaN samples; "
                f"index_max{name} and max{name} remain NaN."
            )
    result.update(
        index_min1=float(middle - 1), min1=float(values[middle - 1]),
        index_min2=float(middle), min2=float(values[middle]),
    )
    h0 = (result["max1"] + result["max2"]) / 2
    hv = float(nanmedian(values[valley_bins], dtype=np.float64))
    if np.isfinite(result["index_max1"]) and np.isfinite(result["index_max2"]):
        start = int(result["index_max1"])
        stop = int(result["index_max2"]) + 1
        result["A0"] = h0 * (stop - start)
        result["A1"] = float(np.sum(values[start:stop]))
        result["A2"] = hv * (stop - start)
        if result["A0"] != 0:
            result["R1"] = result["A1"] / result["A0"]
            result["R2"] = result["A2"] / result["A0"]
    return result


def _supported_bin_domain(values: np.ndarray) -> bool:
    return values.ndim == 1 and values.size >= 2


def _valley_neighborhood_indices(length: int | None, bin_count: int) -> np.ndarray:
    if length is None:
        length = max(2, 2 * int(bin_count * DEFAULT_VALLEY_NEIGHBORHOOD_FRACTION / 2 + 0.5))
    if (
        isinstance(length, bool) or not isinstance(length, Integral)
        or length < 2 or length > bin_count or length % 2
    ):
        raise ValueError(
            f"valley_neighborhood_length must be an even integer from 2 through {bin_count}."
        )
    middle = bin_count // 2
    return np.arange(middle - length // 2, middle + length // 2)


def _peak(values: np.ndarray, start: int, stop: int) -> tuple[float, float]:
    half = values[start:stop]
    valid = np.flatnonzero(~np.isnan(half))
    if valid.size == 0:
        return float("nan"), float("nan")
    index = start + int(valid[np.argmax(half[valid])])
    return float(index), float(values[index])


def pack_spectrum_metrics(
    spectrum, *, source_attrs: dict[str, object],
    valley_neighborhood_length: int | None = None,
    phase_spectrum=None, phase_attrs: dict[str, object] | None = None,
) -> dict[str, object]:
    """Pack scalar metrics with units and references to the aggregate spectrum."""
    schema = EyeFlowOutputPaths.active()
    result = analyze_spectrum(spectrum, valley_neighborhood_length=valley_neighborhood_length)
    values = np.asarray(spectrum)
    valley_bins = (
        _valley_neighborhood_indices(valley_neighborhood_length, values.size)
        if _supported_bin_domain(values) else np.array([], dtype=int)
    )
    attrs = {
        "source_path": "/" + schema.spectrograms_masked_time_aggregated_median,
        "source_bin_count": int(np.size(spectrum)),
        "minimum_bin_count": 2,
        "bin_domain_split_index": values.size // 2 if _supported_bin_domain(values) else 0,
        "peak_half_valid_bin_counts": (
            [
                int(np.count_nonzero(~np.isnan(half)))
                for half in np.array_split(values, [values.size // 2])
            ]
            if _supported_bin_domain(values) else []
        ),
        "dimDesc": [],
        "peak_selection": "maximum_per_bin_domain_half_ignore_NaN_lowest_index_on_ties",
        "interval_endpoints": "both_included",
        "frequency_bin_order": source_attrs.get("frequency_bin_order", "HD_export_order"),
        "nonfinite_value_policy": "missing_required_values_and_interval_NaNs_propagate",
        "valley_neighborhood_length": len(valley_bins),
        "valley_neighborhood_indices": valley_bins,
        "valley_neighborhood_selection": (
            "quarter_of_bin_domain" if valley_neighborhood_length is None else "explicit_length"
        ),
        "valley_neighborhood_default_rounding": "nearest_even_ties_up_minimum_2",
        "valley_aggregation_method": "nanmedian",
        "valley_nonfinite_value_policy": "ignore_NaN_all_NaN_neighborhood_returns_NaN",
    }
    output = {}
    for name, value in result.items():
        unit = (
            "1"
            if name.startswith("index_") or name in {"R1", "R2"}
            else source_attrs.get("unit", "a.u.")
        )
        metric_attrs = {**attrs, "unit": unit, "definition": SPECTRUM_METRIC_DEFINITIONS[name]}
        if name.startswith("index_"):
            metric_attrs["index_base"] = 0
            metric_attrs["index_storage"] = "integer_valued_float32_with_NaN_for_missing_index"
        output[f"{schema.spectrograms_masked_metrics_root}/{name}"] = (
            np.float32(value),
            metric_attrs,
        )
    if phase_spectrum is not None:
        if phase_attrs is None:
            raise ValueError("Phase-spectrum attributes must provide phase_samples.")
        if np.ndim(phase_spectrum) != 2 or np.shape(phase_spectrum)[1] != values.size:
            raise ValueError("Phase and pooled spectra must share the same frequency-bin domain.")
        output.update(_pack_phase_metrics(
            phase_spectrum, fixed_metrics=result, source_attrs=phase_attrs,
            metric_attrs=attrs, valley_neighborhood_length=valley_neighborhood_length,
        ))
    return output


def analyze_phase_spectrum(
    spectrum, *, index_max1: float, index_max2: float,
    valley_neighborhood_length: int | None = None,
) -> dict[str, np.ndarray]:
    """Calculate phase curves with peak bins fixed from the pooled spectrum.

    Input axes are (phase, freq). A0 uses the mean height at those fixed bins;
    A1 is their inclusive interval sum; A2 uses the central neighborhood median.
    Missing peak indexes leave all curves NaN. At each phase, missing peak
    heights invalidate A0, interval NaNs invalidate A1, and A2 ignores valley
    NaNs unless the whole neighborhood is missing. Zero A0 yields NaN ratios.
    """
    values = np.asarray(spectrum, dtype=np.float64)
    if values.ndim != 2 or any(size == 0 for size in values.shape):
        raise ValueError("Phase spectra must have nonempty (phase, freq) axes.")
    result = {name: np.full(values.shape[0], np.nan) for name in PHASE_METRIC_NAMES}
    if values.shape[1] < 2 or np.isnan(index_max1) or np.isnan(index_max2):
        return result
    if not (
        np.isfinite(index_max1) and np.isfinite(index_max2)
        and index_max1 == int(index_max1) and index_max2 == int(index_max2)
        and 0 <= index_max1 < index_max2 < values.shape[1]
    ):
        raise ValueError("Fixed peak indices must be ordered integer bin indices within the spectrum.")
    start, stop = int(index_max1), int(index_max2) + 1
    width = stop - start
    valley_bins = _valley_neighborhood_indices(valley_neighborhood_length, values.shape[1])
    result["A0"] = width * (values[:, start] + values[:, stop - 1]) / 2
    result["A1"] = np.sum(values[:, start:stop], axis=1)
    result["A2"] = width * nanmedian(values[:, valley_bins], axis=1, dtype=np.float64)
    for name, numerator in (("R1", "A1"), ("R2", "A2")):
        np.divide(
            result[numerator], result["A0"], out=result[name],
            where=(result["A0"] != 0) & ~np.isnan(result["A0"]),
        )
    return result


def _pack_phase_metrics(
    spectrum, *, fixed_metrics: dict[str, float], source_attrs: dict[str, object],
    metric_attrs: dict[str, object], valley_neighborhood_length: int | None,
) -> dict[str, object]:
    schema = EyeFlowOutputPaths.active()
    root = schema.spectrograms_masked_metrics_root
    phase = np.asarray(source_attrs["phase_samples"])
    if phase.dtype.kind != "f":
        phase = phase.astype(np.float64)
    if phase.shape != (np.shape(spectrum)[0],):
        raise ValueError("Phase samples must be strictly increasing coordinates in [0, 1).")
    windows = early_late_phase_windows(phase)
    curves = analyze_phase_spectrum(
        spectrum, index_max1=fixed_metrics["index_max1"], index_max2=fixed_metrics["index_max2"],
        valley_neighborhood_length=valley_neighborhood_length,
    )
    attrs = {
        **{key: value for key, value in source_attrs.items() if not key.startswith("aggregation_")},
        **metric_attrs,
        "source_path": "/" + schema.spectrograms_masked_time_aggregated_median_beat,
        "peak_source_path": "/" + schema.spectrograms_masked_time_aggregated_median,
        "fixed_peak_index_source_paths": [f"/{root}/index_max1", f"/{root}/index_max2"],
        "fixed_peak_indices": [fixed_metrics["index_max1"], fixed_metrics["index_max2"]],
        "peak_selection": "fixed_pooled_spectrum_indices_at_all_phases",
        "calculation_method": "fixed_peak_bin_geometry_on_beat_median_spectrum",
        "nonfinite_value_policy": "missing_peak_heights_and_interval_NaNs_propagate_valley_ignores_NaN",
        "phase_samples": phase,
        "phase_endpoint": "excluded",
        "dimDesc": ["phase"],
    }
    output = {}
    for name, curve in curves.items():
        unit = "1" if name in {"R1", "R2"} else source_attrs.get("unit", "a.u.")
        curve_path = f"{root}/{name}_t"
        curve_attrs = {
            **attrs, "unit": unit,
            "definition": PHASE_METRIC_DEFINITIONS[name],
        }
        output[curve_path] = (curve.astype(np.float32), curve_attrs)
        scalar_attrs = {
            key: value for key, value in curve_attrs.items()
            if key not in {"phase_samples", "phase_endpoint"}
        }
        scalar_attrs.update(source_path="/" + curve_path, dimDesc=[])
        aggregates = {}
        for window, (mask, interval) in windows.items():
            aggregate = float(nanmedian(curve[mask], dtype=np.float64))
            aggregates[window] = aggregate
            output[f"{root}/{name}_{window}"] = (np.float32(aggregate), {
                **scalar_attrs,
                "definition": f"NaN-aware median of {name}_t over the {window} phase interval",
                "aggregation_method": "nanmedian",
                "aggregation_axes": ["phase"],
                "aggregation_weighting": "equal_weight_per_valid_phase_sample",
                "phase_interval": interval,
                "phase_interval_endpoints": "both_included",
                "source_phase_sample_count": int(np.count_nonzero(mask)),
                "valid_sample_count": int(np.count_nonzero(~np.isnan(curve[mask]))),
                "nonfinite_value_policy": "ignore_NaN_all_NaN_window_returns_NaN",
            })
        output[f"{root}/delta_{name}"] = (np.float32(aggregates["early"] - aggregates["late"]), {
            **scalar_attrs,
            "definition": f"{name}_early minus {name}_late",
            "source_paths": [f"/{root}/{name}_early", f"/{root}/{name}_late"],
            "calculation_method": "early_median_minus_late_median",
            "nonfinite_value_policy": "missing_early_or_late_median_propagates_NaN",
        })
    return output
