"""Peak, valley, and inclusive-area metrics for the 32-bin gap spectrum."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from utils.logger import Logger

SPECTRUM_METRIC_DEFINITIONS = {
    "index_max1": "zero-based index of the maximum in bins 0 through 15",
    "max1": "maximum spectrum value in bins 0 through 15, ignoring NaNs",
    "index_max2": "zero-based index of the maximum in bins 16 through 31",
    "max2": "maximum spectrum value in bins 16 through 31, ignoring NaNs",
    "index_min1": "fixed zero-based valley index 15",
    "min1": "spectrum value at zero-based index 15",
    "index_min2": "fixed zero-based valley index 16",
    "min2": "spectrum value at zero-based index 16",
    "h0": "arithmetic mean of max1 and max2; missing inputs propagate NaN",
    "h1": "arithmetic mean of min1 and min2; missing inputs propagate NaN",
    "rho_h": "h1 divided by h0; zero or missing h0 produces NaN",
    "A0": "h0 times the number of bins from index_max1 through index_max2, inclusive",
    "A1": "unweighted spectrum sum from index_max1 through index_max2, inclusive; NaNs propagate",
    "rho_A": "A1 divided by A0; zero or missing A0 produces NaN",
}


def analyze_spectrum(spectrum) -> dict[str, float]:
    """Analyze a 32-bin spectrum in its stored frequency order.

    Select one maximum in each half, including endpoints and ignoring NaNs.
    Ties select the lowest index. An all-NaN half has a NaN index and value.
    Valley indexes remain 15 and 16 for a valid shape, even if values are NaN.
    Means require both inputs. Areas include both peak bins without wrapping,
    interpolation, or frequency-spacing weights. A1 propagates interval NaNs.
    Unsupported shapes return all fourteen metrics as NaN with a warning.
    """
    values = np.asarray(spectrum)
    result = {name: float("nan") for name in SPECTRUM_METRIC_DEFINITIONS}
    if values.shape != (32,):
        Logger.log_warning(
            f"Gap spectrum metrics require exactly 32 bins; received shape {values.shape}. "
            "Returning NaN for all spectrum metrics."
        )
        return result
    values = np.asarray(values, dtype=np.float64)
    result["index_max1"], result["max1"] = _peak(values, 0, 16)
    result["index_max2"], result["max2"] = _peak(values, 16, 32)
    result.update(index_min1=15.0, min1=float(values[15]), index_min2=16.0, min2=float(values[16]))
    result["h0"] = (result["max1"] + result["max2"]) / 2
    result["h1"] = (result["min1"] + result["min2"]) / 2
    if result["h0"] != 0:
        result["rho_h"] = result["h1"] / result["h0"]
    if np.isfinite(result["index_max1"]) and np.isfinite(result["index_max2"]):
        start = int(result["index_max1"])
        stop = int(result["index_max2"]) + 1
        result["A0"] = result["h0"] * (stop - start)
        result["A1"] = float(np.sum(values[start:stop]))
        if result["A0"] != 0:
            result["rho_A"] = result["A1"] / result["A0"]
    return result


def _peak(values: np.ndarray, start: int, stop: int) -> tuple[float, float]:
    half = values[start:stop]
    valid = np.flatnonzero(~np.isnan(half))
    if valid.size == 0:
        return float("nan"), float("nan")
    index = start + int(valid[np.argmax(half[valid])])
    return float(index), float(values[index])


def pack_spectrum_metrics(spectrum, *, source_attrs: dict[str, object]) -> dict[str, object]:
    """Pack scalar metrics with units and references to the aggregate spectrum."""
    schema = EyeFlowOutputPaths.active()
    result = analyze_spectrum(spectrum)
    attrs = {
        "source_path": "/" + schema.spectrograms_masked_time_aggregated_median,
        "source_bin_count": int(np.size(spectrum)),
        "expected_bin_count": 32,
        "dimDesc": [],
        "peak_selection": "maximum_per_16_bin_half_ignore_NaN_lowest_index_on_ties",
        "interval_endpoints": "both_included",
        "frequency_bin_order": source_attrs.get("frequency_bin_order", "HD_export_order"),
        "nonfinite_value_policy": "missing_required_values_and_interval_NaNs_propagate",
    }
    output = {}
    for name, value in result.items():
        unit = (
            "1"
            if name.startswith("index_") or name in {"rho_h", "rho_A"}
            else source_attrs.get("unit", "a.u.")
        )
        metric_attrs = {**attrs, "unit": unit, "definition": SPECTRUM_METRIC_DEFINITIONS[name]}
        if name.startswith("index_"):
            metric_attrs["index_base"] = 0
            metric_attrs["index_storage"] = "integer_valued_float32_with_NaN_for_missing_index"
        output[f"{schema.spectrograms_masked_time_aggregated_metrics_root}/{name}"] = (
            np.float32(value),
            metric_attrs,
        )
    return output
