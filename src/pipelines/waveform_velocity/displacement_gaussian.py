"""Two-Gaussian fitting for time-meaned transverse displacement profiles."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from pipeline_engine.base import DatasetValue


_FWHM_FACTOR = np.float32(2.355)
_SQRT_TWO_PI = np.float32(np.sqrt(2.0 * np.pi))
_GAUSSIAN_METRIC_NAMES = (
    "Gaussian_Baseline",
    "Gaussian_A_L",
    "Gaussian_A_R",
    "Gaussian_Mu_L",
    "Gaussian_Mu_R",
    "Gaussian_Sigma_L",
    "Gaussian_Sigma_R",
    "Gaussian_FWHM_L",
    "Gaussian_FWHM_R",
    "Gaussian_Area_L",
    "Gaussian_Area_R",
    "Gaussian_Peak_Separation",
    "Gaussian_RMSE",
)


@dataclass(frozen=True)
class _CurveFitResult:
    fitted: np.ndarray
    metrics: tuple[np.float32, ...]
    success: bool
    initialization_complete: bool


def fit_two_gaussian_profiles(
    profile: DatasetValue,
    max_x_position: DatasetValue,
    max_y_position: DatasetValue,
) -> tuple[DatasetValue, dict[str, DatasetValue]]:
    """Fit two ordered Gaussian components to every time-meaned profile."""

    values, peak_x, peak_y, spatial_axis = _validated_inputs(
        profile,
        max_x_position,
        max_y_position,
    )
    fitted_profiles = np.full(values.shape, np.nan, dtype=np.float32)
    metric_shape = values.shape[1:]
    metric_data = {
        name: np.full(metric_shape, np.nan, dtype=np.float32)
        for name in _GAUSSIAN_METRIC_NAMES
    }
    success = np.zeros(metric_shape, dtype=bool)
    initialization_complete = np.zeros(metric_shape, dtype=bool)

    for beat_index, branch_index, radius_index in np.ndindex(metric_shape):
        curve = values[:, beat_index, branch_index, radius_index]
        result = _fit_curve(
            curve,
            peak_x[:, beat_index, branch_index, radius_index],
            peak_y[:, beat_index, branch_index, radius_index],
        )
        fitted_profiles[:, beat_index, branch_index, radius_index] = result.fitted
        for name, metric in zip(_GAUSSIAN_METRIC_NAMES, result.metrics):
            metric_data[name][beat_index, branch_index, radius_index] = metric
        success[beat_index, branch_index, radius_index] = result.success
        initialization_complete[
            beat_index, branch_index, radius_index
        ] = result.initialization_complete

    profile_attrs = dict(profile.attrs or {})
    fitted_dataset = DatasetValue(
        data=fitted_profiles,
        attrs={
            "unit": profile_attrs.get("unit", "pixels"),
            "dimDesc": [spatial_axis, "beat", "branch", "radius"],
            "model": (
                "B + A_L*exp(-(r-mu_L)^2/(2*sigma_L^2)) + "
                "A_R*exp(-(r-mu_R)^2/(2*sigma_R^2))"
            ),
            "fit_axis": spatial_axis,
            "nan_mask": "preserved_from_source_profile",
        },
        h5_options=_profile_h5_options(fitted_profiles.shape),
    )

    metric_units = {
        "Gaussian_Baseline": "pixels",
        "Gaussian_A_L": "pixels",
        "Gaussian_A_R": "pixels",
        "Gaussian_Mu_L": "pixels",
        "Gaussian_Mu_R": "pixels",
        "Gaussian_Sigma_L": "pixels",
        "Gaussian_Sigma_R": "pixels",
        "Gaussian_FWHM_L": "pixels",
        "Gaussian_FWHM_R": "pixels",
        "Gaussian_Area_L": "pixels^2",
        "Gaussian_Area_R": "pixels^2",
        "Gaussian_Peak_Separation": "pixels",
        "Gaussian_RMSE": "pixels",
    }
    metrics = {
        name: _metric_dataset(data, unit=metric_units[name])
        for name, data in metric_data.items()
    }
    metrics["Gaussian_Fit_Success"] = _metric_dataset(success, unit="1")
    metrics["Gaussian_Initialization_Complete"] = _metric_dataset(
        initialization_complete,
        unit="1",
    )
    return fitted_dataset, metrics


def _validated_inputs(
    profile: DatasetValue,
    max_x_position: DatasetValue,
    max_y_position: DatasetValue,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    values = np.asarray(profile.data, dtype=np.float32)
    profile_dims = list((profile.attrs or {}).get("dimDesc", ()))
    peak_x = np.asarray(max_x_position.data, dtype=np.float32)
    peak_y = np.asarray(max_y_position.data, dtype=np.float32)
    peak_dims = ["peak", "beat", "branch", "radius"]
    if (
        values.ndim != 4
        or len(profile_dims) != 4
        or profile_dims[1:] != ["beat", "branch", "radius"]
    ):
        raise ValueError(
            "Gaussian source profile must have dimensions "
            "(spatial_sample, beat, branch, radius)."
        )
    expected_peak_shape = (2, *values.shape[1:])
    if (
        peak_x.shape != expected_peak_shape
        or list((max_x_position.attrs or {}).get("dimDesc", ())) != peak_dims
        or peak_y.shape != expected_peak_shape
        or list((max_y_position.attrs or {}).get("dimDesc", ())) != peak_dims
    ):
        raise ValueError(
            "Gaussian peak initialization must have dimensions "
            "(two_peaks, beat, branch, radius)."
        )
    return values, peak_x, peak_y, profile_dims[0]


def _fit_curve(
    curve: np.ndarray,
    peak_x: np.ndarray,
    peak_y: np.ndarray,
) -> _CurveFitResult:
    values = np.asarray(curve, dtype=np.float32)
    fitted = np.full(values.shape, np.nan, dtype=np.float32)
    empty_metrics = tuple(np.float32(np.nan) for _ in _GAUSSIAN_METRIC_NAMES)
    finite = np.isfinite(values)
    valid_x = np.flatnonzero(finite).astype(np.float64)
    valid_y = values[finite].astype(np.float64)
    initialization_complete = False
    if (
        valid_x.size < 10
        or not np.all(np.isfinite(peak_x))
        or not np.all(np.isfinite(peak_y))
    ):
        return _CurveFitResult(fitted, empty_metrics, False, False)

    r_min = float(valid_x[0])
    r_max = float(valid_x[-1])
    span = r_max - r_min
    mu_l, mu_r = (float(peak_x[0]), float(peak_x[1]))
    if span <= 0.0 or not r_min <= mu_l < mu_r <= r_max:
        return _CurveFitResult(fitted, empty_metrics, False, False)

    baseline = max(0.0, float(np.percentile(valid_y, 5.0)))
    amplitude_l = max(0.0, float(peak_y[0]) - baseline)
    amplitude_r = max(0.0, float(peak_y[1]) - baseline)
    sigma_l, complete_l = _initial_sigma(values, mu_l, float(peak_y[0]), baseline)
    sigma_r, complete_r = _initial_sigma(values, mu_r, float(peak_y[1]), baseline)
    initialization_complete = complete_l and complete_r

    position_logits = _ordered_position_logits(mu_l, mu_r, r_min, r_max)
    initial = np.asarray(
        [
            baseline,
            amplitude_l,
            amplitude_r,
            position_logits[0],
            position_logits[1],
            sigma_l,
            sigma_r,
        ],
        dtype=np.float64,
    )
    minimum_sigma = max(span * 1.0e-5, 1.0e-3)
    lower = np.asarray([0.0, 0.0, 0.0, -20.0, -20.0, minimum_sigma, minimum_sigma])
    upper = np.asarray(
        [np.inf, np.inf, np.inf, 20.0, 20.0, np.inf, np.inf]
    )
    initial = np.minimum(np.maximum(initial, lower), upper)

    def residual(parameters: np.ndarray) -> np.ndarray:
        return _model_from_optimized(parameters, valid_x, r_min, r_max) - valid_y

    try:
        optimized = least_squares(
            residual,
            initial,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            max_nfev=2000,
        )
    except (FloatingPointError, RuntimeError, ValueError):
        return _CurveFitResult(
            fitted,
            empty_metrics,
            False,
            initialization_complete,
        )
    if not optimized.success or not np.all(np.isfinite(optimized.x)):
        return _CurveFitResult(
            fitted,
            empty_metrics,
            False,
            initialization_complete,
        )

    parameters = optimized.x
    fitted_valid = _model_from_optimized(parameters, valid_x, r_min, r_max)
    if not np.all(np.isfinite(fitted_valid)):
        return _CurveFitResult(
            fitted,
            empty_metrics,
            False,
            initialization_complete,
        )
    fitted[finite] = fitted_valid.astype(np.float32)
    baseline, amplitude_l, amplitude_r = parameters[:3]
    mu_l, mu_r = _ordered_means(parameters[3], parameters[4], r_min, r_max)
    sigma_l, sigma_r = parameters[5:7]
    fwhm_l = float(_FWHM_FACTOR) * sigma_l
    fwhm_r = float(_FWHM_FACTOR) * sigma_r
    area_l = amplitude_l * sigma_l * float(_SQRT_TWO_PI)
    area_r = amplitude_r * sigma_r * float(_SQRT_TWO_PI)
    rmse = float(np.sqrt(np.mean(np.square(fitted_valid - valid_y))))
    metrics = tuple(
        np.float32(value)
        for value in (
            baseline,
            amplitude_l,
            amplitude_r,
            mu_l,
            mu_r,
            sigma_l,
            sigma_r,
            fwhm_l,
            fwhm_r,
            area_l,
            area_r,
            mu_r - mu_l,
            rmse,
        )
    )
    return _CurveFitResult(fitted, metrics, True, initialization_complete)


def _model_from_optimized(
    parameters: np.ndarray,
    x_values: np.ndarray,
    r_min: float,
    r_max: float,
) -> np.ndarray:
    baseline, amplitude_l, amplitude_r = parameters[:3]
    mu_l, mu_r = _ordered_means(parameters[3], parameters[4], r_min, r_max)
    sigma_l, sigma_r = parameters[5:7]
    return (
        baseline
        + amplitude_l * np.exp(-np.square(x_values - mu_l) / (2.0 * sigma_l**2))
        + amplitude_r * np.exp(-np.square(x_values - mu_r) / (2.0 * sigma_r**2))
    )


def _ordered_means(
    left_logit: float,
    separation_logit: float,
    r_min: float,
    r_max: float,
) -> tuple[float, float]:
    logits = np.asarray([left_logit, separation_logit, 0.0], dtype=np.float64)
    weights = np.exp(logits - np.max(logits))
    fractions = weights / np.sum(weights)
    span = r_max - r_min
    mu_l = r_min + span * float(fractions[0])
    mu_r = mu_l + span * float(fractions[1])
    return mu_l, mu_r


def _ordered_position_logits(
    mu_l: float,
    mu_r: float,
    r_min: float,
    r_max: float,
) -> tuple[float, float]:
    span = r_max - r_min
    fractions = np.asarray(
        [mu_l - r_min, mu_r - mu_l, r_max - mu_r],
        dtype=np.float64,
    ) / span
    fractions = np.maximum(fractions, 1.0e-6)
    fractions /= np.sum(fractions)
    return (
        float(np.log(fractions[0] / fractions[2])),
        float(np.log(fractions[1] / fractions[2])),
    )


def _initial_sigma(
    curve: np.ndarray,
    peak_position: float,
    peak_value: float,
    baseline: float,
) -> tuple[float, bool]:
    values = np.asarray(curve, dtype=np.float64)
    peak_index = int(round(peak_position))
    finite_indexes = np.flatnonzero(np.isfinite(values))
    finite_span = (
        float(finite_indexes[-1] - finite_indexes[0])
        if finite_indexes.size >= 2
        else 1.0
    )
    fallback = max(finite_span / 10.0, 0.5)
    fallback_sigma = fallback / float(_FWHM_FACTOR)
    if (
        peak_index < 0
        or peak_index >= values.size
        or not np.isfinite(peak_value)
        or peak_value <= baseline
    ):
        return fallback_sigma, False

    half_height = baseline + (peak_value - baseline) / 2.0
    left = _nearest_level_intersection(values, peak_index, half_height, -1)
    right = _nearest_level_intersection(values, peak_index, half_height, 1)
    if left is None or right is None or right <= left:
        return fallback_sigma, False
    return max((right - left) / float(_FWHM_FACTOR), 1.0e-3), True


def _nearest_level_intersection(
    values: np.ndarray,
    peak_index: int,
    level: float,
    direction: int,
) -> float | None:
    if direction < 0:
        pair_starts = range(peak_index - 1, -1, -1)
    else:
        pair_starts = range(peak_index, values.size - 1)
    for start in pair_starts:
        left_x = start
        right_x = start + 1
        left_y = values[left_x]
        right_y = values[right_x]
        if not np.isfinite(left_y) or not np.isfinite(right_y):
            continue
        left_delta = left_y - level
        right_delta = right_y - level
        if left_delta == 0.0:
            return float(left_x)
        if right_delta == 0.0:
            return float(right_x)
        if left_delta * right_delta < 0.0:
            fraction = (level - left_y) / (right_y - left_y)
            return float(left_x + fraction)
    return None


def _metric_dataset(data: np.ndarray, *, unit: str) -> DatasetValue:
    values = np.asarray(data)
    return DatasetValue(
        data=values,
        attrs={
            "unit": unit,
            "dimDesc": ["beat", "branch", "radius"],
        },
        h5_options=_profile_h5_options(values.shape),
    )


def _profile_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    options: dict[str, object] = {
        "compression": "gzip",
        "compression_opts": 4,
        "shuffle": True,
    }
    if len(shape) == 4 and all(shape):
        options["chunks"] = (shape[0], 1, 1, 1)
    return options
