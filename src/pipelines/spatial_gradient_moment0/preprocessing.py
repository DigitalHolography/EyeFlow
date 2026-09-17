"""Independent float-preserving stages for vessel-wall preprocessing.

Arrays use ``(time, y, x)`` ordering throughout.  The crop stage is an
intentional identity here: the existing velocity-geometry crop is applied by
the caller before this runner, and the existing resize/rotation helpers are
used together by the movable ``rotate`` stage.
"""

from __future__ import annotations

import json
import math
import warnings
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy import ndimage

from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    _center_pad_for_rotation,
    _resize_subimage_stack,
    _rotate_stack_with_nan,
)

from .config import PreprocessingConfig


@dataclass
class PipelineResult:
    data: np.ndarray
    intermediates: list[tuple[str, np.ndarray]] = field(default_factory=list)


STAGE_NAMES = (
    "crop",
    "gaussian3d",
    "temporal_filter",
    "spatial_filter",
    "clahe",
    "dog",
    "non_local_means",
    "bilateral",
    "longitudinal_filter",
    "gradient",
    "rotate",
    "post_gradient_filter",
)


def run_preprocessing_pipeline(
    cropped_stack,
    *,
    rotation_angle: float,
    config: PreprocessingConfig,
    capture_intermediates: bool | None = None,
    rotation_operation: Callable[[np.ndarray, float, str], np.ndarray] | None = None,
) -> PipelineResult:
    """Execute configured stages in order without quantizing scientific data."""

    validate_preprocessing_config(config)
    data = _as_stack(cropped_stack).copy()
    capture = config.output.save_intermediates if capture_intermediates is None else bool(
        capture_intermediates
    )
    intermediates: list[tuple[str, np.ndarray]] = []
    rotate_stage = rotation_operation or (
        lambda values, angle, interpolation: rotate_with_existing_geometry(
            values, angle, interpolation=interpolation
        )
    )
    stage_functions: dict[str, Callable[[np.ndarray], np.ndarray]] = {
        "crop": lambda values: values,
        "gaussian3d": lambda values: gaussian3d(values, config),
        "temporal_filter": lambda values: temporal_filter(values, config),
        "spatial_filter": lambda values: spatial_filter(values, config),
        "clahe": lambda values: clahe(values, config),
        "dog": lambda values: difference_of_gaussians(values, config),
        "non_local_means": lambda values: non_local_means(values, config),
        "bilateral": lambda values: bilateral(values, config),
        "longitudinal_filter": lambda values: longitudinal_filter(values, config),
        "gradient": lambda values: gradient(values, config),
        "rotate": lambda values: rotate_stage(
            values, rotation_angle, config.rotation.interpolation
        ),
        "post_gradient_filter": lambda values: post_gradient_filter(values, config),
    }
    for stage_name in config.pipeline:
        if not _stage_enabled(config, stage_name):
            continue
        data = np.asarray(stage_functions[stage_name](data), dtype=np.float32)
        if data.ndim != 3:
            raise RuntimeError(
                f"Preprocessing stage {stage_name!r} returned shape {data.shape}; "
                "expected (time, y, x)."
            )
        if capture:
            intermediates.append((stage_name, data.copy()))
    return PipelineResult(data=data, intermediates=intermediates)


def validate_preprocessing_config(config: PreprocessingConfig) -> None:
    """Reject ambiguous orderings and invalid stage parameters early."""

    unknown = [stage for stage in config.pipeline if stage not in STAGE_NAMES]
    if unknown:
        raise ValueError(
            f"Unknown preprocessing stage(s): {', '.join(unknown)}. "
            f"Available stages: {', '.join(STAGE_NAMES)}."
        )
    duplicates = sorted({stage for stage in config.pipeline if config.pipeline.count(stage) > 1})
    if duplicates:
        raise ValueError(
            "Stages may occur only once; use the dedicated post-gradient stage instead: "
            + ", ".join(duplicates)
        )
    if not config.pipeline or config.pipeline[0] != "crop" or not config.crop.enabled:
        raise ValueError("The existing velocity-geometry crop must be the first enabled stage.")
    if config.crop.source != "velocity_geometry":
        raise ValueError("crop.source must remain 'velocity_geometry' in this pipeline.")
    if "gradient" not in config.pipeline or not config.gradient.enabled:
        raise ValueError("Exactly one enabled gradient stage is required.")
    if "rotate" not in config.pipeline or not config.rotation.enabled:
        raise ValueError(
            "Exactly one enabled rotate stage is required so profiles match existing geometry."
        )
    if config.rotation.interpolation not in {"existing_bilinear", "bilinear"}:
        raise ValueError(
            "The shared geometry implementation currently supports only its existing "
            "bilinear interpolation."
        )
    if config.gradient.method not in {"sobel", "scharr", "derivative_gaussian"}:
        raise ValueError("gradient.method must be sobel, scharr, or derivative_gaussian.")
    if config.gradient.direction not in {"x", "magnitude"}:
        raise ValueError("gradient.direction must be 'x' or 'magnitude'.")
    if (
        config.gradient.method == "derivative_gaussian"
        and config.gradient.derivative_sigma_x <= 0.0
    ):
        raise ValueError("A derivative-of-Gaussian Gx requires derivative_sigma_x > 0.")
    if (
        config.gradient.method == "derivative_gaussian"
        and config.gradient.direction == "magnitude"
        and config.gradient.derivative_sigma_y <= 0.0
    ):
        raise ValueError(
            "Derivative-of-Gaussian magnitude requires derivative_sigma_y > 0."
        )
    _validate_sigmas(
        (config.gaussian3d.sigma_t, config.gaussian3d.sigma_y, config.gaussian3d.sigma_x),
        "gaussian3d sigma",
    )
    _validate_sigmas(
        (
            config.gradient.derivative_sigma_t,
            config.gradient.derivative_sigma_y,
            config.gradient.derivative_sigma_x,
        ),
        "derivative-of-Gaussian sigma",
    )
    for name, method, allowed in (
        ("temporal_filter", config.temporal_filter.method, {"none", "gaussian", "mean", "median"}),
        ("spatial_filter", config.spatial_filter.method, {"none", "gaussian", "mean", "median"}),
        (
            "post_gradient_filter",
            config.post_gradient_filter.method,
            {"none", "gaussian", "mean", "median"},
        ),
        (
            "longitudinal_filter",
            config.longitudinal_filter.method,
            {"mean", "median", "trimmed_mean"},
        ),
    ):
        if method not in allowed:
            raise ValueError(f"{name}.method must be one of {sorted(allowed)}.")
    for name, value in (
        ("temporal_filter.window", config.temporal_filter.window),
        ("spatial_filter.kernel_x", config.spatial_filter.kernel_x),
        ("spatial_filter.kernel_y", config.spatial_filter.kernel_y),
        ("post_gradient_filter.kernel_x", config.post_gradient_filter.kernel_x),
        ("post_gradient_filter.kernel_y", config.post_gradient_filter.kernel_y),
        ("post_gradient_filter.kernel_t", config.post_gradient_filter.kernel_t),
        ("longitudinal_filter.kernel_y", config.longitudinal_filter.kernel_y),
    ):
        if int(value) != value or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    if not 0.0 <= config.longitudinal_filter.trim_fraction < 0.5:
        raise ValueError("longitudinal_filter.trim_fraction must be in [0, 0.5).")
    for name, truncate in (
        ("gaussian3d.truncate", config.gaussian3d.truncate),
        ("temporal_filter.truncate", config.temporal_filter.truncate),
        ("spatial_filter.truncate", config.spatial_filter.truncate),
        ("gradient.truncate", config.gradient.truncate),
        ("dog.truncate", config.dog.truncate),
        ("post_gradient_filter.truncate", config.post_gradient_filter.truncate),
    ):
        if not math.isfinite(float(truncate)) or float(truncate) <= 0:
            raise ValueError(f"{name} must be finite and greater than zero.")


def gaussian3d(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.gaussian3d
    return _nan_gaussian_filter(
        stack,
        sigma=(settings.sigma_t, settings.sigma_y, settings.sigma_x),
        truncate=settings.truncate,
    )


def temporal_filter(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.temporal_filter
    if settings.method == "none":
        return stack.copy()
    if settings.method == "gaussian":
        _validate_sigmas((settings.sigma,), "temporal_filter.sigma")
        return _nan_gaussian_filter(
            stack, sigma=(settings.sigma, 0.0, 0.0), truncate=settings.truncate
        )
    size = (settings.window, 1, 1)
    if settings.method == "mean":
        return _nan_uniform_filter(stack, size=size)
    return _nan_median_filter(stack, size=size)


def spatial_filter(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.spatial_filter
    if settings.method == "none":
        return stack.copy()
    if settings.method == "gaussian":
        _validate_sigmas((settings.sigma_y, settings.sigma_x), "spatial_filter sigma")
        return _nan_gaussian_filter(
            stack,
            sigma=(0.0, settings.sigma_y, settings.sigma_x),
            truncate=settings.truncate,
        )
    size = (1, settings.kernel_y, settings.kernel_x)
    if settings.method == "mean":
        return _nan_uniform_filter(stack, size=size)
    return _nan_median_filter(stack, size=size)


def gradient(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.gradient
    values = _as_stack(stack)
    finite = np.isfinite(values)
    filled = np.where(finite, values, np.float32(0.0))
    if settings.method == "sobel":
        gx = np.stack(
            [ndimage.sobel(frame, axis=1, mode="nearest") for frame in filled]
        )
        gy = (
            np.stack(
                [ndimage.sobel(frame, axis=0, mode="nearest") for frame in filled]
            )
            if settings.direction == "magnitude"
            else None
        )
    elif settings.method == "scharr":
        kernel_x = np.asarray(
            [[-3.0, 0.0, 3.0], [-10.0, 0.0, 10.0], [-3.0, 0.0, 3.0]],
            dtype=np.float32,
        )
        gx = ndimage.correlate(filled, kernel_x[None, :, :], mode="nearest")
        gy = (
            ndimage.correlate(filled, kernel_x.T[None, :, :], mode="nearest")
            if settings.direction == "magnitude"
            else None
        )
    else:
        sigma = (
            settings.derivative_sigma_t,
            settings.derivative_sigma_y,
            settings.derivative_sigma_x,
        )
        gx = _nan_gaussian_derivative(
            values,
            sigma=sigma,
            axis=2,
            truncate=settings.truncate,
        )
        gy = (
            _nan_gaussian_derivative(
                values,
                sigma=sigma,
                axis=1,
                truncate=settings.truncate,
            )
            if settings.direction == "magnitude"
            else None
        )
    result = np.hypot(gx, gy, dtype=np.float32) if gy is not None else gx.astype(np.float32)
    if settings.direction == "x" and not settings.signed:
        result = np.abs(result, dtype=np.float32)
    # Missing crop/rotation padding must never become a strong synthetic wall.
    valid_neighborhood = ndimage.minimum_filter(finite, size=(1, 3, 3), mode="nearest")
    result[~valid_neighborhood] = np.nan
    return result.astype(np.float32, copy=False)


def difference_of_gaussians(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.dog
    small_sigma = (
        settings.sigma_small_t,
        settings.sigma_small_y,
        settings.sigma_small_x,
    )
    large_sigma = (
        settings.sigma_large_t,
        settings.sigma_large_y,
        settings.sigma_large_x,
    )
    _validate_sigmas(small_sigma, "dog small sigma")
    _validate_sigmas(large_sigma, "dog large sigma")
    small = _nan_gaussian_filter(stack, sigma=small_sigma, truncate=settings.truncate)
    large = _nan_gaussian_filter(stack, sigma=large_sigma, truncate=settings.truncate)
    return np.subtract(small, large, dtype=np.float32)


def clahe(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    from skimage import exposure

    settings = config.clahe
    result = np.full(stack.shape, np.nan, dtype=np.float32)
    for index, frame in enumerate(stack):
        finite = np.isfinite(frame)
        if not np.any(finite):
            continue
        low = float(np.min(frame[finite]))
        high = float(np.max(frame[finite]))
        if high <= low:
            result[index, finite] = 0.0
            continue
        normalized = np.zeros(frame.shape, dtype=np.float32)
        normalized[finite] = (frame[finite] - low) / (high - low)
        equalized = exposure.equalize_adapthist(
            normalized,
            kernel_size=(settings.tile_y, settings.tile_x),
            clip_limit=settings.clip_limit,
        )
        result[index, finite] = np.asarray(equalized, dtype=np.float32)[finite]
    return result


def non_local_means(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    from skimage.restoration import denoise_nl_means

    settings = config.non_local_means
    return _per_frame_preserving_nan(
        stack,
        lambda frame: denoise_nl_means(
            frame,
            patch_size=settings.patch_size,
            patch_distance=settings.patch_distance,
            h=settings.h,
            fast_mode=settings.fast_mode,
            preserve_range=True,
            channel_axis=None,
        ),
    )


def bilateral(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    from skimage.restoration import denoise_bilateral

    settings = config.bilateral
    return _per_frame_preserving_nan(
        stack,
        lambda frame: denoise_bilateral(
            frame,
            win_size=settings.window_size,
            sigma_color=settings.sigma_color,
            sigma_spatial=settings.sigma_spatial,
            mode="edge",
            channel_axis=None,
        ),
    )


def longitudinal_filter(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.longitudinal_filter
    size = (1, settings.kernel_y, 1)
    if settings.method == "mean":
        return _nan_uniform_filter(stack, size=size)
    if settings.method == "median":
        return _nan_median_filter(stack, size=size)

    def trimmed(values):
        finite = np.sort(values[np.isfinite(values)])
        if finite.size == 0:
            return np.nan
        trim = int(math.floor(settings.trim_fraction * finite.size))
        kept = finite[trim : finite.size - trim] if trim else finite
        return float(np.mean(kept, dtype=np.float64))

    filtered = ndimage.generic_filter(
        _as_stack(stack), trimmed, size=size, mode="nearest"
    )
    filtered[~np.isfinite(stack)] = np.nan
    return filtered.astype(np.float32, copy=False)


def post_gradient_filter(stack: np.ndarray, config: PreprocessingConfig) -> np.ndarray:
    settings = config.post_gradient_filter
    if settings.method == "none":
        return stack.copy()
    if settings.method == "gaussian":
        sigma = (settings.sigma_t, settings.sigma_y, settings.sigma_x)
        _validate_sigmas(sigma, "post_gradient_filter sigma")
        return _nan_gaussian_filter(stack, sigma=sigma, truncate=settings.truncate)
    size = (settings.kernel_t, settings.kernel_y, settings.kernel_x)
    if settings.method == "mean":
        return _nan_uniform_filter(stack, size=size)
    return _nan_median_filter(stack, size=size)


def rotate_with_existing_geometry(
    stack: np.ndarray,
    angle: float,
    *,
    interpolation: str = "existing_bilinear",
) -> np.ndarray:
    """Reuse the exact crop-resize-pad-rotate implementation used by EyeFlow."""

    if interpolation not in {"existing_bilinear", "bilinear"}:
        raise ValueError("Only the existing bilinear geometry interpolation is supported.")
    resized = _resize_subimage_stack(_as_stack(stack))
    padded = _center_pad_for_rotation(resized, np.nan)
    return _rotate_stack_with_nan(padded, float(angle)).astype(np.float32, copy=False)


def processing_metadata(config: PreprocessingConfig) -> dict[str, object]:
    gradient_settings = config.gradient
    operator = gradient_settings.method
    if operator == "derivative_gaussian":
        operator = (
            "derivative of Gaussian "
            f"(sigma_t={gradient_settings.derivative_sigma_t:g}, "
            f"sigma_x={gradient_settings.derivative_sigma_x:g}, "
            f"sigma_y={gradient_settings.derivative_sigma_y:g})"
        )
    elif operator == "sobel":
        operator = "3x3 Sobel"
    else:
        operator = "3x3 Scharr"
    if gradient_settings.direction == "magnitude":
        operator += " magnitude"
    elif gradient_settings.signed:
        operator += " signed Gx"
    else:
        operator += " absolute Gx"
    return {
        "measurement": (
            "spatial_gradient_magnitude"
            if gradient_settings.direction == "magnitude"
            else "spatial_gradient_x"
        ),
        "preprocessing_preset": config.name,
        "processing_order": "_".join(config.pipeline),
        "spatial_operator": operator,
        "gradient_method": gradient_settings.method,
        "gradient_direction": gradient_settings.direction,
        "gradient_signed": bool(gradient_settings.signed),
        "interpolation": config.rotation.interpolation,
        "preprocessing_config_json": json.dumps(
            config.to_dict(), sort_keys=True, separators=(",", ":")
        ),
        "interpolation_order": (
            "gradient_before_rotation"
            if config.pipeline.index("gradient") < config.pipeline.index("rotate")
            else "rotation_before_gradient"
        ),
    }


def _stage_enabled(config: PreprocessingConfig, stage: str) -> bool:
    settings_name = "rotation" if stage == "rotate" else stage
    return bool(getattr(config, settings_name).enabled)


def _nan_gaussian_filter(stack, *, sigma, truncate: float) -> np.ndarray:
    values = _as_stack(stack)
    sigma = tuple(float(item) for item in sigma)
    _validate_sigmas(sigma, "Gaussian sigma")
    if not any(sigma):
        return values.copy()
    finite = np.isfinite(values)
    kwargs = dict(sigma=sigma, mode="nearest", truncate=float(truncate))
    filtered_values = ndimage.gaussian_filter(np.where(finite, values, 0.0), **kwargs)
    weights = ndimage.gaussian_filter(finite.astype(np.float32), **kwargs)
    result = np.full(values.shape, np.nan, dtype=np.float32)
    np.divide(filtered_values, weights, out=result, where=finite & (weights > 0.0))
    return result


def _nan_uniform_filter(stack, *, size) -> np.ndarray:
    values = _as_stack(stack)
    finite = np.isfinite(values)
    size = tuple(int(item) for item in size)
    filtered_values = ndimage.uniform_filter(
        np.where(finite, values, 0.0), size=size, mode="nearest"
    )
    weights = ndimage.uniform_filter(finite.astype(np.float32), size=size, mode="nearest")
    result = np.full(values.shape, np.nan, dtype=np.float32)
    np.divide(filtered_values, weights, out=result, where=finite & (weights > 0.0))
    return result


def _nan_gaussian_derivative(stack, *, sigma, axis: int, truncate: float) -> np.ndarray:
    """Differentiate a NaN-normalized Gaussian convolution along one axis."""

    values = _as_stack(stack)
    finite = np.isfinite(values)
    weights_input = finite.astype(np.float32)
    filled = np.where(finite, values, np.float32(0.0))
    order = [0, 0, 0]
    order[axis] = 1
    kwargs = dict(sigma=sigma, mode="nearest", truncate=float(truncate))
    smooth_values = ndimage.gaussian_filter(filled, order=0, **kwargs)
    smooth_weights = ndimage.gaussian_filter(weights_input, order=0, **kwargs)
    derivative_values = ndimage.gaussian_filter(filled, order=order, **kwargs)
    derivative_weights = ndimage.gaussian_filter(weights_input, order=order, **kwargs)
    numerator = derivative_values * smooth_weights - smooth_values * derivative_weights
    denominator = np.square(smooth_weights)
    result = np.full(values.shape, np.nan, dtype=np.float32)
    valid = finite & (denominator > np.finfo(np.float32).eps)
    np.divide(numerator, denominator, out=result, where=valid)
    return result


def _nan_median_filter(stack, *, size) -> np.ndarray:
    values = _as_stack(stack)
    size = tuple(int(item) for item in size)
    if not np.any(~np.isfinite(values)):
        return ndimage.median_filter(values, size=size, mode="nearest").astype(
            np.float32, copy=False
        )
    finite = np.isfinite(values)
    if size[1:] == (1, 1) and np.all(finite == finite[0:1]):
        filtered = ndimage.median_filter(
            np.where(finite, values, np.float32(0.0)),
            size=size,
            mode="nearest",
        ).astype(np.float32, copy=False)
        filtered[~finite] = np.nan
        return filtered
    if size[1:] == (1, 1):
        radius_before = size[0] // 2
        radius_after = size[0] - radius_before - 1
        padded = np.pad(
            values,
            ((radius_before, radius_after), (0, 0), (0, 0)),
            mode="edge",
        )
        result = np.empty_like(values)
        for index in range(len(values)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                result[index] = np.nanmedian(padded[index : index + size[0]], axis=0)
        result[~finite] = np.nan
        return result
    filtered = ndimage.generic_filter(
        values,
        lambda samples: np.nanmedian(samples) if np.any(np.isfinite(samples)) else np.nan,
        size=size,
        mode="nearest",
    )
    filtered[~finite] = np.nan
    return filtered.astype(np.float32, copy=False)


def _per_frame_preserving_nan(stack, operation) -> np.ndarray:
    values = _as_stack(stack)
    result = np.full(values.shape, np.nan, dtype=np.float32)
    for index, frame in enumerate(values):
        finite = np.isfinite(frame)
        if not np.any(finite):
            continue
        fill_value = float(np.median(frame[finite]))
        filtered = operation(np.where(finite, frame, fill_value))
        result[index, finite] = np.asarray(filtered, dtype=np.float32)[finite]
    return result


def _validate_sigmas(sigmas, name: str) -> None:
    for sigma in sigmas:
        if not math.isfinite(float(sigma)) or float(sigma) < 0.0:
            raise ValueError(f"{name} values must be finite and non-negative.")


def _as_stack(stack) -> np.ndarray:
    values = np.asarray(stack, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError(f"Preprocessing input must have shape (time, y, x), got {values.shape}.")
    return values


__all__ = [
    "PipelineResult",
    "STAGE_NAMES",
    "processing_metadata",
    "run_preprocessing_pipeline",
    "validate_preprocessing_config",
]
