"""Photometric preprocessing for retinal registration."""

from __future__ import annotations

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None

from .parameters import PhotometricConfig


def minmax_float01(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(image)
    if not finite.any():
        return np.zeros(image.shape, np.float32)

    minimum = float(np.min(image[finite]))
    maximum = float(np.max(image[finite]))
    if maximum <= minimum + 1e-12:
        return np.zeros(image.shape, np.float32)

    output = (image - minimum) / (maximum - minimum)
    return np.clip(output, 0.0, 1.0).astype(np.float32)


def gray_float01(image: np.ndarray) -> np.ndarray:
    if image.ndim == 3:
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    else:
        gray = np.asarray(image)
    if gray.dtype == np.uint8:
        return gray.astype(np.float32) / 255.0
    arr = np.asarray(gray, np.float32)
    finite = np.isfinite(arr)
    if not finite.any():
        return np.zeros(arr.shape, np.float32)
    fill = float(np.nanmedian(arr[finite]))
    arr = np.nan_to_num(arr, nan=fill, posinf=fill, neginf=fill)
    if float(arr.min()) >= 0.0 and float(arr.max()) <= 1.0 + 1e-6:
        return np.clip(arr, 0.0, 1.0).astype(np.float32)
    return minmax_float01(arr)


def robust_signed_to_unit(
    values: np.ndarray,
    valid_mask: np.ndarray | None,
    clip: float,
) -> np.ndarray:
    arr = np.asarray(values, np.float32)
    finite = np.isfinite(arr)
    if valid_mask is not None:
        finite &= np.asarray(valid_mask, bool)
    sample = arr[finite]
    if sample.size == 0:
        return np.full(arr.shape, 0.5, np.float32)
    median = float(np.median(sample))
    mad = float(np.median(np.abs(sample - median)))
    scale = max(1.4826 * mad, float(np.std(sample)) * 0.05, 1e-6)
    z = np.clip((arr - median) / scale, -clip, clip)
    output = 0.5 + 0.5 * z / clip
    output = np.nan_to_num(output, nan=0.5, posinf=1.0, neginf=0.0)
    if valid_mask is not None:
        output = output.astype(np.float32, copy=False)
        output[~np.asarray(valid_mask, bool)] = 0.5
    return np.clip(output, 0.0, 1.0).astype(np.float32)


def photometric_normalize(
    image: np.ndarray,
    config: PhotometricConfig,
    valid_mask: np.ndarray | None,
) -> np.ndarray:
    gray = gray_float01(image)
    if config.mode == "none":
        return robust_signed_to_unit(gray, valid_mask, config.clip)

    eps = np.float32(1.0 / 255.0)
    homomorphic: np.ndarray | None = None
    local_contrast: np.ndarray | None = None

    if config.mode in {"homomorphic", "hybrid"}:
        log_image = np.log(np.maximum(gray, eps)).astype(np.float32)
        if config.illumination_sigma > 0:
            low_frequency = cv2.GaussianBlur(
                log_image,
                (0, 0),
                config.illumination_sigma,
                borderType=cv2.BORDER_REFLECT101,
            )
            homomorphic = log_image - low_frequency
        else:
            homomorphic = log_image

    if config.mode in {"local_contrast", "hybrid"}:
        sigma = max(config.local_contrast_sigma, 0.5)
        local_mean = cv2.GaussianBlur(gray, (0, 0), sigma, borderType=cv2.BORDER_REFLECT101)
        local_second = cv2.GaussianBlur(
            gray * gray, (0, 0), sigma, borderType=cv2.BORDER_REFLECT101
        )
        local_std = np.sqrt(np.maximum(local_second - local_mean * local_mean, 0.0))
        noise_floor = max(float(np.percentile(local_std, 20.0)) * 0.35, 1.0 / 255.0)
        local_contrast = (gray - local_mean) / np.maximum(local_std, noise_floor)

    if config.mode == "homomorphic":
        assert homomorphic is not None
        return robust_signed_to_unit(homomorphic, valid_mask, config.clip)
    if config.mode == "local_contrast":
        assert local_contrast is not None
        return robust_signed_to_unit(local_contrast, valid_mask, config.clip)

    assert homomorphic is not None and local_contrast is not None
    homo_unit = robust_signed_to_unit(homomorphic, valid_mask, config.clip)
    local_unit = robust_signed_to_unit(local_contrast, valid_mask, config.clip)
    combined = 0.55 * (homo_unit - 0.5) + 0.45 * (local_unit - 0.5)
    output = np.clip(0.5 + combined, 0.0, 1.0).astype(np.float32)
    if valid_mask is not None:
        output[~np.asarray(valid_mask, bool)] = 0.5
    return output


def structural_registration_signal(
    image: np.ndarray,
    config: PhotometricConfig,
    valid_mask: np.ndarray,
    edge_weight: float,
) -> np.ndarray:
    structural = photometric_normalize(image, config, valid_mask)
    grad_x = cv2.Sobel(structural, cv2.CV_32F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(structural, cv2.CV_32F, 0, 1, ksize=3)
    edges = np.hypot(grad_x, grad_y).astype(np.float32)
    edge_unit = robust_signed_to_unit(edges, valid_mask, config.clip)

    signal = (1.0 - edge_weight) * structural + edge_weight * edge_unit
    signal = np.clip(signal, 0.0, 1.0).astype(np.float32)
    signal[~np.asarray(valid_mask, bool)] = 0.5
    return signal


def photometric_confidence(
    reference: np.ndarray,
    valid_mask: np.ndarray,
    config: PhotometricConfig,
) -> np.ndarray:
    gray = gray_float01(reference)
    valid = np.asarray(valid_mask, bool)
    if not valid.any():
        return np.zeros(gray.shape, np.float32)
    sigma_bg = max(config.illumination_sigma, 1.0)
    background = cv2.GaussianBlur(gray, (0, 0), sigma_bg, borderType=cv2.BORDER_REFLECT101)
    sigma_local = max(config.local_contrast_sigma, 1.0)
    mean = cv2.GaussianBlur(gray, (0, 0), sigma_local, borderType=cv2.BORDER_REFLECT101)
    second = cv2.GaussianBlur(gray * gray, (0, 0), sigma_local, borderType=cv2.BORDER_REFLECT101)
    local_std = np.sqrt(np.maximum(second - mean * mean, 0.0))

    bg_values = background[valid]
    dark = float(np.percentile(bg_values, config.dark_percentile))
    median = float(np.percentile(bg_values, 50.0))
    brightness = np.clip((background - dark) / max(median - dark, 1e-6), 0.0, 1.0)
    std_reference = max(float(np.percentile(local_std[valid], 70.0)), 1e-6)
    texture = np.clip(local_std / std_reference, 0.0, 1.0)
    quality = np.sqrt(np.clip(brightness * (0.25 + 0.75 * texture), 0.0, 1.0))
    confidence = config.confidence_floor + (1.0 - config.confidence_floor) * quality
    confidence = confidence.astype(np.float32)
    confidence[~valid] = 0.0
    return np.clip(confidence, 0.0, 1.0)


def soft_mask(binary_mask: np.ndarray, sigma: float) -> np.ndarray:
    mask = (binary_mask > 0).astype(np.float32)
    if sigma > 0:
        mask = cv2.GaussianBlur(mask, (0, 0), sigma, borderType=cv2.BORDER_REFLECT101)
    peak = float(mask.max())
    if peak > 0:
        mask /= peak
    return np.clip(mask, 0.0, 1.0).astype(np.float32)
