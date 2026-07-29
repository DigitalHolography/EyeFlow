"""Chunked DopplerView-compatible flat-field correction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi


@dataclass(frozen=True)
class FlatFieldParameters:
    gaussian_width: float
    border_amount: float
    input_min: float
    input_max: float
    normalize_input: bool
    output_scale: float


def fit_flat_field_parameters(
    volume,
    *,
    gaussian_width: float,
    border_amount: float = 0.15,
    frame_chunk_size: int = 8,
) -> FlatFieldParameters:
    """Fit the global normalization used by DopplerView without loading a volume."""

    input_min, input_max = _volume_range(volume, frame_chunk_size)
    normalize_input = input_min < 0.0 or input_max > 1.0
    height, width = (int(volume.shape[-2]), int(volume.shape[-1]))
    y_slice, x_slice = _border_slices(height, width, border_amount)
    source_sum = np.float64(0.0)
    corrected_sum = np.float64(0.0)

    for start, stop in _frame_ranges(int(volume.shape[0]), frame_chunk_size):
        source = _normalized_chunk(
            volume[start:stop],
            input_min,
            input_max,
            normalize_input,
        )
        ratio = _spatial_flat_field_ratio(source, gaussian_width)
        source_sum += np.sum(source[:, y_slice, x_slice], dtype=np.float64)
        corrected_sum += np.sum(ratio[:, y_slice, x_slice], dtype=np.float64)

    output_scale = (
        float(source_sum / corrected_sum)
        if np.isfinite(corrected_sum) and corrected_sum != 0
        else 1.0
    )
    return FlatFieldParameters(
        gaussian_width=float(gaussian_width),
        border_amount=float(border_amount),
        input_min=float(input_min),
        input_max=float(input_max),
        normalize_input=bool(normalize_input),
        output_scale=output_scale,
    )


def corrected_flat_field_chunk(volume, frame_slice: slice, parameters) -> np.ndarray:
    source = _normalized_chunk(
        volume[frame_slice],
        parameters.input_min,
        parameters.input_max,
        parameters.normalize_input,
    )
    corrected = (
        np.float64(parameters.output_scale)
        * _spatial_flat_field_ratio(source, parameters.gaussian_width)
    )
    if parameters.normalize_input:
        corrected = parameters.input_min + (
            parameters.input_max - parameters.input_min
        ) * corrected
    return corrected.astype(np.float32)


def _volume_range(volume, frame_chunk_size: int) -> tuple[float, float]:
    minimum = np.inf
    maximum = -np.inf
    for start, stop in _frame_ranges(int(volume.shape[0]), frame_chunk_size):
        chunk = np.asarray(volume[start:stop])
        minimum = min(minimum, float(np.nanmin(chunk)))
        maximum = max(maximum, float(np.nanmax(chunk)))
    return minimum, maximum


def _normalized_chunk(
    values,
    input_min: float,
    input_max: float,
    normalize_input: bool,
) -> np.ndarray:
    source = np.asarray(values, dtype=np.float64)
    if not normalize_input:
        return source
    if input_max <= input_min:
        return np.zeros_like(source)
    return (source - input_min) / (input_max - input_min)


def _spatial_flat_field_ratio(values: np.ndarray, gaussian_width: float) -> np.ndarray:
    blurred = ndi.gaussian_filter(
        values,
        sigma=(0.0, float(gaussian_width), float(gaussian_width)),
        mode="reflect",
        truncate=2.0,
    )
    return values / (blurred + 1e-8)


def _border_slices(
    height: int,
    width: int,
    border_amount: float,
) -> tuple[slice, slice]:
    if border_amount == 0:
        return slice(0, height), slice(0, width)
    y0 = int(np.ceil(height * border_amount))
    y1 = int(np.floor(height * (1.0 - border_amount)))
    x0 = int(np.ceil(width * border_amount))
    x1 = int(np.floor(width * (1.0 - border_amount)))
    return slice(y0, y1), slice(x0, x1)


def _frame_ranges(frame_count: int, chunk_size: int):
    size = max(1, int(chunk_size))
    for start in range(0, frame_count, size):
        yield start, min(start + size, frame_count)
