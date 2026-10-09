"""CPU/GPU spatial resampling and mask transforms for sampled segments."""

from __future__ import annotations

from time import perf_counter

import numpy as np
from scipy import ndimage as ndi

try:
    import cv2
except ImportError:
    cv2 = None

from calculations.compute_backend import optional_cupy_backend
from utils.logger import Logger

INTERPOLATED_SEGMENT_SIDE = 128


def interpolate_segments(
    segment_maps: np.ndarray,
    output_side_pixels: int = INTERPOLATED_SEGMENT_SIDE,
) -> np.ndarray:
    """Interpolate segment values to one square size while preserving NaNs.

    Only the final two axes are resized; every preceding axis is retained.
    """

    values = np.asarray(segment_maps, dtype=np.float32)
    _assert_spatial_array(values, output_side_pixels)
    if values.shape[-2:] == (output_side_pixels, output_side_pixels):
        return values.copy()
    if values.shape[-2] == 0 or values.shape[-1] == 0:
        return np.full(
            (*values.shape[:-2], output_side_pixels, output_side_pixels),
            np.nan,
            dtype=np.float32,
        )
    return _interpolate_values(values, output_side_pixels)


def resample_rotate_segment(
    segment_maps: np.ndarray,
    rotation_degrees: float,
    output_side_pixels: int = INTERPOLATED_SEGMENT_SIDE,
    *,
    return_device: bool = False,
):
    """Resize and rotate one segment stack with a single affine resampling."""

    values = np.asarray(segment_maps, dtype=np.float32)
    _assert_spatial_array(values, output_side_pixels)
    if values.shape[-2] != values.shape[-1]:
        raise ValueError("native segment arrays must be square.")
    canvas_side = _rotation_canvas_side(output_side_pixels)
    output_shape = (*values.shape[:-2], canvas_side, canvas_side)
    if not np.isfinite(rotation_degrees):
        return np.full(output_shape, np.nan, dtype=np.float32)
    if values.shape[-1] == 0:
        return np.full(output_shape, np.nan, dtype=np.float32)
    return _resample_rotate_values(
        values,
        float(rotation_degrees),
        output_side_pixels,
        canvas_side,
        return_device=return_device,
    )


def interpolate_segment_masks(
    segment_masks: np.ndarray,
    output_side_pixels: int = INTERPOLATED_SEGMENT_SIDE,
) -> np.ndarray:
    """Interpolate Boolean segment masks with nearest-neighbor sampling."""

    masks = np.asarray(segment_masks, dtype=bool)
    _assert_spatial_array(masks, output_side_pixels)
    if masks.shape[-2:] == (output_side_pixels, output_side_pixels):
        return masks.copy()
    output_shape = (*masks.shape[:-2], output_side_pixels, output_side_pixels)
    if masks.shape[-2] == 0 or masks.shape[-1] == 0:
        return np.zeros(output_shape, dtype=bool)

    zoom = _spatial_zoom(masks, output_side_pixels)
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            resized = backend.ndimage.zoom(
                backend.cupy.asarray(masks, dtype=backend.cupy.float32),
                zoom,
                order=0,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            )
            return backend.cupy.asnumpy(resized) >= np.float32(0.5)
        except Exception:
            pass
    return ndi.zoom(
        masks.astype(np.float32),
        zoom,
        order=0,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ) >= np.float32(0.5)


def dilate_segment_masks(
    segment_masks: np.ndarray,
    *,
    iterations: int,
    exclusion_masks: np.ndarray | None = None,
    horizontal_only: bool = False,
) -> np.ndarray:
    """Expand masks, optionally along x only, excluding competing fringe."""

    masks = np.asarray(segment_masks, dtype=bool)
    exclusions = None
    if exclusion_masks is not None:
        exclusions = np.asarray(exclusion_masks, dtype=bool)
        if exclusions.shape != masks.shape:
            raise ValueError("exclusion_masks must match segment_masks shape.")
    if masks.ndim < 2:
        raise ValueError("segment_masks must end with spatial (y, x) axes.")
    if iterations < 0:
        raise ValueError("iterations must be non-negative.")
    if iterations == 0:
        return masks.copy()
    if cv2 is None:
        raise RuntimeError("OpenCV is required to dilate segment masks.")

    flat_masks = masks.reshape((-1, *masks.shape[-2:]))
    dilated = np.empty_like(flat_masks)
    kernel = np.ones((1 if horizontal_only else 3, 3), dtype=np.uint8)
    for mask_index, mask in enumerate(flat_masks):
        dilated[mask_index] = cv2.dilate(
            mask.astype(np.uint8),
            kernel,
            iterations=iterations,
        ).astype(bool)
    result = dilated.reshape(masks.shape)
    if exclusions is not None:
        result = masks | (result & ~exclusions)
    return result


def rotate_segments(
    interpolated_maps: np.ndarray,
    rotation_degrees: np.ndarray,
) -> np.ndarray:
    """Rotate uniformly sized segment maps on a diagonal-sized canvas.

    ``interpolated_maps`` must use ``(annulus, branch, ..., y, x)`` ordering.
    The rotation array must have matching ``(annulus, branch)`` dimensions.
    """

    values = np.asarray(interpolated_maps, dtype=np.float32)
    _assert_segment_array(values, rotation_degrees)
    canvas_side = _rotation_canvas_side(values.shape[-1])
    rotated = np.full(
        (*values.shape[:-2], canvas_side, canvas_side),
        np.nan,
        dtype=np.float32,
    )
    valid_indexes = np.argwhere(np.isfinite(rotation_degrees))
    progress_step = max(1, len(valid_indexes) // 10)
    started = perf_counter()
    Logger.log(
        f"Starting rotation of {len(valid_indexes)} segment stacks onto "
        f"{canvas_side}x{canvas_side} canvases."
    )
    for work_index, (annulus_index, branch_index) in enumerate(valid_indexes, start=1):
        segment = _pad_for_rotation(values[annulus_index, branch_index], np.nan)
        rotated[annulus_index, branch_index] = _rotate_values(
            segment,
            float(rotation_degrees[annulus_index, branch_index]),
        )
        if work_index % progress_step == 0 or work_index == len(valid_indexes):
            Logger.log(
                f"Segment rotation progress: {work_index}/{len(valid_indexes)} "
                f"in {perf_counter() - started:.2f}s."
            )
    return rotated


def rotate_segment_masks(
    interpolated_masks: np.ndarray,
    rotation_degrees: np.ndarray,
) -> np.ndarray:
    """Rotate uniformly sized Boolean masks on a diagonal-sized canvas."""

    masks = np.asarray(interpolated_masks, dtype=bool)
    _assert_segment_array(masks, rotation_degrees)
    if masks.ndim != 4:
        raise ValueError("interpolated_masks must have (annulus, branch, y, x) shape.")
    canvas_side = _rotation_canvas_side(masks.shape[-1])
    rotated = np.zeros((*masks.shape[:-2], canvas_side, canvas_side), dtype=bool)
    for annulus_index, branch_index in np.argwhere(np.isfinite(rotation_degrees)):
        segment = _pad_for_rotation(masks[annulus_index, branch_index], False)
        rotated[annulus_index, branch_index] = ndi.rotate(
            segment.astype(np.float32),
            float(rotation_degrees[annulus_index, branch_index]),
            reshape=False,
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        ) >= np.float32(0.5)
    return rotated


def _resample_rotate_values(
    values: np.ndarray,
    angle_degrees: float,
    output_side_pixels: int,
    canvas_side: int,
    *,
    return_device: bool = False,
):
    matrix, offset = _fused_affine_mapping(
        values.ndim,
        values.shape[-1],
        output_side_pixels,
        canvas_side,
        angle_degrees,
    )
    output_shape = (*values.shape[:-2], canvas_side, canvas_side)
    valid = np.isfinite(values)
    shared_validity = _shared_spatial_validity(valid)
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            gpu_values = backend.cupy.asarray(values)
            gpu_valid = backend.cupy.asarray(
                shared_validity if shared_validity is not None else valid
            )
            resampled_values = backend.ndimage.affine_transform(
                backend.cupy.where(
                    gpu_valid,
                    gpu_values,
                    backend.cupy.float32(0.0),
                ),
                backend.cupy.asarray(matrix),
                backend.cupy.asarray(offset),
                output_shape=output_shape,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
            )
            weight_matrix = matrix[-2:, -2:] if shared_validity is not None else matrix
            weight_offset = offset[-2:] if shared_validity is not None else offset
            weight_shape = (
                (canvas_side, canvas_side)
                if shared_validity is not None
                else output_shape
            )
            resampled_weights = backend.ndimage.affine_transform(
                gpu_valid.astype(backend.cupy.float32),
                backend.cupy.asarray(weight_matrix),
                backend.cupy.asarray(weight_offset),
                output_shape=weight_shape,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
            )
            finite_output = resampled_weights >= backend.cupy.float32(0.5)
            backend.cupy.divide(
                resampled_values,
                resampled_weights,
                out=resampled_values,
            )
            backend.cupy.copyto(
                resampled_values,
                backend.cupy.nan,
                where=~finite_output,
            )
            if return_device:
                return resampled_values
            return backend.cupy.asnumpy(resampled_values)
        except Exception as exc:
            _log_gpu_fallback("fused segment transform", exc)

    filled = np.where(valid, values, np.float32(0.0))
    resampled_values = ndi.affine_transform(
        filled,
        matrix,
        offset,
        output_shape=output_shape,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
    ).astype(np.float32, copy=False)
    weight_values = shared_validity if shared_validity is not None else valid
    weight_matrix = matrix[-2:, -2:] if shared_validity is not None else matrix
    weight_offset = offset[-2:] if shared_validity is not None else offset
    weight_shape = (
        (canvas_side, canvas_side)
        if shared_validity is not None
        else output_shape
    )
    resampled_weights = ndi.affine_transform(
        weight_values.astype(np.float32),
        weight_matrix,
        weight_offset,
        output_shape=weight_shape,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
    ).astype(np.float32, copy=False)
    finite_output = resampled_weights >= np.float32(0.5)
    np.divide(
        resampled_values,
        resampled_weights,
        out=resampled_values,
        where=finite_output,
    )
    np.copyto(resampled_values, np.nan, where=~finite_output)
    return resampled_values


def _fused_affine_mapping(
    ndim: int,
    source_side: int,
    interpolated_side: int,
    canvas_side: int,
    angle_degrees: float,
) -> tuple[np.ndarray, np.ndarray]:
    radians = np.deg2rad(angle_degrees)
    cosine = float(np.cos(radians))
    sine = float(np.sin(radians))
    rotation = np.asarray(((cosine, sine), (-sine, cosine)), dtype=np.float64)
    zoom = float(interpolated_side) / float(source_side)
    canvas_center = np.full(2, (canvas_side - 1.0) / 2.0, dtype=np.float64)
    padding_before = (canvas_side - interpolated_side) // 2
    spatial_offset = (
        canvas_center
        - rotation @ canvas_center
        - float(padding_before)
        + 0.5
    ) / zoom - 0.5
    matrix = np.eye(ndim, dtype=np.float64)
    matrix[-2:, -2:] = rotation / zoom
    offset = np.zeros(ndim, dtype=np.float64)
    offset[-2:] = spatial_offset
    return matrix, offset


def _shared_spatial_validity(valid: np.ndarray) -> np.ndarray | None:
    if valid.ndim == 2:
        return valid
    flattened = valid.reshape((-1, *valid.shape[-2:]))
    first = flattened[0]
    return first if np.all(flattened == first) else None


def _interpolate_values(values: np.ndarray, output_side_pixels: int) -> np.ndarray:
    zoom = _spatial_zoom(values, output_side_pixels)
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            gpu_values = backend.cupy.asarray(values)
            valid = backend.cupy.isfinite(gpu_values)
            resized_values = backend.ndimage.zoom(
                backend.cupy.where(valid, gpu_values, backend.cupy.float32(0.0)),
                zoom,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            )
            resized_weights = backend.ndimage.zoom(
                valid.astype(backend.cupy.float32),
                zoom,
                order=1,
                mode="grid-constant",
                cval=0.0,
                prefilter=False,
                grid_mode=True,
            )
            valid_output = resized_weights > backend.cupy.float32(1e-6)
            safe_weights = backend.cupy.where(
                valid_output,
                resized_weights,
                backend.cupy.float32(1.0),
            )
            resized = backend.cupy.empty(
                resized_values.shape,
                dtype=backend.cupy.float32,
            )
            backend.cupy.divide(
                resized_values,
                safe_weights,
                out=resized,
            )
            resized[~valid_output] = backend.cupy.nan
            return backend.cupy.asnumpy(resized)
        except Exception:
            pass

    valid = np.isfinite(values)
    resized_values = ndi.zoom(
        np.where(valid, values, np.float32(0.0)),
        zoom,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ).astype(np.float32, copy=False)
    resized_weights = ndi.zoom(
        valid.astype(np.float32),
        zoom,
        order=1,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ).astype(np.float32, copy=False)
    resized = np.full(resized_values.shape, np.nan, dtype=np.float32)
    np.divide(
        resized_values,
        resized_weights,
        out=resized,
        where=resized_weights > np.float32(1e-6),
    )
    return resized


def _rotate_values(values: np.ndarray, angle_degrees: float) -> np.ndarray:
    backend = optional_cupy_backend()
    if backend is not None:
        try:
            gpu_values = backend.cupy.asarray(values, dtype=backend.cupy.float32)
            valid = backend.cupy.isfinite(gpu_values)
            rotated_values = backend.ndimage.rotate(
                backend.cupy.where(valid, gpu_values, backend.cupy.float32(0.0)),
                angle_degrees,
                axes=(-2, -1),
                reshape=False,
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=False,
            )
            rotated_weights = backend.ndimage.rotate(
                valid.astype(backend.cupy.float32),
                angle_degrees,
                axes=(-2, -1),
                reshape=False,
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=False,
            )
            valid_output = rotated_weights >= backend.cupy.float32(0.5)
            safe_weights = backend.cupy.where(
                valid_output,
                rotated_weights,
                backend.cupy.float32(1.0),
            )
            rotated = backend.cupy.empty(
                rotated_values.shape,
                dtype=backend.cupy.float32,
            )
            backend.cupy.divide(
                rotated_values,
                safe_weights,
                out=rotated,
            )
            rotated[~valid_output] = backend.cupy.nan
            return backend.cupy.asnumpy(rotated)
        except Exception:
            pass

    valid = np.isfinite(values)
    rotated_values = ndi.rotate(
        np.where(valid, values, np.float32(0.0)),
        angle_degrees,
        axes=(-2, -1),
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    )
    rotated_weights = ndi.rotate(
        valid.astype(np.float32),
        angle_degrees,
        axes=(-2, -1),
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    )
    rotated = np.full(rotated_values.shape, np.nan, dtype=np.float32)
    np.divide(
        rotated_values,
        rotated_weights,
        out=rotated,
        where=rotated_weights >= np.float32(0.5),
    )
    return rotated.astype(np.float32, copy=False)


def _spatial_zoom(values: np.ndarray, output_side_pixels: int) -> tuple[float, ...]:
    return (1.0,) * (values.ndim - 2) + (
        output_side_pixels / values.shape[-2],
        output_side_pixels / values.shape[-1],
    )


def _rotation_canvas_side(interpolated_side_pixels: int) -> int:
    return int(np.ceil(np.sqrt(2.0) * (interpolated_side_pixels - 1))) + 1


def _pad_for_rotation(
    values: np.ndarray,
    fill_value: float | bool,
) -> np.ndarray:
    canvas_side = _rotation_canvas_side(values.shape[-1])
    total_padding = canvas_side - values.shape[-1]
    padding_before = total_padding // 2
    padding_after = total_padding - padding_before
    padding = [(0, 0)] * values.ndim
    padding[-2] = (padding_before, padding_after)
    padding[-1] = (padding_before, padding_after)
    return np.pad(values, padding, mode="constant", constant_values=fill_value)


def _log_gpu_fallback(operation: str, exc: Exception) -> None:
    Logger.log_debug(
        f"CuPy {operation} failed; using CPU fallback: "
        f"{type(exc).__name__}: {exc}"
    )


def _assert_spatial_array(values: np.ndarray, output_side_pixels: int) -> None:
    if values.ndim < 2:
        raise ValueError("segment arrays must end with spatial (y, x) axes.")
    if output_side_pixels <= 0:
        raise ValueError("output_side_pixels must be positive.")


def _assert_segment_array(
    values: np.ndarray,
    rotation_degrees: np.ndarray,
) -> None:
    if values.ndim < 4:
        raise ValueError(
            "segment arrays must have (annulus, branch, ..., y, x) shape."
        )
    if values.shape[-2] != values.shape[-1]:
        raise ValueError("interpolated segment arrays must be square.")
    if tuple(rotation_degrees.shape) != tuple(values.shape[:2]):
        raise ValueError(
            "rotation_degrees must match the segment annulus and branch axes."
        )
