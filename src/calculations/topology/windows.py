"""Native window geometry shared by topology construction and sampling."""

from __future__ import annotations


def centered_window_bounds(
    image_shape: tuple[int, int],
    center_xy: tuple[int, int],
    side_pixels: int,
) -> tuple[int, int, int, int]:
    """Clip an odd, centered square to native image bounds (x0, x1, y0, y1)."""

    assert side_pixels > 0 and side_pixels % 2 == 1
    half_width = side_pixels // 2
    center_x, center_y = center_xy
    return (
        max(center_x - half_width, 0),
        min(center_x + half_width + 1, int(image_shape[1])),
        max(center_y - half_width, 0),
        min(center_y + half_width + 1, int(image_shape[0])),
    )


def window_target_slices(
    bounds_xyxy,
    center_xy,
    side_pixels: int,
) -> tuple[slice, slice]:
    """Place clipped native bounds inside their conceptual padded square."""

    x_start, x_stop, y_start, y_stop = (int(value) for value in bounds_xyxy)
    center_x, center_y = (int(value) for value in center_xy)
    conceptual_x_start = center_x - side_pixels // 2
    conceptual_y_start = center_y - side_pixels // 2
    target_x_start = x_start - conceptual_x_start
    target_y_start = y_start - conceptual_y_start
    return (
        slice(target_y_start, target_y_start + y_stop - y_start),
        slice(target_x_start, target_x_start + x_stop - x_start),
    )
