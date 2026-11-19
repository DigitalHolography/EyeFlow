"""Format-safe color composition helpers."""

from __future__ import annotations


def opaque_on_background(
    color, *, opacity: float, background="white",
) -> tuple[float, float, float]:
    """Return the opaque RGB result of compositing a color on a background."""
    from matplotlib.colors import to_rgb

    opacity = float(opacity)
    if not 0.0 <= opacity <= 1.0:
        raise ValueError("opacity must be between zero and one.")
    foreground_rgb = to_rgb(color)
    background_rgb = to_rgb(background)
    return tuple(
        opacity * foreground + (1.0 - opacity) * base
        for foreground, base in zip(foreground_rgb, background_rgb)
    )


__all__ = ["opaque_on_background"]
