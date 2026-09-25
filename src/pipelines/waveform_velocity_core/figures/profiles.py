"""Diagnostics for raw transverse velocity profiles."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

from .common import _plt

PROFILE_FOLDER = "velocityProfiles"


def export_cross_section_profile_artifacts(
    writer,
    context,
    *,
    max_gif_frames: int = 100,
) -> list[Path]:
    """Export raw profile maps without obsolete centering or Poiseuille fits."""

    del max_gif_frames
    paths: list[Path] = []
    for vessel_name, segments in (
        ("artery", getattr(context, "artery_segment_result", None)),
        ("vein", getattr(context, "vein_segment_result", None)),
    ):
        if segments is None:
            continue
        raw = np.asarray(
            segments.transverse_velocity_profiles_masked,
            dtype=np.float32,
        )
        if raw.ndim != 4 or raw.shape[-1] == 0 or not np.any(np.isfinite(raw)):
            continue
        raw_aggregate = _hierarchical_profile_median(raw)
        raw_x_pixels = (
            np.arange(raw.shape[-1], dtype=np.float32)
            - np.float32((raw.shape[-1] - 1) / 2.0)
        )
        paths.append(
            _save_profile_map(
                writer,
                raw_aggregate,
                raw_x_pixels,
                float(context.source_data.timing.dt_seconds),
                vessel_name,
            )
        )
    return paths


def _save_profile_map(
    writer,
    aggregate: np.ndarray,
    x_pixels: np.ndarray,
    dt_seconds: float,
    vessel_name: str,
) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    finite = aggregate[np.isfinite(aggregate)]
    color_min, color_max = _positive_focused_limits(finite)
    extent = (
        float(x_pixels[0]),
        float(x_pixels[-1]),
        max((aggregate.shape[0] - 1) * dt_seconds, dt_seconds),
        0.0,
    )
    image = ax.imshow(
        aggregate,
        aspect="auto",
        interpolation="nearest",
        extent=extent,
        cmap="coolwarm",
        vmin=color_min,
        vmax=color_max,
    )
    fig.colorbar(image, ax=ax, label="Velocity (mm/s)")
    ax.set(
        xlabel="Transverse sample (pixels)",
        ylabel="Time (s)",
        title=f"{vessel_name.capitalize()} transverse velocity profiles (raw)",
    )
    return _save_figure(writer, fig, f"{vessel_name}_velocity_profile_map.png")


def _hierarchical_profile_median(profiles: np.ndarray) -> np.ndarray:
    values = np.asarray(profiles, dtype=np.float32)
    if values.ndim != 4:
        raise ValueError(
            "profile aggregation requires (radius, branch, frame, sample) values."
        )
    if values.shape[0] == 0:
        return np.full(values.shape[2:], np.nan, dtype=np.float32)
    per_radius = np.empty((values.shape[0], *values.shape[2:]), dtype=np.float32)
    for radius_index in range(values.shape[0]):
        per_radius[radius_index] = _nanmedian(values[radius_index], axis=0)
    return _nanmedian(per_radius, axis=0)


def _positive_focused_limits(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return -1.0, 1.0
    positive = finite[finite > 0]
    if not positive.size:
        limit = max(np.max(np.abs(finite)) * 1.12, 1.0)
        return -limit, limit
    upper = max(np.max(positive) * 1.12, np.finfo(float).eps)
    return min(0.0, max(np.min(finite), -0.25 * upper)), upper


def _artifact_path(writer, suffix: str) -> Path:
    from input_output.output_manager import OutputType

    filename = f"{PROFILE_FOLDER}/{writer.stem}_{suffix}"
    path = writer.output.path_for(OutputType.PNG, filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _save_figure(writer, fig, suffix: str) -> Path:
    path = _artifact_path(writer, suffix)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    _plt().close(fig)
    return path


def _nanmedian(values: np.ndarray, *, axis):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        try:
            return np.nanmedian(values, axis=axis)
        except IndexError:
            # Some NumPy releases fail in the partition implementation for
            # very sparse arrays. The finite-only fallback is slower but this
            # path is diagnostic-only and keeps raw-profile export reliable.
            return _finite_median(values, axis=axis)


def _finite_median(values: np.ndarray, *, axis):
    values = np.asarray(values)
    axis = int(axis)
    if axis < 0:
        axis += values.ndim
    if axis < 0 or axis >= values.ndim:
        raise np.AxisError(axis, values.ndim)
    moved = np.moveaxis(values, axis, -1)
    result = np.full(moved.shape[:-1], np.nan, dtype=np.float32)
    for index in np.ndindex(result.shape):
        finite = moved[index]
        finite = finite[np.isfinite(finite)]
        if finite.size:
            result[index] = np.median(finite)
    return result
