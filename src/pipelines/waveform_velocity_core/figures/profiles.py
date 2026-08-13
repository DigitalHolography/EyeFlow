"""MATLAB-style diagnostics for transverse velocity profiles."""

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
    paths: list[Path] = []
    for vessel_name, color, segments in (
        ("artery", "#b2182b", getattr(context, "artery_segment_result", None)),
        ("vein", "#2166ac", getattr(context, "vein_segment_result", None)),
    ):
        if segments is None:
            continue
        raw = np.asarray(
            segments.transverse_velocity_profiles_masked,
            dtype=np.float32,
        )
        if raw.ndim == 4 and raw.shape[-1] > 0 and np.any(np.isfinite(raw)):
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
        if _has_centered_profiles(segments):
            centered = np.asarray(
                segments.centered_velocity_profiles,
                dtype=np.float32,
            )
            centered_x_um = np.asarray(
                segments.centered_profile_x_micrometers,
                dtype=np.float32,
            )
            aggregate = _hierarchical_profile_median(centered)
            aggregate_x_um = _nanmedian(centered_x_um, axis=(0, 1))
            paths.append(
                _save_mean_profile(
                    writer,
                    aggregate,
                    aggregate_x_um,
                    vessel_name,
                    color,
                )
            )
            gif_path = _save_profile_gif(
                writer,
                aggregate,
                aggregate_x_um,
                float(context.source_data.timing.dt_seconds),
                vessel_name,
                color,
                max_frames=max_gif_frames,
            )
            if gif_path is not None:
                paths.append(gif_path)
        paths.extend(_save_poiseuille_profiles(writer, segments, vessel_name))
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
    return _save_figure(
        writer,
        fig,
        f"{vessel_name}_velocity_profile_map.png",
    )


def _save_mean_profile(
    writer,
    aggregate: np.ndarray,
    x_um: np.ndarray,
    vessel_name: str,
    color: str,
) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    mean_profile = _nanmedian(aggregate, axis=0)
    ax.plot(x_um, mean_profile, color=color, linewidth=2.2)
    ax.axhline(0.0, color="0.25", linestyle="--", linewidth=0.8)
    y_min, y_max = _positive_focused_limits(mean_profile)
    ax.set(
        ylim=(y_min, y_max),
        xlabel="Centered transverse position (µm)",
        ylabel="Velocity (mm/s)",
        title=f"{vessel_name.capitalize()} cycle-mean centered profile",
    )
    ax.grid(alpha=0.18)
    return _save_figure(
        writer,
        fig,
        f"{vessel_name}_velocity_profile_cycle_mean.png",
    )


def _save_profile_gif(
    writer,
    aggregate: np.ndarray,
    x_um: np.ndarray,
    dt_seconds: float,
    vessel_name: str,
    color: str,
    *,
    max_frames: int,
) -> Path | None:
    from PIL import Image

    finite = aggregate[np.isfinite(aggregate)]
    if finite.size == 0:
        return None
    stride = max(int(np.ceil(aggregate.shape[0] / max(max_frames, 1))), 1)
    indexes = np.arange(0, aggregate.shape[0], stride, dtype=np.int32)
    y_min, y_max = _positive_focused_limits(finite)
    frames: list[Image.Image] = []
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.4, 3.8), dpi=90)
    (line,) = ax.plot(x_um, aggregate[indexes[0]], color=color, linewidth=2.2)
    ax.axhline(0.0, color="0.25", linestyle="--", linewidth=0.8)
    for frame_index in indexes:
        line.set_ydata(aggregate[frame_index])
        ax.set(
            xlim=(float(x_um[0]), float(x_um[-1])),
            ylim=(y_min, y_max),
            xlabel="Centered transverse position (µm)",
            ylabel="Velocity (mm/s)",
            title=(
                f"{vessel_name.capitalize()} profile — "
                f"t={frame_index * dt_seconds:.3f} s"
            ),
        )
        ax.grid(alpha=0.18)
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba(), dtype=np.uint8)
        frames.append(Image.fromarray(rgba[..., :3].copy()))
    plt.close(fig)
    path = _artifact_path(
        writer,
        f"{vessel_name}_velocity_profiles.gif",
        gif=True,
    )
    duration_ms = max(int(round(1000.0 * dt_seconds * stride)), 20)
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
    )
    return path


def _save_poiseuille_profiles(writer, segments, vessel_name: str) -> list[Path]:
    raw = np.asarray(
        segments.transverse_velocity_profiles_masked,
        dtype=np.float32,
    )
    x_um = np.asarray(segments.profile_x_micrometers, dtype=np.float32)
    coefficients = np.asarray(segments.poiseuille_coefficients, dtype=np.float32)
    origins = np.asarray(
        segments.poiseuille_origin_micrometers,
        dtype=np.float32,
    )
    roots = np.asarray(
        segments.poiseuille_roots_micrometers,
        dtype=np.float32,
    )
    spatial_std = np.asarray(
        segments.poiseuille_profile_spatial_std,
        dtype=np.float32,
    )
    branch_ids = np.asarray(segments.branch_ids, dtype=np.int32)
    paths: list[Path] = []
    for radius_index, branch_index in np.ndindex(raw.shape[:2]):
        coeff = coefficients[radius_index, branch_index]
        origin = float(origins[radius_index, branch_index])
        fit_roots = roots[radius_index, branch_index]
        if not (
            np.all(np.isfinite(coeff))
            and np.isfinite(origin)
            and np.all(np.isfinite(fit_roots))
        ):
            continue
        profile = _nanmean(raw[radius_index, branch_index], axis=0)
        profile_std = spatial_std[radius_index, branch_index]
        threshold = float(np.polyval(coeff, fit_roots[0]))
        selected = np.isfinite(profile) & (profile > threshold)
        plot_x_um = x_um - np.float32(origin)
        fit_x_relative = np.linspace(
            float(plot_x_um[0]),
            float(plot_x_um[-1]),
            256,
            dtype=np.float32,
        )
        fit_y = np.polyval(coeff, fit_x_relative)
        plt = _plt()
        fig, ax = plt.subplots(figsize=(7.4, 4.2))
        ax.fill_between(
            plot_x_um,
            profile - profile_std,
            profile + profile_std,
            color="0.7",
            edgecolor="none",
            label="spatial standard deviation",
        )
        ax.plot(plot_x_um[selected], profile[selected], "xk", label="fit samples")
        ax.plot(
            fit_x_relative,
            fit_y,
            color="black",
            linewidth=1.5,
            label="Poiseuille fit",
        )
        ax.axhline(0.0, color="0.25", linestyle="--", linewidth=0.8)
        ax.axhline(threshold, color="0.25", linestyle=":", linewidth=0.8)
        ax.plot(fit_roots, [-2.0, -2.0], color="black", linewidth=1.5)
        ax.plot(fit_roots, [-2.0, -2.0], "k|", markersize=10)
        plot_values = np.concatenate(
            (profile, profile - profile_std, profile + profile_std, fit_y, [threshold])
        )
        y_min, y_max = _positive_focused_limits(plot_values)
        ax.set(
            ylim=(y_min, y_max),
            xlabel="Position (µm)",
            ylabel="Velocity (mm/s)",
            title=(
                f"{vessel_name.capitalize()} branch {int(branch_ids[branch_index])}, "
                f"circle {radius_index + 1}"
            ),
        )
        ax.legend(frameon=False)
        ax.grid(alpha=0.18)
        prefix = "A" if vessel_name == "artery" else "V"
        suffix = (
            f"{prefix}{int(branch_ids[branch_index])}_C{radius_index + 1}"
            "_poiseuille_profile.png"
        )
        paths.append(_save_figure(writer, fig, suffix))
    return paths


def _has_centered_profiles(segments) -> bool:
    profiles = np.asarray(
        getattr(segments, "centered_velocity_profiles", np.empty((0,))),
    )
    return profiles.ndim == 4 and profiles.shape[-1] > 0 and np.any(np.isfinite(profiles))


def _hierarchical_profile_median(profiles: np.ndarray) -> np.ndarray:
    return _nanmedian(_nanmedian(profiles, axis=1), axis=0)


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


def _artifact_path(writer, suffix: str, *, gif: bool = False) -> Path:
    from input_output.output_manager import OutputType

    output_type = OutputType.GIF if gif else OutputType.PNG
    filename = f"{PROFILE_FOLDER}/{writer.stem}_{suffix}"
    path = writer.output.path_for(output_type, filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _save_figure(writer, fig, suffix: str) -> Path:
    path = _artifact_path(writer, suffix)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    _plt().close(fig)
    return path


def _nanmean(values: np.ndarray, *, axis):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def _nanmedian(values: np.ndarray, *, axis):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(values, axis=axis)
