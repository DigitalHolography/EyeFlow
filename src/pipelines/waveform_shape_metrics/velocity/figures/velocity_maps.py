"""Velocity-map, histogram, and combined visualization PNG exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_flow_velocity.context_builders.spectrum import (
    histogram_matrix as _histogram_matrix,
)
from calculations.blood_flow_velocity.signal_analysis.waveform import rescale as _rescale
from input_output.writers.png import PngArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _array_or_none,
    _log,
    _plt,
    display_velocity as _display_velocity,
)
from .plotting import (
    _colorbar,
    _image_map,
    _matlab_vessel_colormap,
    _style_axes,
    _velocity_gradient_values,
)


def _export_final_visualizations(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    velocity_map = _array_or_none(ctx.analysis.get("retinal_vessel_velocity"))
    velocity_avg = _array_or_none(ctx.analysis.get("velocity_map_avg"))
    if velocity_avg is None and velocity_map is not None:
        velocity_avg = np.nanmean(velocity_map, axis=0)
    if velocity_avg is None:
        _log(ctx, "Skipping final velocity visualizations; velocity map is unavailable.")
        return []
    velocity_map_display = _display_velocity(velocity_map) if velocity_map is not None else None
    velocity_avg_display = _display_velocity(velocity_avg)
    vmax = _velocity_colorbar_vmax(velocity_avg_display, ctx.section_mask, velocity_map_display)
    paths = [
        _image_map(
            writer,
            np.where(ctx.artery_section_mask, velocity_avg_display, 0.0),
            "map_v_Artery.png",
            cmap=_matlab_vessel_colormap("artery"),
            colorbar=True,
            label="mm/s",
        ),
        _image_map(
            writer,
            np.where(ctx.vein_section_mask, velocity_avg_display, 0.0),
            "map_v_Vein.png",
            cmap=_matlab_vessel_colormap("vein"),
            colorbar=True,
            label="mm/s",
        ),
        _image_map(
            writer,
            np.where(ctx.vessel_section_mask, velocity_avg_display, 0.0),
            "map_v_Vessel.png",
            cmap="turbo",
            colorbar=True,
            label="mm/s",
        ),
        _colorbar(
            writer,
            "colorbar_v_Artery.png",
            cmap=_matlab_vessel_colormap("artery"),
            vmin=0.0,
            vmax=vmax,
            label="mm/s",
        ),
        _colorbar(
            writer,
            "colorbar_v_Vein.png",
            cmap=_matlab_vessel_colormap("vein"),
            vmin=0.0,
            vmax=vmax,
            label="mm/s",
        ),
    ]
    flow_rgb = _flow_rgb(ctx, velocity_avg_display, velocity_map_display)
    writer.save_image(flow_rgb, "v_map.png")
    paths.append(writer.path("v_map.png"))
    hist_artery = _histogram_plot(
        writer,
        ctx,
        velocity_map_display,
        ctx.artery_section_mask,
        "histogramVelocityartery.png",
        "artery",
    )
    hist_vein = _histogram_plot(
        writer,
        ctx,
        velocity_map_display,
        ctx.vein_section_mask,
        "histogramVelocityvein.png",
        "vein",
    )
    paths.extend([hist_artery, hist_vein])
    paths.append(_combined_plot(writer, ctx, flow_rgb, velocity_map_display))
    return paths

def _histogram_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    velocity_map: np.ndarray | None,
    mask: np.ndarray,
    suffix: str,
    vessel: str,
) -> Path:
    if velocity_map is None:
        velocity_map = _display_velocity(
            np.asarray(ctx.analysis["velocity_map_avg"], dtype=np.float32)
        )[None, :, :]
    histo = _histogram_matrix(velocity_map, mask)
    fig, ax = _plt().subplots(figsize=(6.0, 2.75))
    extent = [0, velocity_map.shape[0] * ctx.dt_seconds, histo.vmin, histo.vmax]
    ax.imshow(
        histo.counts,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap=_vessel_histogram_colormap(vessel),
        interpolation="nearest",
        vmin=0,
        vmax=histo.count_max,
    )
    ax.set_facecolor("black")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    _style_axes(ax)
    return writer.savefig(fig, suffix)

def _vessel_histogram_colormap(vessel: str):
    return _matlab_vessel_colormap(vessel)

def _velocity_colorbar_vmax(
    velocity_avg: np.ndarray,
    section_mask: np.ndarray,
    velocity_map: np.ndarray | None,
) -> float:
    del velocity_map
    section_values = np.asarray(velocity_avg, dtype=np.float32)[section_mask]
    finite = np.isfinite(section_values)
    if not np.any(finite):
        return 1.0
    vmax = float(np.nanmax(np.abs(section_values[finite])))
    return vmax if vmax > 0.0 else 1.0

def _combined_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    flow_rgb: np.ndarray,
    velocity_map: np.ndarray | None,
) -> Path:
    fig = _plt().figure(figsize=(8.8, 5.5), constrained_layout=True)
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.25, 1.0),
        height_ratios=(1.0, 1.0),
        hspace=0.18,
    )
    map_ax = fig.add_subplot(grid[:, 0])
    hist_axes = (fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1]))
    map_ax.imshow(flow_rgb)
    map_ax.axis("off")
    for ax, mask, title, vessel, xlabel in (
        (hist_axes[0], ctx.artery_section_mask, "Artery", "artery", ""),
        (hist_axes[1], ctx.vein_section_mask, "Vein", "vein", "Time (s)"),
    ):
        histogram_source = (
            velocity_map
            if velocity_map is not None
            else _display_velocity(np.asarray(ctx.analysis["velocity_map_avg"]))[None, :, :]
        )
        histo = _histogram_matrix(
            histogram_source,
            mask,
        )
        extent = [0, histo.counts.shape[1] * ctx.dt_seconds, histo.vmin, histo.vmax]
        ax.imshow(
            histo.counts,
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=_vessel_histogram_colormap(vessel),
            interpolation="nearest",
            vmin=0,
            vmax=histo.count_max,
        )
        ax.set_facecolor("black")
        ax.set_title(title, pad=3)
        ax.set_xlabel(xlabel, labelpad=2)
        ax.set_ylabel("Velocity (mm/s)", labelpad=2)
        _style_axes(ax)
    return writer.savefig(fig, "AVGflowVideoCombined.png", dpi=150)

def _flow_rgb(
    ctx: PulseFigureContext,
    velocity_avg: np.ndarray,
    velocity_map: np.ndarray | None = None,
) -> np.ndarray:
    background = _rescale(ctx.moment0_avg)
    rgb = np.dstack([background, background, background])
    value = _velocity_gradient_values(velocity_avg, ctx.section_mask, velocity_map)
    artery_color = _matlab_vessel_colormap("artery")(value)[..., :3]
    vein_color = _matlab_vessel_colormap("vein")(value)[..., :3]
    vessel_color = _matlab_vessel_colormap("av")(value)[..., :3]
    mask_av = ctx.artery_mask & ctx.vein_mask
    artery_only = ctx.artery_mask & ~mask_av
    vein_only = ctx.vein_mask & ~mask_av
    rgb[artery_only] = artery_color[artery_only]
    rgb[vein_only] = vein_color[vein_only]
    rgb[mask_av] = vessel_color[mask_av]
    return np.clip(rgb, 0, 1)
