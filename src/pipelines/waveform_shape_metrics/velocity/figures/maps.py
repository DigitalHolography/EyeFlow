"""Frequency and delta-frequency map PNG exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from input_output.writers.png import PngArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _array_or_none,
    display_frequency as _display_frequency,
    display_velocity as _display_velocity,
)
from .plotting import (
    _heatmap_with_colorbar,
    _image_map,
    _matlab_vessel_colormap,
    _mask_background_map,
    _velocity_gradient_values,
)


def _export_maps(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    f_bkg_avg = _array_or_none(ctx.analysis.get("fRMS_bkg_avg"))
    f_avg = _array_or_none(ctx.analysis.get("fRMS_avg"))
    delta = _array_or_none(ctx.analysis.get("deltafRMS"))
    if f_bkg_avg is not None:
        f_bkg_avg = _display_frequency(f_bkg_avg)
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                f_bkg_avg,
                "f_bkg_map.png",
                "f_bkg_colorBar.png",
                cmap="gray",
                label="background RMS frequency (kHz)",
            )
        )
    if delta is not None:
        df_mean = _display_frequency(np.nanmean(delta, axis=0))
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                np.where(ctx.vessel_section_mask, df_mean, np.nan),
                "df_map_vessel.png",
                "df_colorBar_vessel.png",
                cmap="gray",
                label="Delta Doppler RMS frequency (kHz)",
            )
        )
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                df_mean,
                "df_map.png",
                "df_colorBar.png",
                cmap="gray",
                label="Delta Doppler RMS frequency (kHz)",
            )
        )
        paths.append(
            _image_map(
                writer,
                np.where(ctx.vessel_section_mask, df_mean, 0.0),
                "map_df_vessel.png",
                cmap="turbo",
                colorbar=True,
                label="kHz",
            )
        )
    if f_avg is not None:
        f_avg = _display_frequency(f_avg)
        velocity_map = _array_or_none(ctx.analysis.get("retinal_vessel_velocity"))
        velocity_avg = _array_or_none(ctx.analysis.get("velocity_map_avg"))
        if velocity_avg is None and velocity_map is not None:
            velocity_avg = np.nanmean(velocity_map, axis=0)
        velocity_values = None
        if velocity_avg is not None:
            velocity_values = _velocity_gradient_values(
                _display_velocity(velocity_avg),
                ctx.section_mask,
                _display_velocity(velocity_map) if velocity_map is not None else None,
            )
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                f_avg,
                "f_map.png",
                "f_colorBar.png",
                cmap="gray",
                label="RMS frequency (kHz)",
            )
        )
        paths.extend(
            [
                _mask_background_map(
                    writer,
                    f_avg,
                    ctx.artery_mask,
                    "artery_seg_map_bkg.png",
                    color="red",
                    cmap=_matlab_vessel_colormap("artery"),
                    values=velocity_values,
                ),
                _mask_background_map(
                    writer,
                    f_avg,
                    ctx.vein_mask,
                    "vein_seg_map_bkg.png",
                    color="blue",
                    cmap=_matlab_vessel_colormap("vein"),
                    values=velocity_values,
                ),
            ]
        )
    return paths
