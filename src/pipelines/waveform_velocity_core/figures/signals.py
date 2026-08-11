"""Signal time-series PNG exporters for waveform velocity analysis."""

from __future__ import annotations

from pathlib import Path

from .signal_inputs import (
    masked_video_signal as _masked_signal,
)
from input_output.writers.png import FigureArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _log,
    _vector,
    display_frequency as _display_frequency,
    display_velocity as _display_velocity,
)
from .plotting import _line_plot


def _export_signal_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    f_video = ctx.analysis.get("fRMS")
    f_bkg = ctx.analysis.get("fRMS_bkg")
    delta = ctx.analysis.get("deltafRMS")
    f_artery = _summary_or_masked_signal(
        ctx,
        "retinal_artery_fRMS_signal",
        f_video,
        ctx.artery_section_mask,
    )
    f_artery_bkg = _summary_or_masked_signal(
        ctx,
        "retinal_artery_fRMS_bkg_signal",
        f_bkg,
        ctx.artery_section_mask,
    )
    f_vein = _summary_or_masked_signal(
        ctx,
        "retinal_vein_fRMS_signal",
        f_video,
        ctx.vein_section_mask,
    )
    f_vein_bkg = _summary_or_masked_signal(
        ctx,
        "retinal_vein_fRMS_bkg_signal",
        f_bkg,
        ctx.vein_section_mask,
    )
    f_vessel_bkg = _summary_or_masked_signal(
        ctx,
        "retinal_vessel_fRMS_bkg_signal",
        f_bkg,
        ctx.vessel_section_mask,
    )
    if all(
        value is not None
        for value in (f_artery, f_artery_bkg, f_vein, f_vein_bkg, f_vessel_bkg)
    ):
        f_artery = _display_frequency(f_artery)
        f_artery_bkg = _display_frequency(f_artery_bkg)
        f_vein = _display_frequency(f_vein)
        f_vein_bkg = _display_frequency(f_vein_bkg)
        f_vessel_bkg = _display_frequency(f_vessel_bkg)
        paths.append(
            _line_plot(
                writer,
                "f_artery_graph.png",
                ctx.time,
                [(f_artery, "-", "tab:red", "arteries"), (f_artery_bkg, "--", "k", "background")],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
        paths.append(
            _line_plot(
                writer,
                "f_vein_graph.png",
                ctx.time,
                [(f_vein, "-", "tab:blue", "veins"), (f_vein_bkg, "--", "k", "background")],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
        paths.append(
            _line_plot(
                writer,
                "f_vessel_graph.png",
                ctx.time,
                [
                    (f_artery, "-", "tab:red", "arteries"),
                    (f_vein, "-", "tab:blue", "veins"),
                    (f_vessel_bkg, "--", "k", "background"),
                ],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
    else:
        _log(ctx, "Skipping fRMS signal PNGs; fRMS intermediates are unavailable.")

    delta_artery = _summary_or_masked_signal(
        ctx,
        "retinal_artery_deltafRMS_signal",
        delta,
        ctx.artery_section_mask,
    )
    delta_vein = _summary_or_masked_signal(
        ctx,
        "retinal_vein_deltafRMS_signal",
        delta,
        ctx.vein_section_mask,
    )
    if delta_artery is not None and delta_vein is not None:
        paths.append(
            _line_plot(
                writer,
                "df_vessel_graph.png",
                ctx.time,
                [
                    (
                        _display_frequency(delta_artery),
                        "-",
                        "tab:red",
                        "arteries",
                    ),
                    (
                        _display_frequency(delta_vein),
                        "-",
                        "tab:blue",
                        "veins",
                    ),
                ],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
    else:
        _log(ctx, "Skipping df_vessel_graph.png; deltafRMS is unavailable.")

    paths.append(
        _line_plot(
            writer,
            "v_vessel_graph.png",
            ctx.time,
            [
                (
                    _display_velocity(
                        _vector(ctx.analysis["retinal_artery_velocity_signal_filtered"])
                    ),
                    "-",
                    "tab:red",
                    "arteries",
                ),
                (
                    _display_velocity(
                        _vector(ctx.analysis["retinal_vein_velocity_signal_filtered"])
                    ),
                    "-",
                    "tab:blue",
                    "veins",
                ),
            ],
            title="average velocity in arteries and veins",
            xlabel="Time(s)",
            ylabel="Velocity (mm/s)",
        )
    )
    return paths


def _summary_or_masked_signal(ctx, key: str, video, mask):
    summary = ctx.analysis.get(key)
    if summary is not None:
        return _vector(summary)
    if video is None:
        return None
    return _masked_signal(video, mask)
