"""Signal time-series PNG exporters for velocity analysis analysis."""

from __future__ import annotations

from pathlib import Path

from input_output.writers.png import FigureArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _vector,
    display_frequency as _display_frequency,
    display_velocity as _display_velocity,
)
from .plotting import _line_plot


def _export_signal_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    velocity = ctx.velocity
    f_artery = _display_frequency(_vector(velocity.artery.frms))
    f_artery_bkg = _display_frequency(_vector(velocity.artery.frms_background))
    f_vein = _display_frequency(_vector(velocity.vein.frms))
    f_vein_bkg = _display_frequency(_vector(velocity.vein.frms_background))
    f_vessel_bkg = _display_frequency(_vector(velocity.vessel_frms_background))
    paths.append(
        _line_plot(
            writer,
            "f_artery_graph.png",
            ctx.time,
            [
                (f_artery, "-", "tab:red", "arteries"),
                (f_artery_bkg, "--", "k", "background"),
            ],
            xlabel="Time(s)",
            ylabel="frequency (kHz)",
        )
    )
    paths.append(
        _line_plot(
            writer,
            "f_vein_graph.png",
            ctx.time,
            [
                (f_vein, "-", "tab:blue", "veins"),
                (f_vein_bkg, "--", "k", "background"),
            ],
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
    paths.append(
        _line_plot(
            writer,
            "df_vessel_graph.png",
            ctx.time,
            [
                (
                    _display_frequency(velocity.artery.delta_frms),
                    "-",
                    "tab:red",
                    "arteries",
                ),
                (
                    _display_frequency(velocity.vein.delta_frms),
                    "-",
                    "tab:blue",
                    "veins",
                ),
            ],
            xlabel="Time(s)",
            ylabel="frequency (kHz)",
        )
    )

    paths.append(
        _line_plot(
            writer,
            "v_vessel_graph.png",
            ctx.time,
            [
                (
                    _display_velocity(_vector(velocity.continuous("artery"))),
                    "-",
                    "tab:red",
                    "arteries",
                ),
                (
                    _display_velocity(_vector(velocity.continuous("vein"))),
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
