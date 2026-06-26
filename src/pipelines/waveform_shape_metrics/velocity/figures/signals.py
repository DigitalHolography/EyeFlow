"""Signal time-series PNG exporters for waveform-shape velocity analysis."""

from __future__ import annotations

from pathlib import Path

from calculations.blood_flow_velocity.context_builders.spectrum import (
    masked_video_signal as _masked_signal,
)
from input_output.writers.png import PngArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _array_or_none,
    _log,
    _vector,
    display_frequency as _display_frequency,
    display_velocity as _display_velocity,
)
from .plotting import _line_plot


def _export_signal_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    f_video = _array_or_none(ctx.analysis.get("fRMS"))
    f_bkg = _array_or_none(ctx.analysis.get("fRMS_bkg"))
    delta = _array_or_none(ctx.analysis.get("deltafRMS"))
    if f_video is not None and f_bkg is not None:
        f_artery = _display_frequency(_masked_signal(f_video, ctx.artery_section_mask))
        f_artery_bkg = _display_frequency(_masked_signal(f_bkg, ctx.artery_section_mask))
        f_vein = _display_frequency(_masked_signal(f_video, ctx.vein_section_mask))
        f_vein_bkg = _display_frequency(_masked_signal(f_bkg, ctx.vein_section_mask))
        f_vessel_bkg = _display_frequency(_masked_signal(f_bkg, ctx.vessel_section_mask))
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

    if delta is not None:
        paths.append(
            _line_plot(
                writer,
                "df_vessel_graph.png",
                ctx.time,
                [
                    (
                        _display_frequency(_masked_signal(delta, ctx.artery_section_mask)),
                        "-",
                        "tab:red",
                        "arteries",
                    ),
                    (
                        _display_frequency(_masked_signal(delta, ctx.vein_section_mask)),
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
                    _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"])),
                    "-",
                    "tab:red",
                    "arteries",
                ),
                (
                    _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"])),
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
