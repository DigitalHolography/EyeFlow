"""Public PNG export entrypoint for waveform-shape velocity figures."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from input_output.writers.png import PngArtifactWriter

from .common import PulseFigureContext, _matplotlib, _output_stem, _section_mask
from .correlation import _export_correlation_plots
from .maps import _export_maps
from .signals import _export_signal_plots
from .spectral import _export_spectral_plots
from .systole import _export_systole_plots
from .velocity_maps import _export_final_visualizations
from .waveforms import _export_ri_pi_plots, _export_waveform_plots

if TYPE_CHECKING:
    from calculations.blood_flow_velocity import PerBeatAnalysisResult
    from pipelines.waveform_shape_metrics.runner import WaveformShapeMetricsContext


PULSE_PNG_SUFFIXES = (
    "f_artery_graph.png",
    "f_vein_graph.png",
    "f_vessel_graph.png",
    "df_vessel_graph.png",
    "v_vessel_graph.png",
    "f_bkg_map.png",
    "f_bkg_colorBar.png",
    "df_map_vessel.png",
    "df_colorBar_vessel.png",
    "df_map.png",
    "df_colorBar.png",
    "f_map.png",
    "f_colorBar.png",
    "artery_seg_map_bkg.png",
    "vein_seg_map_bkg.png",
    "map_df_vessel.png",
    "find_systoles_indices_artery.png",
    "find_systoles_indices_vein.png",
    "RI_v_artery.png",
    "PI_v_artery.png",
    "RI_v_vein.png",
    "PI_v_vein.png",
    "ArterialWaveformAnalysis_v_artery.png",
    "VenousWaveformAnalysis_v_vein.png",
    "ArterialSpectralAnalysis_v_artery.png",
    "ArterialSpectralAnalysis_Phase_v_artery.png",
    "VenousSpectralAnalysis_v_vein.png",
    "VenousSpectralAnalysis_Phase_v_vein.png",
    "syntheticSpectralAnalysis.png",
    "arterial_venous_msc.png",
    "arterial_venous_correlation.png",
    "detrended_signals.png",
    "lags.png",
    "Transfer_function_Velocity_AV_mod.png",
    "Transfer_function_Velocity_AV_phase.png",
    "ArterialVenous_Delay.png",
    "map_v_Artery.png",
    "map_v_Vein.png",
    "map_v_Vessel.png",
    "colorbar_v_Artery.png",
    "colorbar_v_Vein.png",
    "v_map.png",
    "histogramVelocityartery.png",
    "histogramVelocityvein.png",
    "AVGflowVideoCombined.png",
)

EXPORTERS = (
    _export_signal_plots,
    _export_maps,
    _export_systole_plots,
    _export_ri_pi_plots,
    _export_waveform_plots,
    _export_spectral_plots,
    _export_correlation_plots,
    _export_final_visualizations,
)


def export_pulse_pngs(
    output,
    context: WaveformShapeMetricsContext,
    per_beat_result: PerBeatAnalysisResult,
    *,
    log=None,
) -> list[str]:
    """Export core pulse-analysis PNGs for a waveform-shape run."""

    if not getattr(output, "available", False):
        return []
    _matplotlib()
    source_data = context.source_data
    analysis = context.dopplerview_analysis
    frame_count = int(np.asarray(analysis["retinal_artery_velocity_signal"]).size)
    pulse_context = PulseFigureContext(
        output=output,
        stem=_output_stem(output),
        time=np.arange(frame_count, dtype=np.float32) * np.float32(source_data.timing.dt_seconds),
        dt_seconds=float(source_data.timing.dt_seconds),
        moment0_avg=np.mean(source_data.moment0, axis=0, dtype=np.float32),
        artery_mask=np.asarray(source_data.retinal_artery_mask, dtype=bool),
        vein_mask=np.asarray(source_data.retinal_vein_mask, dtype=bool),
        section_mask=_section_mask(analysis, source_data.retinal_artery_mask.shape),
        analysis=analysis,
        per_beat_result=per_beat_result,
        log=log,
    )
    writer = PngArtifactWriter(output, pulse_context.stem)
    paths: list[Path] = []
    for exporter in EXPORTERS:
        paths.extend(exporter(writer, pulse_context))
    return [str(path) for path in paths]
