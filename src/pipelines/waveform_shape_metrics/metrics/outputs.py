"""Output packing for global and by-segment waveform-shape metrics."""

from __future__ import annotations

import numpy as np

from pipeline_engine import with_attrs

from .calculator import WaveformShapeMetricsCalculator
from .formulas import LATEX_FORMULAS
from .models import WaveformShapeMetricInputs


def pack_waveform_shape_metric_calculations(
    inputs: WaveformShapeMetricInputs,
    *,
    include_segments: bool = True,
) -> dict[str, object]:
    """Pack the existing global and by-segment waveform output groups."""
    beat_periods = np.asarray(inputs.beat_period_seconds, dtype=np.float32)
    calculator = WaveformShapeMetricsCalculator()
    metrics: dict[str, object] = {}

    for vessel_name, vessel in (
        ("artery", inputs.artery),
        ("vein", inputs.vein),
    ):
        if (
            include_segments
            and vessel.raw_segments is not None
            and vessel.bandlimited_segments is not None
        ):
            _pack_segment_outputs(
                metrics,
                calculator,
                vessel_name,
                vessel.raw_segments,
                vessel.bandlimited_segments,
                beat_periods,
            )
        if vessel.raw_global is not None and vessel.bandlimited_global is not None:
            _pack_global_outputs(
                metrics,
                calculator,
                vessel_name,
                vessel.raw_global,
                vessel.bandlimited_global,
                beat_periods,
            )

    return metrics


def _pack_segment_outputs(
    metrics: dict[str, object],
    calculator: WaveformShapeMetricsCalculator,
    vessel_name: str,
    raw_segments: np.ndarray,
    bandlimited_segments: np.ndarray,
    beat_periods: np.ndarray,
) -> None:
    seg_band, branch_band, global_band, branch_count_band, radius_count_band, note_band = (
        calculator._compute_block_segment(bandlimited_segments, beat_periods)
    )
    seg_raw, branch_raw, global_raw, branch_count_raw, radius_count_raw, _note_raw = (
        calculator._compute_block_segment(raw_segments, beat_periods)
    )

    note = note_band
    if (branch_count_band != branch_count_raw) or (radius_count_band != radius_count_raw):
        note = note_band + " | WARNING: raw/band branch/radius dims differ."

    _pack_group(
        metrics,
        vessel_name,
        "by_segment/bandlimited_segment",
        seg_band,
        {
            "definition": ["per-segment metrics stored as (beat, branch, radius)"],
            "segment_indexing": [note],
        },
    )
    _pack_group(
        metrics,
        vessel_name,
        "by_segment/raw_segment",
        seg_raw,
        {
            "definition": ["per-segment metrics stored as (beat, branch, radius)"],
            "segment_indexing": [note],
        },
    )
    _pack_group(
        metrics,
        vessel_name,
        "by_segment/bandlimited_branch",
        branch_band,
        {"definition": ["median over radii per branch"]},
    )
    _pack_group(
        metrics,
        vessel_name,
        "by_segment/raw_branch",
        branch_raw,
        {"definition": ["median over radii per branch"]},
    )
    _pack_group(
        metrics,
        vessel_name,
        "by_segment/bandlimited_global",
        global_band,
        {"definition": ["median over all branch-radius segment values per beat"]},
    )
    _pack_group(
        metrics,
        vessel_name,
        "by_segment/raw_global",
        global_raw,
        {"definition": ["median over all branch-radius segment values per beat"]},
    )

    _pack_parameters(metrics, vessel_name, "by_segment", calculator)


def _pack_global_outputs(
    metrics: dict[str, object],
    calculator: WaveformShapeMetricsCalculator,
    vessel_name: str,
    raw_global: np.ndarray,
    bandlimited_global: np.ndarray,
    beat_periods: np.ndarray,
) -> None:
    raw = calculator._compute_block_global(raw_global, beat_periods)
    bandlimited = calculator._compute_block_global(bandlimited_global, beat_periods)
    metric_keys = calculator._metric_keys()

    for key, definition, unit, family in metric_keys:
        for signal_type, values in (("raw", raw), ("bandlimited", bandlimited)):
            metrics[f"{vessel_name}/global/{signal_type}/{key}"] = with_attrs(
                values[key],
                {
                    "unit": [unit],
                    "definition": [definition],
                    "metric_family": [family],
                    "latex_formula": [LATEX_FORMULAS[key]],
                },
            )

    _pack_parameters(metrics, vessel_name, "global", calculator)

    diagnostics = {
        "A2_cumsum",
        "A2_m",
        "A2_cumsum_interp",
        "A2_m_interp",
        "m_50",
        "m_80",
        "rho_h",
        "w_h",
        "delta_phi_all",
        "t_phi_n",
        "t_phi_n_over_T",
        "phase_harmonics_used",
        "harmonic_phases",
        "harmonic_weights",
        "harmonic_energies_weights",
    }
    raw_graphics = calculator._compute_graphics_support_block(
        raw_global,
        beat_periods,
    )
    bandlimited_graphics = calculator._compute_graphics_support_block(
        bandlimited_global,
        beat_periods,
    )
    for signal_type, graphics in (
        ("raw", raw_graphics),
        ("bandlimited", bandlimited_graphics),
    ):
        for name, values in graphics.items():
            suffix = f"global/{signal_type}/"
            if name in diagnostics:
                suffix += "diagnostics/"
            metrics[f"{vessel_name}/{suffix}{name}"] = values


def _pack_group(
    metrics: dict[str, object],
    vessel_name: str,
    group_name: str,
    values: dict[str, np.ndarray],
    attrs: dict[str, object],
) -> None:
    for metric_name, value in values.items():
        metrics[f"{vessel_name}/{group_name}/{metric_name}"] = with_attrs(
            value,
            attrs,
        )


def _pack_parameters(
    metrics: dict[str, object],
    vessel_name: str,
    group_name: str,
    calculator: WaveformShapeMetricsCalculator,
) -> None:
    parameters = {
        "ratio_R_VTI": calculator.ratio_R_VTI,
        "ratio_SF_VTI": calculator.ratio_SF_VTI,
        "ratio_vend_start": calculator.ratio_vend_start,
        "ratio_vend_end": calculator.ratio_vend_end,
        "eps": calculator.eps,
        "ratio_W50": calculator.ratio_W50,
        "ratio_W80": calculator.ratio_W80,
        "H_LOW_MAX": calculator.H_LOW_MAX,
        "H_MAX": calculator.H_MAX,
    }
    for name, value in parameters.items():
        metrics[f"{vessel_name}/{group_name}/params/{name}"] = np.asarray(
            value,
            dtype=int if name.startswith("H_") else float,
        )
