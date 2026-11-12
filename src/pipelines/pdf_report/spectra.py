"""Select and present the shared spectra/gaps PDF data."""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from input_output.readers.nodes import is_dataset
from input_output.schema import EyeFlowOutputPaths
from matplotlib.figure import Figure


@dataclass(frozen=True)
class SpectraGapsReportData:
    se_curve: np.ndarray | None = None
    sei: float | None = None
    sei_early: float | None = None
    sei_late: float | None = None
    acquisition_valid_fraction: float | None = None
    valid_se_phase_fraction: float | None = None
    retained_field_fraction: float | None = None


def extract_spectra_gaps_report_data(h5file) -> SpectraGapsReportData | None:
    """Read the once-per-run shared gap-spectra summary."""
    schema = EyeFlowOutputPaths.active()
    root = schema.spectrograms_masked_metrics_root
    qc = schema.spectrograms_masked_qc
    se_dataset = _dataset(h5file, f"{root}/SE")
    se_curve = None if se_dataset is None else np.asarray(se_dataset, dtype=float)
    if se_curve is not None and se_curve.ndim != 1:
        se_curve = None
    data = SpectraGapsReportData(
        se_curve=se_curve,
        sei=_scalar(h5file, f"{root}/SEI"),
        sei_early=_scalar(h5file, f"{root}/SEI_early"),
        sei_late=_scalar(h5file, f"{root}/SEI_late"),
        acquisition_valid_fraction=_scalar(
            h5file, f"{qc.data_completeness}/AcquisitionValidFraction",
        ),
        valid_se_phase_fraction=_scalar(
            h5file, f"{qc.data_completeness}/ValidSEPhaseFraction",
        ),
        retained_field_fraction=_scalar(
            h5file, f"{qc.spatial_support_quality}/RetainedFieldFraction",
        ),
    )
    values = (
        data.sei, data.sei_early, data.sei_late,
        data.acquisition_valid_fraction, data.valid_se_phase_fraction,
        data.retained_field_fraction,
    )
    return None if data.se_curve is None and all(value is None for value in values) else data


def build_spectra_gaps_pages(
    data: SpectraGapsReportData, folder_name: str,
) -> list[Figure]:
    """Build the standalone spectra and gaps report."""
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Spectra and Gaps")
    fig.text(
        0.08, 0.84, "Low-Frequency Gap Spectra",
        fontsize=13, fontweight="bold", va="top",
    )
    text_ax = fig.add_axes((0.08, 0.43, 0.38, 0.36))
    plot_ax = fig.add_axes((0.54, 0.43, 0.38, 0.36))
    text_ax.axis("off")
    _draw_spectral_section(text_ax, plot_ax, data)
    return [fig]


def _draw_spectral_section(text_ax, plot_ax, data: SpectraGapsReportData) -> None:
    def number(value, suffix=""):
        return "unavailable" if value is None else f"{value:.2f}{suffix}"

    lines = [
        f"Median spectral excess (SEI): {number(data.sei)}",
        f"Early / late SEI: {number(data.sei_early)} / {number(data.sei_late)}",
        "",
        "Descriptive quality indicators (no pass/fail threshold):",
        f"Finite acquisition spectra: {number(_percent(data.acquisition_valid_fraction), '%')}",
        f"Valid cardiac phases for SE: {number(_percent(data.valid_se_phase_fraction), '%')}",
        f"Retained gap-mask field: {number(_percent(data.retained_field_fraction), '%')}",
        "",
        "SEI is the median of the phase-resolved SE curve.",
        "Mask coverage is not an anatomical quality score.",
    ]
    text_ax.text(
        0, 0.95, "\n".join(lines), va="top", fontsize=9,
        linespacing=1.5, transform=text_ax.transAxes,
    )
    if data.se_curve is None or not data.se_curve.size:
        plot_ax.text(
            0.5, 0.5, "Phase-resolved SE unavailable",
            ha="center", va="center", transform=plot_ax.transAxes,
        )
    else:
        phase = np.arange(data.se_curve.size) / data.se_curve.size
        plot_ax.plot(phase, data.se_curve, color="black", linewidth=1.2)
        if data.sei is not None:
            plot_ax.axhline(
                data.sei, color="tab:blue", linestyle="--",
                linewidth=0.9, label="Median SEI",
            )
            plot_ax.legend(fontsize=8)
    plot_ax.set_xlim(0, 1)
    plot_ax.set_xlabel("Cardiac phase")
    plot_ax.set_ylabel("Spectral excess SE (unitless)")


def _dataset(h5file, path):
    node = h5file.get(path)
    return node if is_dataset(node) and node.size else None


def _scalar(h5file, path: str) -> float | None:
    dataset = _dataset(h5file, path)
    if dataset is None or dataset.size != 1:
        return None
    value = float(np.asarray(dataset, dtype=float).reshape(-1)[0])
    return value if np.isfinite(value) else None


def _percent(fraction: float | None) -> float | None:
    return None if fraction is None else 100 * fraction


def _add_page_header(fig: Figure, folder_name: str, section: str) -> None:
    fig.text(
        0.5, 0.978, folder_name, fontsize=10, color="0.35",
        ha="center", va="top",
    )
    fig.text(
        0.5, 0.952, section, fontsize=16, fontweight="bold",
        ha="center", va="top",
    )
