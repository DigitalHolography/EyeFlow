"""Export branch lumen-size signals from masked transverse gradient metrics."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.math import nanmean, nanmedian
from input_output.output_manager import OutputType


def export_lumen_size_pngs(
    output,
    branch_signals: np.ndarray,
    branch_ids: np.ndarray,
    *,
    vessel_name: str,
    period_seconds: float,
) -> list[Path]:
    """Plot branch medians and the mean curve of the largest-quartile branches.

    Input axes are (time, branch), using the persisted median over beat and
    radius. Branch selection uses the temporal median of each curve, with
    ties at the 75th percentile included.
    Only NaNs are ignored; no lumen-size quality-control mask is applied.
    Time is shown in seconds across the mean cardiac-cycle duration.
    """
    branch_signals = np.asarray(branch_signals, dtype=np.float32)
    branch_ids = np.asarray(branch_ids).reshape(-1)
    if branch_signals.ndim != 2 or branch_signals.shape[1] != branch_ids.size:
        raise ValueError(
            "branch_signals must have shape (time, branch) with one branch ID per branch."
        )
    if not np.isfinite(period_seconds) or period_seconds <= 0.0:
        raise ValueError("period_seconds must be finite and positive.")

    branch_medians = nanmedian(branch_signals, axis=0)
    finite = np.isfinite(branch_medians)
    selected = np.zeros(branch_ids.size, dtype=bool)
    if np.any(finite):
        threshold = np.percentile(branch_medians[finite], 75.0)
        selected = finite & (branch_medians >= threshold)
    top_quartile_mean = nanmean(branch_signals[:, selected], axis=1)

    vessel = vessel_name.lower()
    return [
        _save_lumen_size_plot(
            output,
            branch_signals,
            [f"Branch {int(branch_id)}" for branch_id in branch_ids],
            title=f"{vessel_name} lumen size by branch",
            filename=f"lumen_size_by_branch_{vessel}.png",
            period_seconds=period_seconds,
        ),
        _save_lumen_size_plot(
            output,
            top_quartile_mean[:, None],
            [f"Mean of top-quartile branches (n={int(np.sum(selected))})"],
            title=f"{vessel_name} lumen size: top-quartile branches",
            filename=f"lumen_size_by_branch_top_quartile_{vessel}.png",
            period_seconds=period_seconds,
        ),
    ]


def _save_lumen_size_plot(
    output,
    signals: np.ndarray,
    labels: list[str],
    *,
    title: str,
    filename: str,
    period_seconds: float,
) -> Path:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(7.4, 4.2))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    # Per-beat Fourier interpolation omits the repeated cycle endpoint.
    time_seconds = np.linspace(0.0, period_seconds, signals.shape[0], endpoint=False)
    for signal, label in zip(signals.T, labels):
        if np.any(np.isfinite(signal)):
            ax.plot(time_seconds, signal, label=label)
    if ax.lines:
        if len(ax.lines) > 1:
            ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize="small")
        else:
            ax.legend(loc="best", fontsize="small")
    else:
        ax.text(
            0.5,
            0.5,
            "No finite lumen-size data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    ax.set(
        title=title,
        xlabel="Time (s)",
        ylabel="Lumen size (pixels)",
        xlim=(0.0, period_seconds),
    )
    ax.grid(True, alpha=0.25)
    path = output.path_for(OutputType.PNG, f"lumen_size/{filename}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    return path
