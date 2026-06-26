"""Figure-facing spectrum analysis orchestration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .pair import PairedSpectrumAnalysisResult, paired_spectrum_analysis
from .signal import SpectrumData, spectrum_signal_analysis


@dataclass(frozen=True)
class VesselSpectrumAnalysis:
    artery: SpectrumData
    vein: SpectrumData


def run_vessel_spectrum_analysis(
    artery_values: np.ndarray,
    vein_values: np.ndarray,
    dt_seconds: float,
) -> VesselSpectrumAnalysis:
    return VesselSpectrumAnalysis(
        artery=spectrum_signal_analysis(artery_values, dt_seconds),
        vein=spectrum_signal_analysis(vein_values, dt_seconds),
    )


def run_paired_vessel_spectrum_analysis(
    artery_values: np.ndarray,
    vein_values: np.ndarray,
    dt_seconds: float,
    beat_indexes: np.ndarray | None = None,
) -> PairedSpectrumAnalysisResult:
    return paired_spectrum_analysis(artery_values, vein_values, dt_seconds, beat_indexes)
