"""Per-beat blood-flow velocity signal calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.math import band_limited_ifft_abs, interpft_real, next_power_of_two

from ._signal_utils import normalize_cycle_boundaries


@dataclass(frozen=True)
class PerBeatSignalAnalysisResult:
    velocity_signal_per_beat: np.ndarray
    velocity_signal_per_beat_fft: np.ndarray
    velocity_signal_per_beat_band_limited: np.ndarray


def per_beat_signal_analysis(
    signal,
    sys_idx_list,
    band_limited_signal_harmonic_count: int,
    *,
    index_base: int | None = None,
) -> PerBeatSignalAnalysisResult:
    signal_array = np.asarray(signal, dtype=np.float32).reshape(-1)
    if signal_array.size == 0:
        raise ValueError("signal must contain at least one sample.")
    if band_limited_signal_harmonic_count < 1:
        raise ValueError("band_limited_signal_harmonic_count must be positive.")

    cycle_boundaries = normalize_cycle_boundaries(
        sys_idx_list,
        signal_array.size,
        index_base=index_base,
    )
    return _analyze_cycles(
        signal_array,
        cycle_boundaries,
        int(band_limited_signal_harmonic_count),
    )


def _analyze_cycles(
    signal_array: np.ndarray,
    cycle_boundaries: np.ndarray,
    harmonic_count: int,
) -> PerBeatSignalAnalysisResult:
    number_of_beats = int(cycle_boundaries.size - 1)
    n_fft = next_power_of_two(int(np.max(np.diff(cycle_boundaries))))
    per_beat, per_beat_fft, band_limited = _empty_outputs(number_of_beats, n_fft)

    for beat_index in range(number_of_beats):
        start = int(cycle_boundaries[beat_index])
        stop = int(cycle_boundaries[beat_index + 1]) + 1
        beat_interp = interpft_real(signal_array[start:stop], n_fft + 1)[:-1]
        beat_fft = np.fft.fft(beat_interp, n=n_fft).astype(np.complex64)
        per_beat[beat_index, :] = beat_interp
        per_beat_fft[beat_index, :] = beat_fft
        band_limited[beat_index, :] = _band_limited_signal(
            beat_fft,
            n_fft,
            harmonic_count,
        )

    return PerBeatSignalAnalysisResult(per_beat, per_beat_fft, band_limited)


def _empty_outputs(
    number_of_beats: int,
    n_fft: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    per_beat = np.full((number_of_beats, n_fft), np.nan, dtype=np.float32)
    per_beat_fft = np.full((number_of_beats, n_fft), np.nan + 0j, dtype=np.complex64)
    band_limited = np.full((number_of_beats, n_fft), np.nan, dtype=np.float32)
    return per_beat, per_beat_fft, band_limited


def _band_limited_signal(
    beat_fft: np.ndarray,
    n_fft: int,
    harmonic_count: int,
) -> np.ndarray:
    return band_limited_ifft_abs(beat_fft, n_fft, harmonic_count)
