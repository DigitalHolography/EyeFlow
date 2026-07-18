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
    interval_mask: np.ndarray | None = None,
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
        interval_mask=interval_mask,
    )


def _analyze_cycles(
    signal_array: np.ndarray,
    cycle_boundaries: np.ndarray,
    harmonic_count: int,
    *,
    interval_mask: np.ndarray | None = None,
) -> PerBeatSignalAnalysisResult:
    interval_pairs = np.column_stack(
        (cycle_boundaries[:-1], cycle_boundaries[1:])
    ).astype(np.int32, copy=False)
    if interval_mask is not None:
        mask = np.asarray(interval_mask, dtype=np.bool_).reshape(-1)
        if mask.size != interval_pairs.shape[0]:
            raise ValueError(
                "interval_mask must contain one value per boundary interval."
            )
        interval_pairs = interval_pairs[mask]
    if interval_pairs.shape[0] == 0:
        raise ValueError("At least one beat interval must be selected.")

    number_of_beats = int(interval_pairs.shape[0])
    n_fft = next_power_of_two(int(np.max(interval_pairs[:, 1] - interval_pairs[:, 0])))
    per_beat, per_beat_fft, band_limited = _empty_outputs(number_of_beats, n_fft)

    for beat_index, (start_index, stop_index) in enumerate(interval_pairs):
        start = int(start_index)
        stop = int(stop_index) + 1
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
