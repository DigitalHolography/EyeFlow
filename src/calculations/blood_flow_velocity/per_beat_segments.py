"""Segment-level per-beat velocity calculations from Tools/exportProfilesToH5.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .per_beat_signal import PerBeatSignalAnalysisResult, per_beat_signal_analysis


@dataclass(frozen=True)
class PerBeatSegmentAnalysisResult:
    velocity_signal_per_beat_per_segment: np.ndarray
    velocity_signal_per_beat_per_segment_fft: np.ndarray
    velocity_signal_per_beat_per_segment_band_limited: np.ndarray


def per_beat_segment_analysis(
    segment_velocity_signals,
    sys_idx_list,
    band_limited_signal_harmonic_count: int,
    *,
    index_base: int | None = None,
) -> PerBeatSegmentAnalysisResult:
    segments = _segment_array(segment_velocity_signals)
    first = per_beat_signal_analysis(
        segments[0, 0, :],
        sys_idx_list,
        band_limited_signal_harmonic_count,
        index_base=index_base,
    )
    outputs = _empty_segment_outputs(segments.shape, first)
    _fill_segment_outputs(outputs, first, branch_index=0, radius_index=0)
    _fill_remaining_segments(
        outputs,
        segments,
        sys_idx_list,
        band_limited_signal_harmonic_count,
        index_base,
    )
    return PerBeatSegmentAnalysisResult(*outputs)


def _segment_array(segment_velocity_signals) -> np.ndarray:
    segments = np.asarray(segment_velocity_signals, dtype=np.float32)
    if segments.ndim != 3:
        raise ValueError(
            "segment_velocity_signals must have shape (radius, branch, frame)."
        )
    if segments.shape[0] == 0 or segments.shape[1] == 0:
        raise ValueError("segment_velocity_signals must include at least one segment.")
    return segments


def _empty_segment_outputs(
    segment_shape: tuple[int, int, int],
    sample: PerBeatSignalAnalysisResult,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_radii, n_branches, _ = segment_shape
    n_beats, n_samples = sample.velocity_signal_per_beat.shape
    shape = (n_samples, n_beats, n_branches, n_radii)
    raw = np.full(shape, np.nan, dtype=np.float32)
    fft = np.full(shape, np.nan + 0j, dtype=np.complex64)
    band_limited = np.full(shape, np.nan, dtype=np.float32)
    return raw, fft, band_limited


def _fill_remaining_segments(
    outputs: tuple[np.ndarray, np.ndarray, np.ndarray],
    segments: np.ndarray,
    sys_idx_list,
    harmonic_count: int,
    index_base: int | None,
) -> None:
    for radius_index in range(segments.shape[0]):
        for branch_index in range(segments.shape[1]):
            if radius_index == 0 and branch_index == 0:
                continue
            result = per_beat_signal_analysis(
                segments[radius_index, branch_index, :],
                sys_idx_list,
                harmonic_count,
                index_base=index_base,
            )
            _fill_segment_outputs(outputs, result, branch_index, radius_index)


def _fill_segment_outputs(
    outputs: tuple[np.ndarray, np.ndarray, np.ndarray],
    result: PerBeatSignalAnalysisResult,
    branch_index: int,
    radius_index: int,
) -> None:
    raw, fft, band_limited = outputs
    raw[:, :, branch_index, radius_index] = result.velocity_signal_per_beat.T
    fft[:, :, branch_index, radius_index] = result.velocity_signal_per_beat_fft.T
    band_limited[:, :, branch_index, radius_index] = (
        result.velocity_signal_per_beat_band_limited.T
    )
