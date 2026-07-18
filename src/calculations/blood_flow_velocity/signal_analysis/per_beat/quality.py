"""Shared beat-quality assessment for paired arterial and venous waveforms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signal import PerBeatSignalAnalysisResult

DURATION_TOO_SHORT = np.uint16(1)
DURATION_TOO_LONG = np.uint16(2)
ARTERY_RESIDUAL = np.uint16(4)
VEIN_RESIDUAL = np.uint16(8)
ARTERY_TEMPLATE = np.uint16(16)
VEIN_TEMPLATE = np.uint16(32)
INSUFFICIENT_WAVEFORM = np.uint16(64)


@dataclass(frozen=True)
class BeatQualitySettings:
    """Conservative, dimensionless thresholds used by beat-level motion QC."""

    raw_bandlimited_residual_limit: float = 0.20
    paired_template_distance_limit: float = 0.30
    minimum_template_beats: int = 3
    minimum_period_ratio: float = 0.55
    maximum_period_ratio: float = 1.60

    def __post_init__(self) -> None:
        if self.raw_bandlimited_residual_limit <= 0:
            raise ValueError("raw_bandlimited_residual_limit must be positive.")
        if self.paired_template_distance_limit <= 0:
            raise ValueError("paired_template_distance_limit must be positive.")
        if self.minimum_template_beats < 1:
            raise ValueError("minimum_template_beats must be positive.")
        if not 0 < self.minimum_period_ratio < self.maximum_period_ratio:
            raise ValueError(
                "Period-ratio limits must be positive and strictly increasing."
            )


@dataclass(frozen=True)
class VesselBeatQualityScores:
    raw_bandlimited_residual: np.ndarray
    template_correlation: np.ndarray


@dataclass(frozen=True)
class BeatQualityResult:
    accepted_mask: np.ndarray
    rejection_flags: np.ndarray
    quality_score: np.ndarray
    period_ratio: np.ndarray
    artery: VesselBeatQualityScores
    vein: VesselBeatQualityScores
    settings: BeatQualitySettings


def assess_beat_quality(
    artery: PerBeatSignalAnalysisResult,
    vein: PerBeatSignalAnalysisResult,
    *,
    period_ratio: np.ndarray | None = None,
    settings: BeatQualitySettings | None = None,
) -> BeatQualityResult:
    """Assess paired beats once and return one mask shared by both vessels."""

    settings = settings or BeatQualitySettings()
    artery_raw, artery_band = _paired_arrays(artery)
    vein_raw, vein_band = _paired_arrays(vein)
    if artery_raw.shape != vein_raw.shape:
        raise ValueError("Artery and vein candidate beat arrays must have the same shape.")

    n_beats = artery_raw.shape[0]
    ratios = _period_ratios(period_ratio, n_beats)
    flags = np.zeros(n_beats, dtype=np.uint16)
    _flag_duration_outliers(flags, ratios, settings)

    artery_residual, artery_usable = _normalized_residual(artery_raw, artery_band)
    vein_residual, vein_usable = _normalized_residual(vein_raw, vein_band)
    usable = artery_usable & vein_usable
    flags[~usable] |= INSUFFICIENT_WAVEFORM

    artery_bad_residual = artery_residual > settings.raw_bandlimited_residual_limit
    vein_bad_residual = vein_residual > settings.raw_bandlimited_residual_limit
    flags[artery_bad_residual] |= ARTERY_RESIDUAL
    flags[vein_bad_residual] |= VEIN_RESIDUAL

    duration_valid = (flags & (DURATION_TOO_SHORT | DURATION_TOO_LONG)) == 0
    artery_reference = duration_valid & usable & ~artery_bad_residual
    vein_reference = duration_valid & usable & ~vein_bad_residual
    artery_correlation = _template_correlation(
        artery_band,
        artery_reference,
        settings.minimum_template_beats,
    )
    vein_correlation = _template_correlation(
        vein_band,
        vein_reference,
        settings.minimum_template_beats,
    )
    paired_distance = np.sqrt(
        np.maximum(0.0, 1.0 - artery_correlation)
        * np.maximum(0.0, 1.0 - vein_correlation)
    ).astype(np.float32)
    paired_bad = (
        np.isfinite(paired_distance)
        & (paired_distance > settings.paired_template_distance_limit)
    )
    flags[paired_bad] |= ARTERY_TEMPLATE | VEIN_TEMPLATE

    template_score = np.nan_to_num(
        paired_distance / np.float32(settings.paired_template_distance_limit),
        nan=0.0,
        posinf=np.inf,
        neginf=0.0,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        duration_score = np.maximum(
            np.float32(settings.minimum_period_ratio) / ratios,
            ratios / np.float32(settings.maximum_period_ratio),
        )
    duration_score = np.nan_to_num(
        duration_score,
        nan=np.inf,
        posinf=np.inf,
        neginf=np.inf,
    )
    quality_score = np.maximum.reduce(
        (
            duration_score,
            artery_residual / np.float32(settings.raw_bandlimited_residual_limit),
            vein_residual / np.float32(settings.raw_bandlimited_residual_limit),
            template_score,
        )
    ).astype(np.float32)
    quality_score[~usable] = np.inf
    accepted = flags == 0
    return BeatQualityResult(
        accepted_mask=accepted,
        rejection_flags=flags,
        quality_score=quality_score,
        period_ratio=ratios,
        artery=VesselBeatQualityScores(artery_residual, artery_correlation),
        vein=VesselBeatQualityScores(vein_residual, vein_correlation),
        settings=settings,
    )


def _paired_arrays(
    result: PerBeatSignalAnalysisResult,
) -> tuple[np.ndarray, np.ndarray]:
    raw = np.asarray(result.velocity_signal_per_beat, dtype=np.float32)
    band = np.asarray(result.velocity_signal_per_beat_band_limited, dtype=np.float32)
    if raw.ndim != 2 or raw.shape != band.shape:
        raise ValueError("Raw and band-limited beat arrays must have matching 2D shapes.")
    return raw, band


def _period_ratios(period_ratio: np.ndarray | None, n_beats: int) -> np.ndarray:
    if period_ratio is None:
        return np.ones(n_beats, dtype=np.float32)
    ratios = np.asarray(period_ratio, dtype=np.float32).reshape(-1)
    if ratios.size != n_beats:
        raise ValueError("period_ratio must contain one value per beat interval.")
    return ratios


def _flag_duration_outliers(
    flags: np.ndarray,
    ratios: np.ndarray,
    settings: BeatQualitySettings,
) -> None:
    flags[~np.isfinite(ratios)] |= INSUFFICIENT_WAVEFORM
    flags[ratios < settings.minimum_period_ratio] |= DURATION_TOO_SHORT
    flags[ratios > settings.maximum_period_ratio] |= DURATION_TOO_LONG


def _normalized_residual(
    raw: np.ndarray,
    band: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    finite = np.all(np.isfinite(raw), axis=1) & np.all(np.isfinite(band), axis=1)
    scale = (
        np.percentile(band, 95, axis=1) - np.percentile(band, 5, axis=1)
    ).astype(np.float32)
    scale_floor = np.maximum(
        np.max(np.abs(band), axis=1).astype(np.float32) * np.float32(1e-6),
        np.float32(np.finfo(np.float32).eps),
    )
    usable = finite & (scale > scale_floor)
    residual = np.full(raw.shape[0], np.inf, dtype=np.float32)
    residual[usable] = (
        np.sqrt(np.mean((raw[usable] - band[usable]) ** 2, axis=1))
        / scale[usable]
    ).astype(np.float32)
    return residual, usable


def _template_correlation(
    band: np.ndarray,
    reference_mask: np.ndarray,
    minimum_template_beats: int,
) -> np.ndarray:
    correlations = np.full(band.shape[0], np.nan, dtype=np.float32)
    if int(np.sum(reference_mask)) < int(minimum_template_beats):
        return correlations

    standardized = _standardize_rows(band)
    template = np.median(standardized[reference_mask], axis=0)
    template = template - np.mean(template)
    template_norm = float(np.linalg.norm(template))
    if template_norm <= np.finfo(np.float32).eps:
        return correlations

    centered = standardized - np.mean(standardized, axis=1, keepdims=True)
    row_norm = np.linalg.norm(centered, axis=1)
    valid = row_norm > np.finfo(np.float32).eps
    correlations[valid] = (
        (centered[valid] @ template) / (row_norm[valid] * template_norm)
    ).astype(np.float32)
    return correlations


def _standardize_rows(values: np.ndarray) -> np.ndarray:
    center = np.median(values, axis=1, keepdims=True)
    scale = (
        np.percentile(values, 95, axis=1) - np.percentile(values, 5, axis=1)
    ).reshape(-1, 1)
    scale = np.maximum(scale, np.finfo(np.float32).eps)
    return ((values - center) / scale).astype(np.float32)
