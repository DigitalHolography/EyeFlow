"""MATLAB-referenced fitting and minimal spatial normalization of profiles."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.math import interpft_real, next_power_of_two


@dataclass(frozen=True)
class ProfileData:
    """A transverse profile and the x coordinates it is sampled on."""

    velocity: np.ndarray
    x_micrometers: np.ndarray


@dataclass(frozen=True)
class ProfileProcessingResult:
    raw_profile: ProfileData
    interpolated_profile: ProfileData
    raw_x_micrometers: np.ndarray
    centered_velocity: np.ndarray
    centered_x_micrometers: np.ndarray
    center_micrometers: np.ndarray
    lumen_edges_micrometers: np.ndarray
    centering_fit_r_squared: np.ndarray
    poiseuille_coefficients: np.ndarray
    poiseuille_origin_micrometers: np.ndarray
    poiseuille_roots_micrometers: np.ndarray
    poiseuille_r_squared: np.ndarray


def process_velocity_profiles(
    profiles: np.ndarray,
    *,
    pixel_size_mm: float,
    velocity_profile_threshold: float,
    interpolation_points: int = 256,
) -> ProfileProcessingResult:
    """Center profiles while retaining both spatial representations.

    ``raw_profile`` contains the input samples on their original transverse
    pixel grid. ``interpolated_profile`` contains the same frames resampled
    onto a common, centered grid between the fitted lumen roots.

    Poiseuille values reproduce ``computeVesselCrossSection.m``: the quadratic
    is fit to samples above the configured fraction of the time-mean maximum,
    and its roots are evaluated at that same velocity level.

    The centered export follows the article rather than the hydrodynamic crop:
    a second quadratic is fit to all positive samples of the time-mean profile,
    its zero roots define the lumen, and every frame is linearly interpolated
    between root anchors. No temporal filtering, velocity clipping, or
    extrapolation beyond those anchors is applied.
    """

    values = np.asarray(profiles, dtype=np.float32)
    if values.ndim != 4:
        raise ValueError(
            "profiles must have shape (radius, branch, frame, transverse_sample)."
        )
    radius_count, branch_count, frame_count, sample_count = values.shape
    if interpolation_points < 2:
        raise ValueError("interpolation_points must be at least 2.")
    raw_x_um = (
        np.arange(sample_count, dtype=np.float32)
        - np.float32((sample_count - 1) / 2.0)
    ) * np.float32(pixel_size_mm * 1000.0)
    segment_shape = (radius_count, branch_count)
    centered = np.full(
        (*segment_shape, frame_count, interpolation_points),
        np.nan,
        dtype=np.float32,
    )
    centered_x = np.full(
        (*segment_shape, interpolation_points),
        np.nan,
        dtype=np.float32,
    )
    center_um = np.full(segment_shape, np.nan, dtype=np.float32)
    lumen_edges_um = np.full((*segment_shape, 2), np.nan, dtype=np.float32)
    centering_r2 = np.full(segment_shape, np.nan, dtype=np.float32)
    poiseuille_coefficients = np.full((*segment_shape, 3), np.nan, dtype=np.float32)
    poiseuille_origin_um = np.full(segment_shape, np.nan, dtype=np.float32)
    poiseuille_roots_um = np.full((*segment_shape, 2), np.nan, dtype=np.float32)
    poiseuille_r2 = np.full(segment_shape, np.nan, dtype=np.float32)

    for segment_index in np.ndindex(segment_shape):
        segment = values[segment_index]
        mean_profile = _nanmean(segment, axis=0)
        poiseuille = _matlab_poiseuille_fit(
            raw_x_um,
            mean_profile,
            velocity_profile_threshold,
        )
        if poiseuille is not None:
            coefficients, origin_um, roots, fit_r2 = poiseuille
            poiseuille_coefficients[segment_index] = coefficients
            poiseuille_origin_um[segment_index] = np.float32(origin_um)
            poiseuille_roots_um[segment_index] = roots
            poiseuille_r2[segment_index] = np.float32(fit_r2)

        centering = _positive_zero_root_fit(raw_x_um, mean_profile)
        if centering is None:
            continue
        roots, fit_r2 = centering
        center = float(np.mean(roots))
        radius = float((roots[1] - roots[0]) / 2.0)
        target_centered = np.linspace(
            -radius,
            radius,
            interpolation_points,
            dtype=np.float32,
        )
        target_absolute = target_centered + np.float32(center)
        for frame_index, profile in enumerate(segment):
            centered[segment_index + (frame_index,)] = _interpolate_with_root_anchors(
                raw_x_um,
                profile,
                target_absolute,
                roots,
            )
        centered_x[segment_index] = target_centered
        center_um[segment_index] = np.float32(center)
        lumen_edges_um[segment_index] = roots
        centering_r2[segment_index] = np.float32(fit_r2)

    return ProfileProcessingResult(
        raw_profile=ProfileData(
            velocity=values,
            x_micrometers=raw_x_um,
        ),
        interpolated_profile=ProfileData(
            velocity=centered,
            x_micrometers=centered_x,
        ),
        raw_x_micrometers=raw_x_um,
        centered_velocity=centered,
        centered_x_micrometers=centered_x,
        center_micrometers=center_um,
        lumen_edges_micrometers=lumen_edges_um,
        centering_fit_r_squared=centering_r2,
        poiseuille_coefficients=poiseuille_coefficients,
        poiseuille_origin_micrometers=poiseuille_origin_um,
        poiseuille_roots_micrometers=poiseuille_roots_um,
        poiseuille_r_squared=poiseuille_r2,
    )


def interpolate_velocity_profiles_per_beat(
    velocity_profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> np.ndarray:
    """Interpolate profiles to ``(x, time, beat, branch, radius)``.

    The temporal interpolation is the same ``interpft`` operation used by the
    standard EyeFlow per-beat signal output. It adds no filtering, clipping,
    rectification, or other profile processing.
    """
    profiles = np.asarray(velocity_profiles, dtype=np.float32)
    if profiles.ndim != 4:
        raise ValueError(
            "velocity_profiles must have shape "
            "(radius, branch, frame, transverse_sample)."
        )

    radius_count, branch_count, frame_count, x_count = profiles.shape
    boundaries = normalize_cycle_boundaries(
        cycle_boundary_indexes,
        frame_count,
        index_base=index_base,
    )
    beat_count = boundaries.size - 1
    time_count = next_power_of_two(int(np.max(np.diff(boundaries))))
    output = np.full(
        (x_count, time_count, beat_count, branch_count, radius_count),
        np.nan,
        dtype=np.float32,
    )

    for radius_index in range(radius_count):
        for branch_index in range(branch_count):
            for beat_index in range(beat_count):
                start = int(boundaries[beat_index])
                stop = int(boundaries[beat_index + 1]) + 1
                beat_profiles = profiles[
                    radius_index,
                    branch_index,
                    start:stop,
                    :,
                ]
                for x_index in range(x_count):
                    output[
                        x_index,
                        :,
                        beat_index,
                        branch_index,
                        radius_index,
                    ] = interpft_real(
                        beat_profiles[:, x_index],
                        time_count + 1,
                    )[:-1]

    return output


def _matlab_poiseuille_fit(
    x_um: np.ndarray,
    profile: np.ndarray,
    threshold_ratio: float,
) -> tuple[np.ndarray, float, np.ndarray, float] | None:
    finite = np.isfinite(profile)
    if not np.any(finite):
        return None
    maximum = float(np.nanmax(profile))
    threshold = float(threshold_ratio) * maximum
    central = np.flatnonzero(finite & (profile > threshold))
    if central.size < 3:
        return None
    center_x = float(np.mean(x_um[central]))
    fit_x = x_um[central].astype(np.float64) - center_x
    fit_y = profile[central].astype(np.float64)
    coefficients = np.polyfit(fit_x, fit_y, 2)
    roots = _quadratic_roots(
        coefficients[0],
        coefficients[1],
        coefficients[2] - threshold,
    )
    if roots is None:
        return None
    return (
        coefficients.astype(np.float32),
        center_x,
        np.asarray(roots, dtype=np.float32),
        _r_squared(fit_y, np.polyval(coefficients, fit_x)),
    )


def _positive_zero_root_fit(
    x_um: np.ndarray,
    profile: np.ndarray,
) -> tuple[np.ndarray, float] | None:
    positive = np.isfinite(profile) & (profile > 0)
    if np.count_nonzero(positive) < 3:
        return None
    fit_x = x_um[positive].astype(np.float64)
    fit_y = profile[positive].astype(np.float64)
    coefficients = np.polyfit(fit_x, fit_y, 2)
    if not np.all(np.isfinite(coefficients)) or coefficients[0] >= 0:
        return None
    roots = _quadratic_roots(*coefficients)
    if roots is None or not roots[0] < roots[1]:
        return None
    center = float(np.mean(roots))
    if not roots[0] < center < roots[1]:
        return None
    return (
        np.asarray(roots, dtype=np.float32),
        _r_squared(fit_y, np.polyval(coefficients, fit_x)),
    )


def _interpolate_with_root_anchors(
    x_um: np.ndarray,
    profile: np.ndarray,
    target_x_um: np.ndarray,
    roots_um: np.ndarray,
) -> np.ndarray:
    finite = np.isfinite(profile)
    inside = finite & (x_um > roots_um[0]) & (x_um < roots_um[1])
    source_x = np.concatenate(
        (
            np.asarray([roots_um[0]], dtype=np.float32),
            x_um[inside],
            np.asarray([roots_um[1]], dtype=np.float32),
        )
    )
    source_y = np.concatenate(
        (
            np.asarray([0.0], dtype=np.float32),
            profile[inside],
            np.asarray([0.0], dtype=np.float32),
        )
    )
    if source_x.size < 3:
        return np.full(target_x_um.shape, np.nan, dtype=np.float32)
    order = np.argsort(source_x)
    return np.interp(
        target_x_um,
        source_x[order],
        source_y[order],
        left=np.nan,
        right=np.nan,
    ).astype(np.float32)


def _quadratic_roots(a: float, b: float, c: float) -> tuple[float, float] | None:
    if not np.all(np.isfinite([a, b, c])) or a == 0:
        return None
    discriminant = b * b - 4.0 * a * c
    if discriminant <= 0:
        return None
    root_delta = np.sqrt(discriminant)
    roots = np.sort(
        np.asarray(
            [
                (-b + root_delta) / (2.0 * a),
                (-b - root_delta) / (2.0 * a),
            ],
            dtype=np.float64,
        )
    )
    return float(roots[0]), float(roots[1])


def _r_squared(observed: np.ndarray, fitted: np.ndarray) -> float:
    residual = float(np.sum((observed - fitted) ** 2))
    total = float(np.sum((observed - np.mean(observed)) ** 2))
    return np.nan if total == 0 else 1.0 - residual / total


def _nanmean(values: np.ndarray, *, axis):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)
