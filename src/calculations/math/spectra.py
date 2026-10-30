"""Array-only descriptors of calibrated power spectra and fractional bin intervals.

Signed spectra retain their stored bin order for interval areas. Probability
quantiles fold calibrated frequencies by magnitude, adding equal-coordinate
masses. Valley normalization divides by the true included spectral minimum;
Shannon entropy subsequently uses probabilities and is expressed in bits.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import xlogy

PERCENTILES = np.arange(101, dtype=np.float64)


def interval_area(values, lower: float, upper: float, spacing: float):
    """Integrate a piecewise-linear spectrum on fractional bin coordinates.

    Only samples supporting the interval participate. Integer endpoints do
    not consult adjacent bins, so NaNs outside an interval cannot poison it.
    """
    values = np.asarray(values, dtype=np.float64)
    if not 0 <= lower < upper <= values.shape[-1] - 1:
        raise ValueError("Integration bounds must be ordered coordinates within the bin domain.")
    interior = np.arange(np.floor(lower) + 1, np.ceil(upper), dtype=int)
    coordinates = np.concatenate(([lower], interior, [upper]))

    def sample(position):
        index = int(np.floor(position))
        fraction = position - index
        if fraction == 0:
            return values[..., index]
        return (1 - fraction) * values[..., index] + fraction * values[..., index + 1]

    samples = np.concatenate((
        sample(lower)[..., None], values[..., interior], sample(upper)[..., None],
    ), axis=-1)
    with np.errstate(invalid="ignore", over="ignore"):
        area = trapezoid(samples, x=coordinates, axis=-1) * spacing
    return np.where(np.isfinite(area), area, np.nan)


def frequency_distribution(frequency_magnitudes, probability):
    """Sort calibrated magnitudes and combine equal coordinates by adding mass."""
    x = np.abs(np.asarray(frequency_magnitudes, dtype=np.float64))
    p = np.asarray(probability, dtype=np.float64)
    if x.ndim != 1 or x.size == 0 or p.ndim < 1 or p.shape[-1] != x.size or not np.isfinite(x).all():
        raise ValueError("Distribution masses must match a finite frequency-coordinate vector.")
    order = np.argsort(x, kind="stable")
    unique, starts = np.unique(x[order], return_index=True)
    return unique, np.add.reduceat(p[..., order], starts, axis=-1)


def inverse_cdf(frequencies, masses, percentiles=PERCENTILES):
    """Discrete weighted inverse CDF; do not invent sub-bin frequency resolution."""
    q = np.asarray(percentiles, dtype=np.float64)
    if q.ndim != 1 or not np.isfinite(q).all() or np.any((q < 0) | (q > 100)):
        raise ValueError("Percentiles must be a finite vector within [0, 100].")
    x, weights = frequency_distribution(frequencies, masses)
    if weights.ndim != 1:
        raise ValueError("inverse_cdf requires one distribution.")
    if not np.isfinite(weights).all() or np.any(weights < 0) or weights.sum() <= 0:
        return np.full(q.shape, np.nan)
    occupied = weights > 0
    x, weights = x[occupied], weights[occupied]
    cumulative = np.cumsum(weights / weights.sum())
    cumulative[-1] = 1.
    return x[np.searchsorted(cumulative, q / 100, side="left")]


def describe_spectra(spectra, frequencies_khz):
    """Describe spectra with axes (..., selected_freq), preserving invalid rows."""
    values = np.asarray(spectra, dtype=np.float64)
    frequencies = np.asarray(frequencies_khz, dtype=np.float64)
    if values.ndim < 1 or values.shape[-1] < 2 or frequencies.shape != values.shape[-1:]:
        raise ValueError("Spectra must share a frequency axis containing at least two bins.")
    if not np.isfinite(frequencies).all():
        raise ValueError("Calibrated frequency coordinates must be finite.")
    rows = values.reshape(-1, values.shape[-1])
    probability = np.full(rows.shape, np.nan)
    normalized = np.full(rows.shape, np.nan)
    valley = np.full(rows.shape[0], np.nan)
    entropies = np.full(rows.shape[0], np.nan)
    quantiles = np.full((rows.shape[0], PERCENTILES.size), np.nan)
    for index, row in enumerate(rows):
        if not np.isfinite(row).all() or np.any(row < 0) or np.max(row) <= 0:
            continue
        valley[index] = np.min(row)
        if valley[index] > 0:
            with np.errstate(over="ignore", invalid="ignore"):
                normalized[index] = row / valley[index]
            if not np.isfinite(normalized[index]).all():
                normalized[index] = np.nan
        # Scaling by the maximum avoids overflow and is algebraically equal
        # to R/sum(R). Probabilities remain defined when the valley is zero.
        scaled = row / np.max(row)
        probability[index] = scaled / np.sum(scaled)
        entropies[index] = -np.sum(xlogy(probability[index], probability[index])) / np.log(2)
        quantiles[index] = inverse_cdf(np.abs(frequencies), probability[index])
    shape = values.shape[:-1]
    return {
        "ValleyMinimum": valley.reshape(shape),
        "ValleyNormalizedSpectrum": normalized.reshape(values.shape),
        "SpectralProbability": probability.reshape(values.shape),
        "SpectralEntropy": entropies.reshape(shape),
        "PowerQuantileFrequency": quantiles.reshape((*shape, PERCENTILES.size)),
    }


