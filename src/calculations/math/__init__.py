"""Reusable mathematical helpers for EyeFlow calculations."""

from .arrays import (
    as_float32_vector,
    as_nonnegative_int_indexes,
    finite_image,
    nan_to_mean,
    rescale,
    standardize,
)
from .bandlimiting import DEFAULT_BAND_LIMITED_HARMONICS, band_limit, reconstruct_band_limited
from .filtering import butter_lowpass_filtfilt, normalized_lowpass_cutoff
from .fourier import (
    harmonic_pack,
    interpft_axis0,
    interpft_real,
    irfft_normalized,
    next_power_of_two,
    rfft_normalized,
    truncate_harmonics,
)
from .periodic import SlidingWindowMethod, centered_sliding_window
from .statistics import (
    nanargmax,
    nanargmin,
    nanmax,
    nanmean,
    nanmean_float32,
    nanmedian,
    nanmin,
    nansum,
)

__all__ = [
    "SlidingWindowMethod",
    "as_float32_vector",
    "as_nonnegative_int_indexes",
    "band_limit",
    "reconstruct_band_limited",
    "DEFAULT_BAND_LIMITED_HARMONICS",
    "butter_lowpass_filtfilt",
    "centered_sliding_window",
    "finite_image",
    "harmonic_pack",
    "interpft_axis0",
    "interpft_real",
    "irfft_normalized",
    "nan_to_mean",
    "nanargmax",
    "nanargmin",
    "nanmax",
    "nanmean",
    "nanmean_float32",
    "nanmedian",
    "nanmin",
    "nansum",
    "next_power_of_two",
    "normalized_lowpass_cutoff",
    "rescale",
    "rfft_normalized",
    "standardize",
    "truncate_harmonics",
]
