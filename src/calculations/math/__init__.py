"""Reusable mathematical helpers for EyeFlow calculations."""

from .arrays import (
    as_float32_vector,
    as_nonnegative_int_indexes,
    finite_image,
    nan_to_mean,
    rescale,
    standardize,
)
from .filtering import butter_lowpass_filtfilt, normalized_lowpass_cutoff
from .fourier import (
    band_limited_ifft_abs,
    harmonic_pack,
    interpft_real,
    irfft_normalized,
    next_power_of_two,
    rfft_normalized,
    truncate_harmonics,
)
from .image import rotate_array_threshold, rotate_image_with_nan
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
    "as_float32_vector",
    "as_nonnegative_int_indexes",
    "band_limited_ifft_abs",
    "butter_lowpass_filtfilt",
    "finite_image",
    "harmonic_pack",
    "interpft_real",
    "irfft_normalized",
    "nanargmax",
    "nanargmin",
    "nanmax",
    "nanmean",
    "nan_to_mean",
    "nanmean_float32",
    "nanmedian",
    "nanmin",
    "nansum",
    "next_power_of_two",
    "normalized_lowpass_cutoff",
    "rescale",
    "rotate_array_threshold",
    "rotate_image_with_nan",
    "standardize",
    "rfft_normalized",
    "truncate_harmonics",
]
