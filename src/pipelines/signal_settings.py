"""Shared velocity and gap Fourier-bandlimiting settings."""

from calculations.math.bandlimiting import DEFAULT_BAND_LIMITED_HARMONICS, validate_harmonics


def band_limited_harmonics(ctx) -> int:
    """Read the shared non-DC count; translate legacy DC-inclusive counts."""
    settings = getattr(ctx, "attrs", {})
    value = settings.get("band_limited_harmonics")
    if value is not None:
        return validate_harmonics(value)
    for key in ("BandLimitedSignalHarmonicCount", "band_limited_signal_harmonic_count"):
        value = settings.get(key)
        if value is not None:
            count = validate_harmonics(value)
            if count < 1:
                raise ValueError(f"{key} must include at least DC.")
            return count - 1
    return DEFAULT_BAND_LIMITED_HARMONICS
