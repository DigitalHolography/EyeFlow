"""Model fits for generic vessel-aligned profiles."""

from .quadratic import (
    COUNT_OUTPUTS,
    DEFAULT_TIME_BLOCK_SIZE,
    DEFAULT_WEIGHT_POWER,
    FLOAT_OUTPUTS,
    border_weights,
    fit_quadratic_profiles,
)

__all__ = [
    "COUNT_OUTPUTS",
    "DEFAULT_TIME_BLOCK_SIZE",
    "DEFAULT_WEIGHT_POWER",
    "FLOAT_OUTPUTS",
    "border_weights",
    "fit_quadratic_profiles",
]
