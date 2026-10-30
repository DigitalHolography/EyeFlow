"""Types and fixed values used by displacement-map estimation."""

from typing import Literal

RegistrationMethod = Literal[
    "classic_demons",
    "symmetric_forces_demons",
    "fast_symmetric_demons",
    "diffeomorphic_demons",
    "level_set_motion",
    "displacement_field",
]
DEFAULT_REGISTRATION_METHOD: RegistrationMethod = "level_set_motion"
PhotometricMode = Literal["none", "homomorphic", "local_contrast", "hybrid"]
PDE_REGISTRATION_METHODS = frozenset(
    {
        "classic_demons",
        "symmetric_forces_demons",
        "fast_symmetric_demons",
        "diffeomorphic_demons",
        "level_set_motion",
    }
)
