"""Types and fixed values used by displacement-map estimation."""

from typing import Literal

RegistrationMethod = Literal[
    "classic_demons",
    "symmetric_forces_demons",
    "fast_symmetric_demons",
    "diffeomorphic_demons",
    "level_set_motion",
    "displacement_field",
    "bspline",
    "affine",
    "curvature_registration",
    "gpu_demons",
]
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
ITK_REGISTRATION_METHODS = frozenset(
    {
        "curvature_registration",
        "gpu_demons",
    }
)
CURVATURE_TIME_STEP = 1.0
CURVATURE_CONSTRAINT_WEIGHT = 0.1
