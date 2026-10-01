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
REGISTRATION_METHOD_OUTPUT_NAMES: dict[str, str] = {
    "classic_demons": "Demons",
    "symmetric_forces_demons": "SymmetricForcesDemons",
    "fast_symmetric_demons": "FastSymmetricForcesDemons",
    "diffeomorphic_demons": "DiffeomorphicDemons",
    "level_set_motion": "LevelSetMotion",
    "displacement_field": "DisplacementField",
}
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


def registration_method_output_name(value: object) -> str:
    """Return the canonical HDF5 path segment for a registration method."""

    method = str(value).strip()
    try:
        return REGISTRATION_METHOD_OUTPUT_NAMES[method]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported displacement registration method: {method!r}."
        ) from exc
