"""Optional numeric acceleration backends used by large batched operations."""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache


@dataclass(frozen=True)
class CuPyBackend:
    cupy: object
    ndimage: object


@lru_cache(maxsize=1)
def optional_cupy_backend() -> CuPyBackend | None:
    """Return a usable CuPy backend unless CPU execution was requested.

    ``EYEFLOW_COMPUTE_BACKEND`` accepts ``auto`` (the default), ``cpu``, or
    ``cupy``. Auto mode falls back silently when CuPy or a CUDA device is not
    available; explicit CuPy mode raises a clear runtime error instead.
    """
    preference = os.environ.get("EYEFLOW_COMPUTE_BACKEND", "auto").strip().lower()
    if preference not in {"auto", "cpu", "cupy"}:
        raise ValueError(
            "EYEFLOW_COMPUTE_BACKEND must be 'auto', 'cpu', or 'cupy'."
        )
    if preference == "cpu":
        return None

    try:
        import cupy
        from cupyx.scipy import ndimage

        if int(cupy.cuda.runtime.getDeviceCount()) < 1:
            raise RuntimeError("no CUDA device is available")
        cupy.empty((1,), dtype=cupy.float32)
    except Exception as exc:
        if preference == "cupy":
            raise RuntimeError(
                "CuPy acceleration was requested but no usable CuPy/CUDA backend "
                "is available."
            ) from exc
        return None
    return CuPyBackend(cupy=cupy, ndimage=ndimage)
