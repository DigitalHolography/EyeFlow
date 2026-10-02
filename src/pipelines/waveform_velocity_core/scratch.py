"""RAM-backed HDF5 workspace for retinal velocity intermediates."""

from __future__ import annotations

from input_output.scratch_h5 import scratch_h5

def velocity_scratch_h5(_ctx):
    """Yield a non-persistent HDF5 file allocated entirely in RAM."""
    return scratch_h5(
        purpose="EyeFlow retinal velocity intermediates",
        filename_prefix="eyeflow-velocity",
    )


__all__ = ["velocity_scratch_h5"]
