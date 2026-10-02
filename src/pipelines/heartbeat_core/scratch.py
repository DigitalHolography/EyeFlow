"""RAM-backed HDF5 workspace for heartbeat intermediates."""

from __future__ import annotations

from input_output.scratch_h5 import scratch_h5


def heartbeat_scratch_h5(_ctx):
    """Yield a non-persistent HDF5 file allocated entirely in RAM."""
    return scratch_h5(
        purpose="EyeFlow heartbeat intermediates",
        filename_prefix="eyeflow-heartbeat",
    )
