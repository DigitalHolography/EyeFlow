"""Keep native topology and array-only profile fits independent of execution."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    ("module", "forbidden"),
    [
        (
            "calculations.topology",
            ("calculations.vessel_segments", "calculations.compute_backend", "pipelines"),
        ),
        (
            "calculations.vessel_segments.profiles.fits.quadratic",
            (
                "calculations.topology",
                "calculations.vessel_segments.sampling",
                "calculations.vessel_segments.measurement",
                "calculations.compute_backend",
                "pipelines",
            ),
        ),
    ],
)
def test_cold_import_does_not_load_downstream_layers(module, forbidden):
    # A fresh process prevents unrelated test imports from concealing cycles.
    source_dir = str(Path(__file__).resolve().parents[1] / "src")
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (source_dir, environment.get("PYTHONPATH")))
    )
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import importlib, sys; "
            f"importlib.import_module({module!r}); "
            f"forbidden = {forbidden!r}; "
            "loaded = [name for name in sys.modules "
            "if any(name == prefix or name.startswith(prefix + '.') "
            "for prefix in forbidden)]; "
            "assert not loaded, loaded",
        ],
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
