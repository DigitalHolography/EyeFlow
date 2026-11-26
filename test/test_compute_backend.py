"""Backend selection shared by coordinator calculations and spawned workers."""

import os

from unittest.mock import patch

import pytest

from calculations.compute_backend import CuPyBackend, optional_cupy_backend
from pipeline_engine.execution.run_service import _batch_compute_backend


def test_cpu_selection_bypasses_an_available_cached_gpu(monkeypatch):
    backend = CuPyBackend(cupy=object(), ndimage=object())
    with patch(
        "calculations.compute_backend._probe_cupy_backend",
        return_value=(backend, None),
    ) as probe:
        monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
        assert optional_cupy_backend() is backend
        monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cpu")
        assert optional_cupy_backend() is None
        monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
        assert optional_cupy_backend() is backend
    assert probe.call_count == 2


def test_explicit_gpu_request_fails_when_startup_probe_fails(monkeypatch):
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cupy")
    with patch(
        "calculations.compute_backend._probe_cupy_backend",
        return_value=(None, RuntimeError("CUDA image kernel unavailable")),
    ):
        with pytest.raises(RuntimeError, match="no usable CuPy/CUDA backend") as error:
            optional_cupy_backend()
    assert "CUDA image kernel unavailable" in str(error.value.__cause__)


def test_batch_preflight_pins_auto_fallback_for_workers(monkeypatch):
    monkeypatch.delenv("EYEFLOW_COMPUTE_BACKEND", raising=False)

    with patch(
        "pipeline_engine.execution.run_service.optional_cupy_backend",
        return_value=None,
    ):
        with _batch_compute_backend():
            assert os.environ["EYEFLOW_COMPUTE_BACKEND"] == "cpu"

    assert "EYEFLOW_COMPUTE_BACKEND" not in os.environ
