"""Fast, isolated filesystem defaults shared by the test suite."""

from __future__ import annotations

import os

import pytest


@pytest.fixture(scope="session", autouse=True)
def _local_scratch_directory(tmp_path_factory):
    previous = os.environ.get("EYEFLOW_SCRATCH_DIR")
    os.environ["EYEFLOW_SCRATCH_DIR"] = str(tmp_path_factory.mktemp("scratch"))
    yield
    if previous is None:
        os.environ.pop("EYEFLOW_SCRATCH_DIR", None)
    else:
        os.environ["EYEFLOW_SCRATCH_DIR"] = previous
