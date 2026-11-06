"""Explicit ownership of pending exports and the final output payload store."""

import json
from pathlib import Path

from input_output.writers.h5 import write_root_attributes


class PendingOutput:
    """Keep result payloads alive until publication, then release them once."""

    def __init__(self, model, artifacts, *, owns_artifacts=False):
        self.model = model
        self.artifacts = artifacts
        self.owns_artifacts = owns_artifacts
        self._closed = False
        self._finished = False

    def finish(self):
        """Join exports, publish their statuses and results, and close storage."""
        if self._finished:
            return Path(self.model.filename)
        if self._closed:
            raise RuntimeError("Cannot publish a closed pending output.")
        try:
            self.artifacts.collect(wait=True)
            write_root_attributes(self.model, {
                "artifact_status": json.dumps(self.artifacts.statuses, sort_keys=True),
                "artifact_failures": json.dumps(self.artifacts.failures, sort_keys=True),
                "artifact_timings": json.dumps(self.artifacts.timings, sort_keys=True),
            })
            self.model.publish()
            self._finished = True
            return Path(self.model.filename)
        finally:
            self.close()

    def close(self):
        """Discard an unfinished result or release a completed one."""
        if self._closed:
            return
        self._closed = True
        try:
            if self.owns_artifacts:
                self.artifacts.__exit__(None, None, None)
        finally:
            self.model.storage.close()
