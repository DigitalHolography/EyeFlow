"""Acquisition I/O ownership, with independent input and retained output lifetimes."""

from contextlib import ExitStack
import multiprocessing
from threading import RLock

from input_output.writers.h5 import open_h5
from utils.logger import Logger

from .input.service import InputService
from .output.service import OutputService


class IOManager:
    """Coordinate independent input/output lanes over borrowed product storage.

    AcquisitionCompletion owns exports, final metadata, publication and storage
    release. Closing these lanes never closes the caller's model storage.
    """

    def __init__(self, model, *, budget=None, context=None, output_policy=None):
        self.model = model
        self.context = context or multiprocessing.get_context("spawn")
        self.budget = budget if budget is not None else model.storage.budget
        self._inputs = ExitStack()
        self._sources = None
        self.input = None
        self.files = {}
        self._input_closed = self._closed = False
        self._lifecycle = RLock()
        self.output = OutputService(model, budget=self.budget, context=self.context, policy=output_policy)

    def open_sources(self, hd_path, dv_path):
        if self._input_closed or self._sources is not None:
            raise RuntimeError("Acquisition inputs have already been opened or closed")
        with Logger.perf_stage("run/input_hdf5_open"):
            self._sources = {
                slot: self._inputs.enter_context(open_h5(path, "r")) if path is not None else None
                for slot, path in (("hd", hd_path), ("dv", dv_path))
            }
        return self._sources["hd"], self._sources["dv"]

    def start_input(self, policy=None):
        if self._input_closed or self._sources is None or self.input is not None:
            raise RuntimeError("Input service cannot be started in this acquisition state")
        self.input = InputService(self._sources, self.model.storage, self.context, policy=policy)
        self._inputs.callback(self.input.close)
        client = self.input.client()
        self._inputs.callback(client.close)
        self.files = self.input.files(client)
        return self.input

    def close_input(self):
        if not self._input_closed:
            self._input_closed = True
            self._inputs.close()

    def check_health(self):
        self.output.check_health()
        self.output.submit(self.model)
        self.output.check_health()

    @property
    def closed(self):
        return self._closed

    def close(self):
        self.output.request_abort()
        with self._lifecycle:
            if self._closed:
                return
            self._closed = True
            with ExitStack() as cleanup:
                cleanup.callback(self.output.close)
                self.close_input()

    abort = close
