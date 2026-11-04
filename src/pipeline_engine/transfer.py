"""Trusted invocation transport with NumPy payloads stored once in local .npy."""

import pickle
from pathlib import Path

import numpy as np
import cloudpickle

from .input_service import InputNode


def dump_result(value, storage):
    # A writable array gets one copy per snapshot, even when it appears twice.
    # ProductStorage additionally reuses immutable arrays across snapshots.
    arrays = {}

    class Pickler(cloudpickle.CloudPickler):
        def persistent_id(self, obj):
            if isinstance(obj, InputNode):
                return "input", obj.file.slot, obj.path
            if (
                isinstance(obj, np.ndarray)
                and obj.nbytes >= 64 * 1024
                and obj.dtype.kind in "biufc"
            ):
                if id(obj) not in arrays:
                    arrays[id(obj)] = obj, storage.store_array(obj)
                return "array", arrays[id(obj)][1]
            return None

    path = storage.path(".pickle")
    with path.open("wb") as stream:
        Pickler(stream, protocol=5).dump(value)
    storage.metrics.add("serializations", 1)
    storage.metrics.add("manifest_bytes", path.stat().st_size)
    return str(path)


def load_result(path, storage, inputs, *, owner_root=None):
    references = {}

    class Unpickler(pickle.Unpickler):
        def persistent_load(self, identity):
            if identity[0] == "input":
                key = tuple(identity)
                if key not in references:
                    references[key] = inputs[identity[1]][identity[2]]
            elif identity[0] == "array":
                from data_manager import DataManager
                DataManager(storage, owner_id=owner_root).validate_storage(identity[1])
                key = "array", identity[1].path
                if key not in references:
                    references[key] = storage.read(identity[1])
            else:
                raise pickle.UnpicklingError("Unknown EyeFlow transfer reference")
            return references[key]

    with Path(path).open("rb") as stream:
        return Unpickler(stream).load()
