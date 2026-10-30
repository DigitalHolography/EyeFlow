"""Trusted invocation transport with NumPy payloads stored once in local .npy."""

import pickle
from pathlib import Path

import numpy as np
from joblib.externals import cloudpickle

from .input_service import InputNode


def dump_result(value, storage):
    class Pickler(cloudpickle.CloudPickler):
        def persistent_id(self, obj):
            if isinstance(obj, InputNode):
                return "input", obj.file.slot, obj.path
            if (
                isinstance(obj, np.ndarray)
                and obj.nbytes >= 64 * 1024
                and obj.dtype.kind in "biufc"
            ):
                return "array", storage.store_array(obj)
            return None

    path = storage.path(".pickle")
    with path.open("wb") as stream:
        Pickler(stream, protocol=5).dump(value)
    return str(path)


def load_result(path, storage, inputs):
    class Unpickler(pickle.Unpickler):
        def persistent_load(self, identity):
            if identity[0] == "input":
                return inputs[identity[1]][identity[2]]
            if identity[0] == "array":
                return storage.read(identity[1])
            raise pickle.UnpicklingError("Unknown EyeFlow transfer reference")

    with Path(path).open("rb") as stream:
        return Unpickler(stream).load()
