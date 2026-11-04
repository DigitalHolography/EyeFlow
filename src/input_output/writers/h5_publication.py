"""Serialize accepted output records and publish one complete HDF5 file."""

import os
import shutil
from contextlib import nullcontext
from pathlib import Path
from uuid import uuid4
from time import perf_counter

import h5py

from .h5_storage import _create_dataset


def _write_attributes(node, record):
    for key, value in record.attrs.items():
        dtype = record.attr_dtypes.get(key)
        if dtype is None:
            node.attrs[key] = value
        else:
            node.attrs.create(key, value, dtype=dtype)


def copy_output_node(model, source, target, name=None):
    """Copy a managed subtree into an isolated HDF5 export snapshot."""
    name = name or source.name.rsplit("/", 1)[-1]
    record = model.records[source.path]
    if record.kind == "group":
        node = target.create_group(name)
        for child_name, child in source.items():
            copy_output_node(model, child, node, child_name)
    else:
        node = _create_dataset(
            target, name, model.storage.read(record.value), h5_options=record.options
        )
    _write_attributes(node, record)


def publish_output_model(model):
    """Serialize locally, then atomically replace the destination filename."""
    started = perf_counter()
    final = Path(model.filename)
    local = model.storage.path(".h5")
    incoming = final.with_name(f".ef-{uuid4().hex[:16]}.h5")
    budget = model.storage.budget
    admission = (
        budget.compute(1, min(16 * 2**20, budget.memory_limit // 8))
        if budget is not None and hasattr(budget, "compute")
        else nullcontext()
    )
    try:
        with admission, h5py.File(local, "w") as output:
            for path, record in sorted(
                model.records.items(), key=lambda pair: (pair[0].count("/"), pair[0])
            ):
                if record.kind == "group":
                    node = output.require_group(path) if path else output
                else:
                    parent, _, name = path.rpartition("/")
                    group = output.require_group(parent) if parent else output
                    node = _create_dataset(
                        group,
                        name,
                        model.storage.read(record.value),
                        h5_options=record.options,
                    )
                _write_attributes(node, record)
        final.parent.mkdir(parents=True, exist_ok=True)
        model.storage.metrics.add("h5_serialization_seconds", perf_counter() - started)
        shutil.copyfile(local, incoming)
        os.replace(incoming, final)
    finally:
        model.storage.metrics.add("final_publication_seconds", perf_counter() - started)
        incoming.unlink(missing_ok=True)
        local.unlink(missing_ok=True)
