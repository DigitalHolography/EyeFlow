"""The single authoritative coordinator acceptance transaction boundary."""

import hashlib
import shutil
import pickle
from dataclasses import fields, is_dataclass, replace
from pathlib import Path
from threading import get_ident
from time import perf_counter
from uuid import uuid4

from data_manager import DataManager
from utils.logger import Logger
from input_output.writers.h5 import write_root_attributes

from ..stage_graph.journal import StageJournal
from ..transfer import load_result
from .identity import ProductId


def published_paths(value, source_root, destination_root):
    if isinstance(value, Path):
        path = value.resolve()
        return destination_root / path.relative_to(source_root) if path.is_relative_to(source_root) else value
    if isinstance(value, dict):
        return {key: published_paths(item, source_root, destination_root) for key, item in value.items()}
    if isinstance(value, list):
        return [published_paths(item, source_root, destination_root) for item in value]
    if isinstance(value, tuple) and type(value) is tuple:
        return tuple(published_paths(item, source_root, destination_root) for item in value)
    if is_dataclass(value) and not isinstance(value, type):
        return replace(value, **{f.name: published_paths(getattr(value, f.name), source_root, destination_root)
                                 for f in fields(value) if f.init})
    return value


def accept_result(result, *, model, state, destination, name, stages=None, error=None):
    """Apply records, bindings, diagnostics and artifacts with existing checkpoints."""
    source_root, destination = result.artifact_root.resolve(), destination.resolve()
    records = {path: record for path, record in result.records.items() if path}
    model.validate_merge(records, owner=name)
    attributes = StageJournal.reconcile(model, result.records[""].attrs, allowed_scope=name)
    changes = published_paths(result.state, source_root, destination)
    artifacts, targets = [], set()
    for path in result.artifacts:
        source = Path(path).resolve()
        if not source.is_relative_to(source_root) or not source.is_file():
            raise ValueError(f"Invalid worker artifact: {source}")
        target = (destination / source.relative_to(source_root)).resolve()
        if not target.is_relative_to(destination):
            raise ValueError(f"Invalid artifact destination: {target}")
        if target.exists() or target in targets:
            raise FileExistsError(f"Pipeline artifact already exists: {target}")
        targets.add(target)
        artifacts.append((source, target))
    published, previous_state = [], dict(state)
    try:
        with model.checkpoint():
            with model.writing_as(name):
                model.merge(records)
            StageJournal.apply(model, {key: value for key, value in attributes.items() if key in StageJournal.keys})
            write_root_attributes(model, {key: value for key, value in attributes.items() if key not in StageJournal.keys})
            for stage, status in (stages or {}).items():
                StageJournal.observe(model, name, stage, status)
            StageJournal.finalize(model, name, error)
            for source, target in artifacts:
                target.parent.mkdir(parents=True, exist_ok=True)
                incoming = target.with_name(f".ef-{uuid4().hex[:16]}{target.suffix}")
                try:
                    shutil.copyfile(source, incoming)
                    if target.exists():
                        raise FileExistsError(f"Pipeline artifact already exists: {target}")
                    incoming.rename(target)
                    published.append(target)
                finally:
                    incoming.unlink(missing_ok=True)
            state.update(changes)
    except BaseException:
        state.clear()
        state.update(previous_state)
        for path in published:
            path.unlink(missing_ok=True)
        raise
    for message in result.messages:
        try:
            Logger.log(message.replace(str(source_root), str(destination)))
        except Exception:
            # Presentation callbacks cannot invalidate a committed transaction.
            pass


class ProductAcceptance:
    """Called only by the coordinator; the receiver thread only acknowledges receipt."""

    def __init__(self, *, model, state, destination, inputs, registry):
        self.model, self.state, self.destination = model, state, destination
        self.inputs, self.registry = inputs, registry
        self.data = DataManager(model.storage)
        self._thread = get_ident()
        self._commits = {}

    def _check_thread(self):
        if get_ident() != self._thread:
            raise RuntimeError("Authoritative acceptance belongs to the coordinator thread.")

    def accept_event(self, event, invocation_id):
        self._check_thread()
        started = perf_counter()
        path = Path(event.result_path).resolve()
        if not path.is_relative_to(self.model.storage.directory.resolve()) or not path.is_file():
            raise ValueError("Invalid worker commit manifest.")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024**2), b""):
                digest.update(chunk)
        fingerprint = digest.digest()
        key = invocation_id, event.commit_id
        if not event.commit_id:
            raise ValueError("A worker commit needs a unique commit ID.")
        if key in self._commits:
            if self._commits[key] != fingerprint:
                raise ValueError("Conflicting duplicate commit delivery.")
            return False
        result = load_result(path, self.model.storage, self.inputs)
        if result.invocation_id != invocation_id or result.commit_id != event.commit_id:
            raise ValueError("Worker commit producer/identity mismatch.")
        if event.handle is not None:
            if (len(result.products) != 1 or result.products[0][0] != event.product
                    or pickle.dumps(result.products[0][1], protocol=5) != pickle.dumps(event.handle, protocol=5)):
                raise ValueError("Product event and committed handle disagree.")
        products = []
        for name, handle in result.products:
            product_id = ProductId(invocation_id, name)
            self.data.acquire(handle)
            self.registry.validate_accept(product_id, handle, result.commit_id)
            products.append((product_id, handle))
        for record in result.records.values():
            if record.kind == "dataset":
                self.data.validate_storage(record.value)
        self.data.validate_storage(result.state)
        accept_result(result, model=self.model, state=self.state, destination=self.destination,
                      name=invocation_id.label)
        # All potentially failing validation/I/O has completed. Registry mutations
        # only install prevalidated metadata and cannot invoke user callbacks.
        for product_id, handle in products:
            self.registry.accept(product_id, handle, result.commit_id)
        self._commits[key] = fingerprint
        metrics = self.model.storage.metrics
        for name, value in (result.metrics or {}).items():
            metrics.add("worker_" + name, value)
        metrics.add("worker_serializations", 1)
        metrics.add("worker_manifest_bytes", path.stat().st_size)
        metrics.add("accepted_commits", 1)
        metrics.add("acceptance_seconds", perf_counter() - started)
        if event.emitted_at is not None:
            metrics.add("emission_to_acceptance_seconds", perf_counter() - event.emitted_at)
            if products:
                metrics.add("product_availability_seconds", perf_counter() - event.emitted_at)
                metrics.add("exposed_products_accepted", len(products))
        return True

    def accept_local(self, invocation_id, stage, name, handle):
        self._check_thread()
        self.data.acquire(handle)
        self.registry.accept(ProductId(invocation_id, name), handle, uuid4().hex)

    def observe(self, invocation_id, event):
        self._check_thread()
        if event.journal is not None:
            if set(event.journal) - set(StageJournal.keys):
                raise ValueError("A graph diagnostic event may contain only stage journal records.")
            StageJournal.apply(self.model, StageJournal.reconcile(self.model, event.journal, allowed_scope=invocation_id.label))
        status = {"StageStarted": "running", "StageFinished": "succeeded",
                  "StageFailed": "failed", "StageSkipped": "skipped"}.get(event.kind)
        if status:
            StageJournal.observe(self.model, invocation_id.label, event.stage, status, event.error)

    def finalize(self, invocation_id, error=None):
        self._check_thread()
        StageJournal.finalize(self.model, invocation_id.label, error)
        self.registry.finalize_producer(invocation_id, error or "Producer did not supply product.")
