"""Spawn-safe pipeline execution with private state, output and artifact staging."""

import os
import traceback
from importlib import import_module
from pathlib import Path
from queue import Empty
from time import perf_counter
from uuid import uuid4

from .protocol import _SINK, InvocationResult, ProductAcceptedEvent, TransportReceipt, WorkerEvent


class _WorkerTransportError(BaseException):
    """A broken sender ends the invocation, rather than a scientific stage."""


def run_worker(spec, events):
    # Environment defaults precede scientific imports in the spawned process.
    from runtime_limits import (
        configure_native_threads, configure_numeric_threads, set_thread_allowance,
    )

    native_threads = spec.threads if spec.native_parallelism else 1
    configure_numeric_threads(native_threads)
    set_thread_allowance(1 if spec.native_parallelism else spec.threads)
    os.environ["EYEFLOW_WORKER_DIRECTORY"] = str(spec.directory)
    from dataclasses import replace

    from input_output.artifact_ownership import (
        ArtifactOwnership,
        ArtifactScope,
        artifact_transaction,
    )
    from input_output.readers.h5 import PipelineInputSource, RawH5SourceReader
    from utils.logger import Logger

    from .base import PIPELINE_REGISTRY
    from .context import PipelineInputs
    from .input_service import InputFile
    from .output_model import OutputModel
    from .invocation import _run_pipeline_descriptor
    from .scratch import RunArrayStorage
    from .storage import ProductStorage
    from .transfer import dump_result, load_result

    storage = ProductStorage(
        directory=Path(spec.directory) / "products", budget=spec.memory_account
    )
    scratch = RunArrayStorage(
        Path(spec.directory), memory_budget_bytes=spec.memory_account.memory_limit
    )
    scratch.budget = spec.memory_account
    files = {
        slot: InputFile(slot, *metadata, spec.input_client) if metadata else None
        for slot, metadata in spec.input_metadata.items()
    }
    variables = {}
    initial_state = {}
    output = OutputModel(storage, Path(spec.directory) / "private.h5")
    manager = replace(
        spec.output_manager,
        layout=replace(spec.output_manager.layout, _root_dir=Path(spec.directory) / "artifacts"),
    )
    ownership = ArtifactOwnership()
    published_artifact_root = spec.output_manager.layout.ef_dir.resolve()

    cursor, sent_state, sent_artifacts = {}, {}, set()
    sent_metrics = {}
    transport_failed = False

    def snapshot(commit_id, products, *, force=False):
        records, next_cursor = output.committed_delta(cursor)
        changes = {key: value for key, value in variables.items()
                   if (value is not sent_state[key] if key in sent_state
                       else key not in initial_state or value is not initial_state[key])}
        artifacts = tuple(path for path in ownership.paths if path not in sent_artifacts)
        metrics = {key: value for key, value in storage.metrics.snapshot().items()
                   if key not in {"serializations", "manifest_bytes"}}
        delta_metrics = {key: value - sent_metrics.get(key, 0) for key, value in metrics.items()}
        candidate = (next_cursor, {**sent_state, **changes}, sent_artifacts | set(artifacts), metrics)
        # Stage failures/skips have status-only events. Stage completion needs a
        # manifest only for real post-handoff records, bindings or artifacts.
        needed = bool(products or changes or artifacts or any(records.keys() - {""}))
        path = dump_result(
            InvocationResult(
                records=records,
                state=changes,
                artifacts=artifacts,
                artifact_root=manager.layout.ef_dir,
                invocation_id=spec.invocation_id,
                commit_id=commit_id,
                products=products,
                metrics=delta_metrics,
            ),
            storage,
        ) if needed or force else None
        if path is None:
            candidate = (*candidate[:3], sent_metrics)
        return path, candidate

    def send(kind, stage=None, detail=None):
        nonlocal cursor, sent_state, sent_artifacts, sent_metrics, transport_failed
        if transport_failed:
            raise _WorkerTransportError("Worker transport has already failed.")
        detail = dict(detail or {})
        detail["emitted_at"] = perf_counter()
        detail["invocation_id"] = spec.invocation_id
        detail["commit_id"] = uuid4().hex
        candidate = None
        try:
            if kind in {"ProductCommit", "StageFinished", "InvocationFinished"}:
                handle = detail.get("handle")
                products = ((detail["product"], handle),) if handle is not None else ()
                # Final completion includes runtime root metadata. Status-only
                # stage events do not need another result serialization.
                path, candidate = snapshot(detail["commit_id"], products, force=kind == "InvocationFinished")
                if path is not None:
                    detail["result_path"] = path
                else:
                    from .stage_graph.journal import StageJournal
                    detail["journal"] = {key: output.attrs[key] for key in StageJournal.keys if key in output.attrs}
            event_type = ProductAcceptedEvent if kind == "ProductCommit" else WorkerEvent
            events.put(event_type(kind, spec.name, stage, **detail))
            receipt = spec.acknowledgements.get(timeout=30)
            if not isinstance(receipt, TransportReceipt) or receipt.commit_id != detail["commit_id"]:
                raise RuntimeError("Invalid worker transport receipt.")
        except Exception as exc:
            transport_failed = True
            reason = "Transport receipt timed out." if isinstance(exc, Empty) else str(exc)
            raise _WorkerTransportError(f"Worker transport failed: {reason}") from exc
        # Only a matching receipt advances the acknowledged delivery boundary.
        # Coordinator acceptance is a separate operation over run-owned backing.
        if candidate is not None:
            cursor, sent_state, sent_artifacts, sent_metrics = candidate

    def forward_log(message):
        staged_artifact_root = manager.layout.ef_dir.resolve()
        send(
            "LogMessage",
            detail={
                "message": str(message).replace(
                    str(staged_artifact_root), str(published_artifact_root)
                )
            },
        )

    token = _SINK.set(send)
    error = None
    outcome = None
    try:
        send("WorkerStarted")
        import_module(spec.module)
        configure_native_threads(native_threads)
        descriptor = PIPELINE_REGISTRY[spec.name]
        initial_state.update(load_result(spec.initial_result, storage, files, owner_root=spec.owner_id))
        variables.update(initial_state)
        configs = spec.settings.configs
        inputs = PipelineInputs(
            **{
                slot: PipelineInputSource(
                    RawH5SourceReader(h5file=files[slot], label=slot.upper()), configs[slot]
                )
                for slot in ("hd", "dv")
            }
        )
        with (
            Logger.capture_messages(forward_log),
            artifact_transaction(ArtifactScope(ownership, manager.layout.ef_dir)),
        ):
            outcome = _run_pipeline_descriptor(
                descriptor,
                work_h5=output,
                output_h5_path=Path(output.filename),
                output_manager=manager,
                holodoppler_h5=files["hd"],
                doppler_vision_h5=files["dv"],
                holodoppler_config=configs["hd"],
                doppler_vision_config=configs["dv"],
                shared_inputs=inputs,
                variables=variables,
                scratch=scratch,
                on_pipeline_success=None,
                on_progress=None,
                product_inputs=spec.product_inputs,
                product_owner_id=spec.owner_id,
                **spec.settings.context_arguments(),
            )
    except BaseException:  # noqa: BLE001 - report partial work at the process boundary
        error = traceback.format_exc()
    finally:
        try:
            if transport_failed:
                events.put(WorkerEvent("InvocationFinished", spec.name, error=error,
                                       invocation_id=spec.invocation_id, status="failed"))
            else:
                send("InvocationFinished", detail={"error": error,
                    "status": "failed" if error else outcome.status,
                    "reason": outcome.reason if outcome is not None else None})
        finally:
            _SINK.reset(token)
            from contextlib import ExitStack
            with ExitStack() as cleanup:
                cleanup.callback(storage.close)
                cleanup.callback(scratch.close)
                cleanup.callback(spec.input_client.close)
