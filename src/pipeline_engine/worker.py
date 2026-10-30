"""Spawn-safe pipeline execution with private state, output and artifact staging."""

import os
import traceback
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path

from .protocol import _SINK, WorkerEvent


@dataclass
class InvocationSpec:
    name: str
    module: str
    directory: str
    output_manager: object
    input_metadata: dict
    input_client: object
    settings: dict
    initial_result: str
    threads: int
    memory_account: object = None
    acknowledgements: object = None
    native_parallelism: bool = False


def run_worker(spec, events):
    # Environment defaults precede scientific imports in the spawned process.
    from .resources import _THREADS, configure_native_threads

    configure_native_threads(spec.threads if spec.native_parallelism else 1)
    _THREADS.set(1 if spec.native_parallelism else spec.threads)
    os.environ["EYEFLOW_WORKER_DIRECTORY"] = spec.directory
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
    from .runtime import _run_pipeline_descriptor
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
    messages = []

    def snapshot():
        return dump_result(
            {
                "records": output.records,
                "state": {
                    key: value
                    for key, value in variables.items()
                    if key not in initial_state or value is not initial_state[key]
                },
                "artifacts": ownership.paths,
                "manager": manager,
                "messages": messages,
            },
            storage,
        )

    def send(kind, stage=None, detail=None):
        detail = dict(detail or {})
        if kind == "ProductCommit":
            detail["result"] = snapshot()
        events.put(WorkerEvent(kind, spec.name, stage, detail))
        # Queue.put only feeds a process-local buffer. Do not enter native
        # calculation (which can crash) until this event is fully received.
        spec.acknowledgements.get()

    token = _SINK.set(send)
    error = None
    try:
        send("WorkerStarted")
        import_module(spec.module)
        descriptor = PIPELINE_REGISTRY[spec.name]
        initial_state.update(load_result(spec.initial_result, storage, files))
        variables.update(initial_state)
        settings = dict(spec.settings)
        configs = settings.pop("configs")
        inputs = PipelineInputs(
            **{
                slot: PipelineInputSource(
                    RawH5SourceReader(h5file=files[slot], label=slot.upper()), configs[slot]
                )
                for slot in ("hd", "dv")
            }
        )
        with (
            Logger.capture_messages(messages.append),
            artifact_transaction(ArtifactScope(ownership, manager.layout.ef_dir)),
        ):
            _run_pipeline_descriptor(
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
                **settings,
            )
    except BaseException:  # noqa: BLE001 - report partial work at the process boundary
        error = traceback.format_exc()
    finally:
        try:
            send("InvocationFinished", detail={"error": error, "result": snapshot()})
        finally:
            _SINK.reset(token)
            spec.input_client.close()
            scratch.close()
            storage.close()
