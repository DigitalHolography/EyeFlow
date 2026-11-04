"""Reproducible synthetic ownership counters, not a representative NAS benchmark.

Run with the editable repository interpreter. --processes additionally measures
the real A-emits/B-consumes/A-crashes fixture (two workers and one sync manager).
"""

import argparse
import json
import multiprocessing
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

from input_output.output_manager import OutputManager, OutputType
from input_output.runs.layout import HoloRunLayout
from pipeline_engine.context import PipelineContext
from pipeline_engine.invocation import _pipeline_checkpoint
from pipeline_engine.output_model import OutputModel
from pipeline_engine.protocol import InvocationResult
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec
from pipeline_engine.stage_graph.executor import StageGraphExecutor
from pipeline_engine.storage import ProductStorage
from pipeline_engine.transfer import dump_result


def local(array_mib, consumers):
    budget = ResourceBudget(cpu_threads=4, memory_bytes=2 * 1024**3)
    with TemporaryDirectory() as temp, budget.activate(), ProductStorage(budget=budget) as storage:
        manager = OutputManager(HoloRunLayout.from_holo(Path(temp) / "scan.holo"))
        model = OutputModel(storage, Path(temp) / "final.h5")
        ctx = PipelineContext(work_h5=model, holodoppler_h5=None, doppler_vision_h5=None,
                              output_manager=manager, pipeline_name="bench")
        count = array_mib * 1024**2 // 4

        def produce(ctx):
            ctx.dag.handoff("x", np.arange(count, dtype=np.float32))

        addresses = []
        def consume(ctx):
            value = ctx.dag.require("x")
            assert value[-1] == np.float32(count - 1)
            addresses.append(value.ctypes.data)
            budget.refresh()

        graph = StageGraph((StageSpec("producer", produce, products=(ProductSpec("x", "Bench/x"),)),
                            *(StageSpec(f"consumer{i}", consume, requires=("x",)) for i in range(consumers))))
        with _pipeline_checkpoint(model, ctx.state.raw, manager) as checkpoint:
            StageGraphExecutor(graph, ctx, checkpoint, _pipeline_checkpoint).run()
        assert len(set(addresses)) == 1, "Consumers received different array buffers"
        cursor, sizes = {}, []
        for _ in range(3):
            records, cursor = model.committed_delta(cursor)
            path = dump_result(InvocationResult(records, {}, (), Path(temp)), storage)
            sizes.append(Path(path).stat().st_size)
        model.publish()
        budget.refresh()
        return {"array_bytes": count * 4, "consumers": consumers,
                "unique_consumer_buffers": len(set(addresses)), "manifest_bytes": sizes,
                "peak_process_tree_rss": budget.peak_rss, "metrics": storage.metrics.snapshot()}


def processes():
    from dataclasses import replace
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "test"))
    from product_exchange_fixtures import PRODUCER, CONSUMER
    from pipeline_engine.base import PipelineDescriptor
    from pipeline_engine.runtime import PipelineRunError, run_pipelines_to_output

    mp = multiprocessing.get_context("spawn")
    with TemporaryDirectory() as temp, mp.Manager() as sync_manager:
        output = OutputManager(HoloRunLayout.from_holo(Path(temp) / "scan.holo"))
        sync = {name: sync_manager.Event() for name in ("received", "crash", "b_started", "b_finished")}
        sync.update(pid=sync_manager.Value("i", 0), backing=sync_manager.Value(str, ""))

        def initial(ctx):
            ctx.state.set("sync", sync)
            ctx.state.set("final_path", str(output.path_for(OutputType.H5)))

        descriptors = [PipelineDescriptor("initial", "initial", True,
            pipeline_factory=lambda: SimpleNamespace(name="initial", run=initial)),
            replace(PRODUCER, dag_requires=("initial",)), replace(CONSUMER, dag_requires=("initial",))]
        budget = ResourceBudget(cpu_threads=8, memory_bytes=2 * 1024**3)
        pending = []
        with budget.activate():
            try:
                run_pipelines_to_output(output_manager=output, pipelines=descriptors,
                    holodoppler_h5=None, doppler_vision_h5=None, copy_hd_passthrough=False,
                    on_finalization_ready=pending.append)
            except PipelineRunError:
                pass
            assert sync["b_finished"].is_set()
            retained = pending[0]
            retained.finish()
            budget.refresh()
            return {"peak_process_tree_rss": budget.peak_rss,
                    "metrics": retained.model.storage.metrics.snapshot(),
                    "cpu_slots_after_finish": budget.cpu_used,
                    "reserved_bytes_after_finish": sum(budget.used.values())}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--array-mib", type=int, default=16)
    parser.add_argument("--consumers", type=int, default=3)
    parser.add_argument("--processes", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"local": local(args.array_mib, args.consumers)}
    if args.processes:
        result["processes"] = processes()
    rendered = json.dumps(result, indent=2)
    if args.output is not None:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
