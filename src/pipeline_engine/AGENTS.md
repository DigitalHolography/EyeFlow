# Pipeline engine map

This package owns pipeline descriptors, DAG resolution, shared run validation,
contexts, scheduling, checkpoints, resource admission, and execution. It does
not own scientific formulas, pipeline-specific HDF5 paths, Tkinter widgets,
or persisted-setting policy.

| Task | Start with |
|---|---|
| Registration metadata and options | `planning/base.py`; pipeline declarations in `src/pipelines/` |
| Dependencies and target closure | `planning/dag.py` |
| GUI/CLI validation and batch policy | `execution/run_service.py` |
| Context, state, and output access | `context.py`, `io/` |
| Invocation lifecycle and statuses | `execution/runtime.py`, `execution/invocation.py` |
| Product-level stage graphs | `stage_graph/`; [stage graphs](../../docs/pipeline_stage_graphs.md) |
| Process work and input replay | `execution/scheduler.py`, `execution/worker.py`, `io/input/`; [parallel execution](../../docs/parallel_pipeline_execution.md) |
| HDF5 model, scratch, resources | `output/output_model.py`, `data/storage.py`, `resources.py` |

`requires` names importable packages; `dag_requires` names produced data.
Options close their own prerequisites before the DAG is resolved. A pipeline
implicitly produces its name; a key with no producer is external and must be
validated by its consumer. `auto_run_after` adds a companion; `run_once` keeps
shared work outside velocity fan-out.

`resolve_run_spec` validates selections, layouts, calibration, and destinations
before execution. `execute_run` replaces the EyeFlow-owned output directory per
attempt. Declared execution variants isolate state and output namespaces;
state consumers need a declared DAG producer. Accepted stage handoffs survive
later stage failure. Keep GUI and CLI on this shared run path.

Use focused DAG, run-service, graph, or coordinator tests for the behavior that
changes. See [test guidance](../../test/AGENTS.md).
