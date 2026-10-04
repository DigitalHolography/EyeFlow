# Pipeline engine context

## Boundary

This package owns pipeline metadata, dependency planning, run validation, the
runtime context, and execution. It does not own scientific calculations,
pipeline-specific HDF5 paths, GUI widgets, or persisted-setting policy.

Read [repository architecture](../../docs/architecture.md) first only when the
task crosses front-end or I/O boundaries.

## Key modules

| Module | Responsibility |
|---|---|
| `base.py` | `PipelineDescriptor`, `PipelineOption`, result payloads, decorators, and the process interface |
| `dag.py` | Semantic-key producer graph, option-dependent requirements, topological order, and target closure |
| `run_service.py` | Shared GUI/CLI run-spec validation, input expansion, destination mapping, and batch execution |
| `runtime.py` | Opens inputs/output, initializes provenance, creates contexts, invokes descriptors, and records completion |
| `context.py` | Namespaced inputs/output, shared state, selected options/order/method, result persistence |
| `errors.py` | User-facing pipeline exception formatting |
| `imports.py` | Convenience import surface used by pipeline declarations |

Pipeline discovery itself is in `src/pipelines/__init__.py`; settings
normalization is in `src/app_settings.py`.

## Dependency semantics

- Import requirements (`requires`) only control availability.
- Data requirements (`dag_requires`) are matched to one producer's pipeline
  name or `dag_produces` key.
- Every pipeline implicitly produces its own name.
- Multiple producers for a key and dependency cycles fail graph construction.
- A key with no producer is considered external and creates no edge.
- Selected option `requires` are closed transitively within one pipeline before
  the DAG is resolved. Selected option `dag_requires` add producer edges.
- Hidden pipelines may enter a plan through dependencies; only visible
  available pipelines are valid user targets.
- Descriptor discovery order breaks ties between otherwise independent nodes.

When modifying scheduling, inspect `base.py`, `dag.py`, the affected pipeline
declarations, `test_pipeline_discovery.py`,
`test_pipeline_library_dependencies.py`, and the scheduling cases in
`test_run_service.py`. Calculation implementations are normally irrelevant.

## Run lifecycle

`resolve_run_spec` validates the velocity method, filters selectable targets,
resolves option closure and the DAG, rejects unavailable required descriptors,
resolves every `.holo` companion layout, maps outputs, and rejects duplicate
destinations. No files are executed at this stage.

`execute_run` iterates `RunRequest`s. Before each attempt it replaces the
EyeFlow output directory, delegates to `run_pipelines_to_output`, records an
output or `RunFailure`, and checks cancellation only between files.

The runtime opens HD/DV sources and the output HDF5 for the full per-file plan.
It creates a new `PipelineContext` for every descriptor over the same output
handle and shared-state dictionary. Pipeline options, execution order, and
directly selected targets are immutable views for the run, along with the
velocity method. A pipeline returning a mapping or `ProcessResult` is persisted;
a pipeline returning `None` owns its direct writes.

## Invariants

- Keep GUI and CLI behavior in the shared run-service path.
- Resolve and validate all plan-level choices before deleting an output.
- Preserve the distinction between selected targets and full execution order;
  both are written as output provenance.
- A state consumer must have a matching declared DAG dependency. The engine
  does not validate arbitrary `ctx.state` keys.
- `input_slot` controls preferred merged-attribute lookup. Normal `.holo` layout
  resolution still requires both HD and DV companions.
- Do not catch scientific exceptions inside the engine. `_run_pipeline_descriptor`
  wraps them with pipeline identity and retains their cause.
- Do not mutate the global catalog during a run. Reload/discovery belongs before
  plan construction.

## Change scopes

| Change | Usually inspect | Usually irrelevant |
|---|---|---|
| Decorator metadata | `base.py`, affected declarations, discovery tests | calculations, UI views |
| DAG behavior | `dag.py`, option declarations, dependency tests | HDF5 writers, figures |
| Run input/output routing | `run_service.py`, `holo_run_layout.py`, run-service tests | scientific formulas |
| Context API/result handling | `context.py`, `runtime.py`, context and writer tests | GUI layout |
| Progress/cancellation | `run_service.py`, `ui/controllers/run.py`, run-service tests | pipeline algorithms unless callback frequency changes |
