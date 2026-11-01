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
| `base.py` | `PipelineDescriptor`, `PipelineOption`, `ExecutionVariant`, result payloads, decorators, and the process interface |
| `dag.py` | Semantic-key producer graph, option-dependent requirements, topological order, and target closure |
| `run_service.py` | Shared GUI/CLI run-spec validation and batch execution; delegates input expansion to `input_output/runs/inputs.py` and destination mapping to `input_output/output_manager.py` |
| `input_access.py` | Context-aware timing and attribute conveniences for pipeline runners |
| `runtime.py`, `invocation.py` | Per-input workflow coordination / one descriptor's execution and checkpoints |
| `scheduler.py` | Dependency readiness, artifact barriers and progress order |
| `process_executor.py`, `worker.py`, `protocol.py` | Worker lifetime, transactional result acceptance and acknowledged messages |
| `input_service.py`, `storage.py`, `transfer.py` | Lazy source replay, immutable payloads and local-file transport |
| `output_model.py`, `finalization.py` | Output tree/journals and retained payload ownership through publication |
| `resources.py`, `src/runtime_limits.py` | Batch admission/memory accounting and process/thread allowances |
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
- `auto_run_after` adds a companion when its named predecessor is scheduled,
  adds an ordering/dependency edge, and preserves directly selected targets.
- Hidden pipelines may enter a plan through dependencies or automatic
  companions; only visible available pipelines are valid user targets.
- Descriptor discovery order breaks ties between otherwise independent nodes.

When modifying scheduling, inspect `base.py`, `dag.py`, the affected pipeline
declarations, `test_pipeline_discovery.py`,
`test_pipeline_library_dependencies.py`, and the scheduling cases in
`test_run_service.py`. Calculation implementations are normally irrelevant.

## Run lifecycle

`resolve_run_spec` validates band-ratio calibration, filters selectable targets,
resolves option closure and the DAG, rejects unavailable required descriptors,
resolves every `.holo` companion layout, maps outputs, and rejects duplicate
destinations. No files are executed at this stage.

`execute_run` iterates `RunRequest`s. Before each attempt it replaces the
EyeFlow output directory, delegates to `run_pipelines_to_output`, records an
output or `RunFailure`, and checks cancellation only between files.

The [execution design](../../docs/parallel_pipeline_execution.md) is the contract
for scheduling, transactional acceptance, resources and storage lifetime.
Base displacement registration and its profile/map companion execute in separate
run-once processes; the companion waits for topology and shared timing. Worker
logs are forwarded by the result receiver while coordinator work runs.
Change the component that owns the behavior; keep transport out of the scheduler
and scientific work out of the executor. HDF5 serialization belongs in
`input_output/writers/h5_publication.py`; local scratch policy in
`src/scratch_location.py`.

`ctx.computation_cache` is a per-input memo shared across execution variants,
separate from transactional `ctx.state`. Use it only for immutable successful
pure calculations with explicit source/settings compatibility keys. Never cache
failures, beat-specific products, or output side effects. A later packing/write
failure does not discard a valid cached calculation; a failed calculation must
remain independently retryable. Each input execution starts with an empty memo.

## Invariants

- Keep GUI and CLI behavior in the shared run-service path.
- Resolve and validate all plan-level choices before deleting an output.
- Preserve the distinction between selected targets and full execution order;
  both are written as output provenance.
- A state consumer must have a matching declared DAG dependency. The engine
  does not validate arbitrary `ctx.state` keys.
- `input_slot` controls preferred merged-attribute lookup. Normal `.holo` layout
  resolution still requires both HD and DV companions.
- Descriptor failures roll back that invocation's dataset replacements,
  attributes, state bindings, and newly created artifacts. Earlier pipeline
  outputs remain available. Skip consumers of failed/skipped active dependencies
  while continuing independent pipelines, including within a failed variant.
  Shared errors or failure of every variant report an incomplete run only after
  safe pipelines finish; its partial HDF5 path remains in `RunResult.outputs`.
  Root `pipeline_status`, `pipeline_failures`, and `pipeline_skips` record outcomes.
  Writers replace datasets; journals do not cover in-place array or
  mutable state-object edits, or overwrites of existing sidecar files.
  `invocation._run_pipeline_descriptor` wraps them with pipeline identity and retains
  their cause.
- Variant producers must declare `produces_execution_variants=True`; publishing
  variants without the declaration, or declaring them without publishing any,
  is an execution error. Nested variant producers are not supported.
- Do not mutate the global catalog during a run. Reload/discovery belongs before
  plan construction.

## Change scopes

| Change | Usually inspect | Usually irrelevant |
|---|---|---|
| Decorator metadata | `base.py`, affected declarations, discovery tests | calculations, UI views |
| DAG behavior | `dag.py`, option declarations, dependency tests | HDF5 writers, figures |
| Run input/output routing | `run_service.py`, `input_output/runs/layout.py`, run-service tests | scientific formulas |
| Context API/result handling | `context.py`, `runtime.py`, context and writer tests | GUI layout |
| Progress/cancellation | `run_service.py`, `ui/controllers/run.py`, run-service tests | pipeline algorithms unless callback frequency changes |
