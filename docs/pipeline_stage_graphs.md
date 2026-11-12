# Pipeline-owned stage graphs

A pipeline may return a `PipelineGraph` of `PipelineNode`s instead of one
result mapping. The pipeline package owns calculations, product definitions,
packers, and exports; the engine executes the graph. Legacy function pipelines
remain supported. `src/pipelines/low_freq_spectra/stages.py` is the production
example.

## Products and dependencies

A node publishes a declared `PipelineProduct` with `ctx.dag.handoff(name,
value)`. A product can name one HDF5 path or use a packer returning standard
`DatasetValue` or `(data, attrs)` payloads. Handoffs validate ownership,
snapshot the value, and commit its datasets atomically. Duplicate product
handoffs or overlapping paths are errors. Graph stages must use handoff for
managed HDF5 writes; direct managed writes are rejected. Keep large source
videos in typed input adapters rather than handing them between nodes.

| Declaration | Meaning |
|---|---|
| `requires` | Product must have been handed off |
| `optional_requires` | Order after producer; read with `ctx.dag.get` if available |
| `after` | Named stage must finish successfully |
| `PipelineProduct(optional=True)` | Stage may omit the product |
| `ctx.dag.skip(reason)` | Expected missing input, recorded without a failure |

The planner checks missing producers, duplicate names, and cycles. Consumers
of a committed product can continue after its producer later fails; consumers
of missing required products skip. Independent nodes continue. External
pipeline dependencies still require whole-pipeline success. Root
`stage_status`, `stage_failures`, `stage_skips`, and `stage_products` attributes
record the outcome and accepted product paths per workflow.

## Artifacts and variants

`artifact=True` nodes may start when their inputs are committed. Managed
figure writers can defer rendering to isolated processes using immutable
snapshots. Reports use `wait_for_artifacts=True` and appropriate graph
dependencies. OutputManager and writers record explicit artifact ownership;
rollback removes only uncommitted files owned by that invocation. Custom
writers must register their paths. A failed export retains numerical results
and makes the run incomplete.

Each velocity workflow has its own context, products, status, and artifact
namespace. The HDF5 adapter maps logical moments paths to band-ratio paths for
that variant. Graph handles are invocation-local; compute threads return
results to the owning pipeline thread for handoff.

See `src/pipeline_engine/graphs.py` and `test/test_pipeline_graphs.py` for the
current API and focused coverage.
