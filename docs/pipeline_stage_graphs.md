# Pipeline-owned stage graphs

A pipeline can return a `StageGraph` containing `StageSpec` nodes. The pipeline
package owns calculations, `ProductSpec` declarations, packers, and exports;
`pipeline_engine/stage_graph/` validates and executes the graph. A simple
pipeline may still return its normal result directly.

## Products and dependencies

A stage publishes each declared product with `ctx.dag.handoff(name, value)`.
A `ProductSpec` can name one HDF5 path or use a packer that returns the usual
path-to-payload mapping. Handoffs snapshot values and atomically commit their
datasets. Duplicate handoffs and overlapping output paths are errors. Managed
graph stages must not bypass handoff for HDF5 writes.

| Declaration | Meaning |
|---|---|
| `requires` | Product must have been handed off |
| `optional_requires` | Wait for the producer, then read with `ctx.dag.get` if available |
| `after` | Named stage must finish successfully |
| `ProductSpec(optional=True)` | Producer may omit the product |
| `ProductSpec(exposed=True)` | Product may cross a pipeline/process boundary |
| `ctx.dag.skip(reason)` | Expected missing input, recorded without a failure |

The planner rejects missing producers, duplicate names, and cycles. Consumers
of an already committed product can continue after its producer later fails;
consumers of a missing required product skip. Independent stages continue.
Stage statuses, failures, skips, timings, and accepted product paths are
recorded in the result.

## Artifacts, processes, and variants

`artifact=True` stages become eligible as soon as their declared products are
ready. Shared figure writers snapshot figures and render them through the
background artifact pool. `wait_for_artifacts=True` creates a barrier for a
consumer such as the report pipeline. Output writers register artifact
ownership so rollback removes only files owned by the failing invocation.

Cross-pipeline values use `ProductRequirement` and an exposed `ProductSpec`;
the product exchange keeps identity and lifetime explicit across process work.
Each velocity workflow has its own context, product store, status, HDF5 root,
and artifact namespace. Graph handles remain invocation-local, and compute
threads return results to the owning pipeline thread before handoff.

See `src/pipeline_engine/stage_graph/`,
`src/pipeline_engine/product_exchange/`, and `test/test_pipeline_graphs.py` for
the current API and focused coverage.
