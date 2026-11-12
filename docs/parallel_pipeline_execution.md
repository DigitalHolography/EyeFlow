# Coordinator and process execution

`pipeline_engine.runtime` owns the resolved run, shared input service, output
model, resource budget, accepted products, and final publication. GUI and CLI
enter through `run_service`.

## Scheduling and inputs

Pipeline declarations may set `execution_mode="process"`; coordinator execution
is the default. Dependency-ready work can overlap across processes and the
coordinator. Process entry points must be importable registered pipelines.
The current process path starts with displacement registration and its later
profile/map analysis. Variant descendants retain coordinator execution.
`pipeline_engine/scheduler.py`, `worker.py`, and `protocol.py` own dispatch
and bounded event queues.

The coordinator opens HD/DV sources once. Workers request data through the
input service rather than reopening source HDF5 files. Requested numeric
datasets are staged in bounded reads to local, read-only `.npy` files and
reused; unrequested datasets are not copied. `EYEFLOW_SCRATCH_DIR` may select
a local scratch directory. Network scratch is rejected, and exhaustion fails
explicitly. Input mappings and worker files remain alive until their consumers
finish.

## Commits, failures, and resources

`OutputModel` holds HDF5 groups, attributes, and payload references before
final serialization. Managed checkpoints journal changed records only.
Stage handoffs commit their declared products atomically; later failures keep
accepted products and discard uncommitted changes. State checkpoints track
bindings, so committed mutable objects must not be edited in place. A worker
publishes closed local snapshots by acknowledged queue event; an abrupt exit
retains its last acknowledged snapshot. Unpublished work cannot be recovered.

A shared batch budget admits compute, input storage, retained arrays, and
rendering. It considers CPU availability and RAM reservations, with local spill
for large retained arrays. Admission estimates are not OS-enforced memory
limits. `resources.py` owns current thresholds; measure representative input
and scratch usage before tuning them.

Reports wait for required artifact jobs and read the assembled output view.
The coordinator writes the accepted HDF5 locally, copies to a destination-side
staging file, and renames it to the final path. Incomplete runs can publish
accepted products with failure status. Analysis of the next input may overlap
the preceding input's export and publication; batch return waits for all
started inputs.

For per-product dependencies and artifact ownership, see
[stage graphs](pipeline_stage_graphs.md). Tests for process spawn, input replay,
crash recovery, and final publication start in `test/test_parallel_coordinator.py`.
