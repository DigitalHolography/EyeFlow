# Contributing to EyeFlow

Start with [AGENTS.md](AGENTS.md) to find the owning subsystem, then read its
local map. [Architecture](docs/architecture.md) describes the boundaries.
The nearest pipeline declaration and its focused tests are the current source
of truth for an implementation.

## Set up and check a change

From the repository root in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .
pip install pytest
python -m pytest test/test_pipeline_discovery.py
```

Install `.[pipelines]` only if the affected pipeline needs optional packages.
Use the smallest relevant test first. The [test guide](test/AGENTS.md)
explains when broader coverage is useful.

## Add or change a pipeline

1. Find a nearby production pipeline in `src/pipelines/`. A small pipeline can
   be one module; a larger one normally has a light `__init__.py` declaration,
   a runner, and output packers. Catalog discovery imports top-level packages,
   so there is no registry list to edit.
2. Declare dependencies and options on `@pipeline` or `@registerPipeline`.
   `requires` names importable Python packages; `dag_requires` names produced
   data. `PipelineOption.requires` names options in the same pipeline;
   `PipelineOption.dag_requires` adds upstream data dependencies. A pipeline
   implicitly produces its own name. A key with no producer is considered
   external and must be validated by its consumer.
3. Use typed HD/DV adapters and existing source loaders. Keep volumes lazy
   until bounded slices are needed. Canonical alignment belongs in
   `input_output/readers/canonical.py`; external aliases belong in the schemas.
   See [data contracts](docs/data-contracts.md) for axes and units.
4. Put reusable numerical work in `src/calculations/`. A pipeline owns its
   product-specific orchestration, state keys, HDF5 packing, and artifacts.
   A `ctx.state` consumer needs a declared DAG producer; stable results belong
   in `ctx.output.h5`, with `unit`, `dimDesc`, and provenance.
5. Use `EyeFlowOutputPaths.active()` for shared output families. Before adding
   or changing a path, search its producer and consumers. Use the shared
   artifact writers and output manager for sidecars.
6. Verify the declaration and the changed contract with focused tests. Check
   both velocity workflows if the product consumes physical velocity.

Hidden pipelines provide shared work through declared dependencies or
`auto_run_after`. Use `run_once=True` only when a product must stay outside
velocity workflow fan-out. The engine owns scheduling and execution variants;
controllers and pipelines do not create a second run path.

Velocity from moments and calibrated frequency bands is physical `mm/s`.
Preserve method, calibration, unit, and timing provenance in derived products;
`pipelines.velocity.semantics` resolves method-specific interpretation. For
product-level handoffs or process execution, see
[stage graphs](docs/pipeline_stage_graphs.md) and
[parallel execution](docs/parallel_pipeline_execution.md).
