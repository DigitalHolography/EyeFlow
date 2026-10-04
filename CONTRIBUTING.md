# Contributing

This guide covers repository setup and the common case of adding or changing a
pipeline. Read [AGENTS.md](AGENTS.md) first for task-specific navigation, then
the local subsystem context. Preserve the boundaries documented in
[architecture](docs/architecture.md).

## Setup and checks

```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .
pip install pytest
python -m pytest
```

Install optional pipeline or GPU dependencies only when the affected code needs
them:

```powershell
pip install -e .[pipelines]
pip install -e .[gpu]
```

Run focused tests before the full suite. The
[production-to-test map](test/AGENTS.md) identifies the right files.

## Pipeline responsibilities

A pipeline is a runtime adapter. It declares dependencies and options, reads
typed inputs or declared upstream state, calls scientific code, and publishes
stable results. Put reusable formulas in `src/calculations/`; put external
HDF5 path logic in `src/input_output/schema/`.

Every pipeline receives a `PipelineContext` with:

- `ctx.inputs.hd` and `ctx.inputs.dv`: raw and typed input adapters;
- `ctx.output.h5`: the shared primary output HDF5 interface;
- `ctx.output`: optional artifact paths/writers;
- `ctx.state`: in-memory state shared by this run's ordered pipelines;
- `ctx.pipeline_options`: resolved selections;
- `ctx.pipeline_order`: the full scheduled order;
- `ctx.pipeline_targets`: pipelines selected directly by the user;
- `ctx.velocity_estimation_method`: validated velocity semantics;
- `ctx.band_ratio_frequency_scale_hz`: validated Hz-per-ratio calibration;
- `ctx.attrs`: merged output/input/config attributes.

Use [the pipeline context](src/pipelines/AGENTS.md) for current shared-state and
scientific dependency flows.

## Add a pipeline

Create a module or package directly below `src/pipelines/` and decorate its
entry point. Discovery is automatic; do not edit a registry list.

```python
import numpy as np

from pipeline_engine import DatasetValue, pipeline


@pipeline(
    name="my_metric",
    description="Compute a mean zeroth-moment signal.",
    requires=["numpy", "h5py"],
    dag_produces=["my_metric"],
    input_slot="hd",
)
def run(ctx):
    ctx.require_inputs("hd")
    moment0 = ctx.inputs.hd.as_holodoppler().moment0_dataset()
    mean_per_frame = np.nanmean(moment0, axis=(1, 2)).astype(np.float32)
    return {
        "Processing/MyMetric/MeanMoment0/value": DatasetValue(
            data=mean_per_frame,
            attrs={"unit": "a.u.", "dimDesc": ["frame"]},
        )
    }
```

`requires` names importable Python packages used to determine availability.
`dag_requires` and `dag_produces` name semantic runtime products. Every
pipeline implicitly produces its pipeline name.

For a larger pipeline, keep the declaration small:

```text
src/pipelines/my_metric/
  __init__.py
  runner.py
  outputs.py
```

`__init__.py` declares metadata and delegates to the runner. Keep source
assembly, orchestration, output packing, and product-specific calculations
separate when that improves testing; do not create layers for trivial code.

Verify discovery:

```powershell
python -c "from pipelines import load_pipeline_catalog; a, m = load_pipeline_catalog(); print([p.name for p in a]); print([p.name for p in m])"
```

## Declare dependencies and options

Use data-oriented DAG keys:

```python
@pipeline(
    name="summarize_signal",
    dag_requires=["prepared_signal"],
    dag_produces=["signal_summary"],
)
def run(ctx):
    signal = ctx.state.get("prepared_signal.state")
    if signal is None:
        raise RuntimeError("Prepared signal state is unavailable.")
    ...
```

The producer must set the state and declare `prepared_signal`. A DAG key with
no producer is treated as external, so consumers must still fail clearly when
its value is absent.

Pipeline options can require other options and add upstream DAG needs:

```python
from pipeline_engine.imports import PipelineOption, pipeline


@pipeline(
    name="my_metric",
    options=[
        PipelineOption("per_beat", "Per beat"),
        PipelineOption(
            "segments",
            "Segments",
            default_enabled=False,
            requires=("per_beat",),
            dag_requires=("prepared_topology",),
        ),
    ],
)
def run(ctx):
    if ctx.option_enabled("segments"):
        ...
```

Use `visibility="hidden"` for shared implementation pipelines that should run
only as dependencies. Test dependency closure and independent option families
in `test_pipeline_library_dependencies.py` or a focused pipeline test.

## Read data safely

Call `ctx.require_inputs(...)` and prefer typed sources:

```python
hd = ctx.inputs.hd.as_holodoppler()
dv = ctx.inputs.dv.as_dopplerview()

moment0 = hd.moment0_dataset()             # lazy (frame, y, x)
artery = dv.retinal_artery_mask()          # eager bool (y, x)
timing = hd.timing()
pixel_pitch = hd.pixel_pitch()
```

For retinal analysis that needs aligned masks, optic-disc geometry, and the
active velocity source, reuse `pipelines.vessel_inputs.load_retinal_source_data`
or an existing pipeline source loader. Do not transpose one mask or hard-code
fallback paths downstream. Exact input contracts are in
[data contracts](docs/data-contracts.md).

When `frequency_bands` is active, `HF / LF` is converted to frequency using the
recorded `band_ratio_frequency_scale_hz`, and velocity remains physical in
`mm/s`. Use `waveform_velocity_core.velocity_semantics` when labeling, plotting,
or packing method-dependent velocity. Do not gate otherwise valid pipelines on
the method. Do not add a dimensionless velocity fallback; incompatible legacy
velocity quantity or unit metadata must fail clearly.

## Publish results

Return a mapping for normal result datasets:

```python
from pipeline_engine import DatasetValue, ProcessResult

metrics = {
    path: DatasetValue(
        data=values,
        attrs={"unit": "1", "dimDesc": ["time", "beat"]},
    )
}
return ProcessResult(metrics=metrics, attrs={"my_metric_version": "1"})
```

Alternatively, write directly with `ctx.output.h5.write(path, value, **attrs)`
or `write_many`; direct-writing pipelines return `None`. Use
`EyeFlowOutputPaths.active()` for existing path families. Before adding a new
literal path, establish one owner and search for downstream consumers.

`ctx.state` is for run-local intermediates, not the only copy of a durable
result. State does not survive files or process restarts.

Artifact locations come from `ctx.output.path_for(...)` or `dir_for(...)` with
an `OutputType`. Supported types are H5, PNG, MP4, AVI, PDF, and EPS. Create
parents before giving a reserved path to a third-party library. The runtime
already owns the primary HDF5; use `ctx.output.open_h5` only for a genuinely
separate auxiliary HDF5 artifact.

## Compatibility and validation checklist

Before finishing a pipeline change:

- confirm the declaration imports during catalog discovery;
- verify option and DAG closure, including hidden dependencies;
- preserve axes, index base, units, sign, NaN, and boundary behavior;
- validate required external data early with clear path/shape errors;
- check both velocity methods if the product consumes velocity;
- update producers, direct consumers, schema documentation, and tests together
  when changing an HDF5 path;
- ensure required results are in the primary HDF5, not only plots or state;
- run the focused tests and then the full suite for cross-cutting changes.

Common mistakes are confusing `requires` with `dag_requires`, registering a
helper as a top-level pipeline, reading whole image volumes when bounded slicing
is available, hard-coding `mm/s`, relying on target order instead of the DAG,
and silently converting missing scientific data to zeros or NaNs.
