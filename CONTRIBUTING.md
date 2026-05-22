# Contributing

This guide explains how to add a new EyeFlow pipeline.

If you only remember one thing, remember this:

1. Create a Python module in `src/pipelines/`.
2. Decorate one `run(ctx)` function with `@pipeline(...)`.
3. Add the module name to `PIPELINE_MODULES` in `src/pipelines/__init__.py`.
4. Read inputs from `ctx`.
5. Write outputs through `ctx`.

## Repository Setup

From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .
```

If your pipeline needs optional pipeline dependencies:

```powershell
pip install -e .[pipelines]
```

## What A Pipeline Is

A pipeline is one analysis step that the GUI or CLI can run.

Each pipeline receives a `PipelineContext` named `ctx`. The context gives your
pipeline access to:

- the Holodoppler input H5 through `ctx.inputs.hd`
- the DopplerView input H5 through `ctx.inputs.dv`
- the EyeFlow work H5 handle through `ctx.runtime.work_h5`
- the output H5 writer and sidecar manager through `ctx.output`
- shared in-memory values through `ctx.state.set(...)` and `ctx.state.get(...)`
- merged input/config attributes through `ctx.attrs`

The pipeline runtime is in `src/pipeline_engine/`.

Useful examples already in the repo:

- `src/pipelines/test1.py`: tiny HD-only pipeline
- `src/pipelines/test2.py`: tiny pipeline that depends on `test1`
- `src/pipelines/dual_input_tutorial.py`: reads HD and DV inputs together
- `src/pipelines/package_pipeline_tutorial/`: package-style pipeline skeleton
- `src/pipelines/waveform_shape_metrics/`: larger real pipeline

## Important Project Context

EyeFlowPython does not create vessel masks. The masks and moment datasets are
expected to already exist in the input H5 files.

The current coded pipeline runner starts from a `.holo` selection. It resolves
that selection into two companion HDF5 inputs:

- HD: Holodoppler data, for example `moment0` and `moment2`
- DV: DopplerView data, for example artery and vein masks

The README also describes the broader input contract:

- moment datasets such as `/moment0`, `/moment1`, `/moment2`
- binary vessel masks such as `/masks/artery` and `/masks/vein`
- explicit timing and spatial calibration metadata

When writing a coded pipeline today, prefer the typed source adapters:

```python
hd = ctx.inputs.hd.as_holodoppler()
dv = ctx.inputs.dv.as_dopplerview()
moment0 = hd.moment0()
artery_mask = dv.retinal_artery_mask()
```

These logical names are mapped to real H5 paths by the schema files in
`src/input_output/schema/`.

## The Smallest Useful Pipeline

Create a new file:

```text
src/pipelines/my_metric.py
```

Add this code:

```python
from pipeline_engine.imports import np, pipeline


@pipeline(
    name="my_metric",
    description="Compute the mean moment0 value for each frame.",
    input_slot="hd",
    dag_produces=["my_metric"],
)
def run(ctx) -> None:
    ctx.require_inputs("hd")

    moment0 = ctx.inputs.hd.array("moment0", dtype=np.float32)
    mean_per_frame = np.nanmean(moment0, axis=(1, 2)).astype(np.float32)

    ctx.output.h5.write(
        "analysis/my_metric/mean_moment0_per_frame",
        mean_per_frame,
        unit="a.u.",
        dimDesc=["frame"],
    )

    ctx.output.h5.set_attr("my_metric_input_path", "moment0")
```

What this does:

- `@pipeline(...)` registers the function as a runnable pipeline.
- `input_slot="hd"` tells the runtime this pipeline mainly uses HD data.
- `ctx.require_inputs("hd")` fails early with a clear error if HD is missing.
- `ctx.inputs.hd.array("moment0")` reads the HD `moment0` dataset.
- `ctx.output.h5.write(...)` creates a dataset in the output EyeFlow H5 file.
- `ctx.output.h5.set_attr(...)` writes metadata to the output H5 file.

## Register The Pipeline

Open:

```text
src/pipelines/__init__.py
```

Add your module name to `PIPELINE_MODULES`:

```python
PIPELINE_MODULES = (
    "dual_input_tutorial",
    "package_pipeline_tutorial",
    "test1",
    "test2",
    "waveform_shape_metrics",
    "my_metric",
)
```

The catalog is explicit on purpose. Helper files should not become runnable
pipelines by accident.

## Reading Inputs From `ctx`

Use `ctx.require_inputs(...)` at the top of your pipeline.

```python
ctx.require_inputs("hd")
ctx.require_inputs("dv")
ctx.require_inputs("hd", "dv")
```

Read HD datasets:

```python
moment0 = ctx.inputs.hd.array("moment0", dtype=np.float32)
moment2 = ctx.inputs.hd.array("moment2", dtype=np.float32)
```

Read DV datasets:

```python
dv = ctx.inputs.dv.as_dopplerview()
artery_mask = dv.retinal_artery_mask()
vein_mask = dv.retinal_vein_mask()
labels = dv.retinal_labeled_vessels()
```

Read a DV config value. This checks the H5 first, then the sidecar JSON config,
then the schema default:

```python
local_background_dist = ctx.inputs.dv.as_dopplerview().local_background_dist()
```

Read Holodoppler timing:

```python
from pipeline_engine.imports import resolve_holodoppler_timing

timing = resolve_holodoppler_timing(ctx)
dt_seconds = timing.dt_seconds
```

## Writing Outputs

Use the output H5 for pipeline metrics and machine-readable values that later
pipelines or downstream tools should consume.

Use `ctx.output.h5.write(...)` for one output:

```python
ctx.output.h5.write(
    "analysis/my_metric/value",
    np.float32(42.0),
    unit="a.u.",
)
```

Use `ctx.output.h5.write_many(...)` for several outputs:

```python
ctx.output.h5.write_many(
    {
        "analysis/my_metric/artery_pixel_count": (
            np.int32(np.count_nonzero(artery_mask)),
            {"unit": "pixels"},
        ),
        "analysis/my_metric/vein_pixel_count": (
            np.int32(np.count_nonzero(vein_mask)),
            {"unit": "pixels"},
        ),
    }
)
```

Dataset paths are relative to the output H5 root. A path like
`analysis/my_metric/value` creates groups automatically.

Prefer output attributes that explain the value:

- `unit`: physical unit, such as `"s"`, `"mm/s"`, `"pixels"`, or `"a.u."`
- `dimDesc`: dimension names, such as `["frame"]` or `["beat", "sample"]`

## Writing Output Files

Each run also has a structured output folder managed by `ctx.output`. Use it
for sidecar files such as JSON summaries, figures, movies, or auxiliary H5
artifacts.

The default layout for a HOLO file named `sample.holo` is:

```text
sample/
    sample_EF/
        h5/sample_EF.h5
        json/
        png/
        mp4/
        avi/
```

The main H5 file is opened by the runtime. Most pipelines should write to that
file with `ctx.output.h5.write(...)`, not by opening another H5 file.

Write a JSON sidecar:

```python
ctx.output.write_json(
    {
        "pipeline": ctx.runtime.pipeline_name,
        "frame_count": int(moment0.shape[-1]),
        "mean_value": float(np.nanmean(moment0)),
    },
    "my_metric_summary.json",
)
```

Reserve a path for an image or movie that another library writes:

```python
from input_output.output_manager import OutputType


figure_path = ctx.output.path_for(OutputType.PNG, "my_metric_overlay.png")
figure_path.parent.mkdir(parents=True, exist_ok=True)

# Example:
# matplotlib_figure.savefig(figure_path, dpi=150)
```

Use `ctx.output.dir_for(...)` when a pipeline writes several files of the same
kind:

```python
png_dir = ctx.output.dir_for(OutputType.PNG)
png_dir.mkdir(parents=True, exist_ok=True)
```

Use a separate H5 file only for an auxiliary artifact that should not be part of
the main result H5:

```python
with ctx.output.open_h5("my_metric_debug.h5") as debug_h5:
    debug_h5.create_dataset("intermediate", data=intermediate_array)
```

Do not store required pipeline metrics only in sidecar files. If a value is part
of the analysis result, write it to the main output H5 with
`ctx.output.h5.write(...)`.

## Reading Values From The Output H5

Inside a later pipeline, read a value already written to the EyeFlow output H5
with `ctx.output.h5.read(...)` or `ctx.output.h5.array(...)`:

```python
mean_per_frame = ctx.output.h5.array(
    "analysis/my_metric/mean_moment0_per_frame",
    dtype=np.float32,
)
```

Outside the pipeline runtime, use `h5py` and the same dataset path:

```python
import h5py

with h5py.File("path/to/output_eyeflow.h5", "r") as h5:
    mean_per_frame = h5["analysis/my_metric/mean_moment0_per_frame"][()]
    unit = h5["analysis/my_metric/mean_moment0_per_frame"].attrs.get("unit")
```

## Sharing Values Between Pipelines

Use `ctx.state.set(...)` for temporary values needed by later pipelines but not
written to the output H5.

Producer:

```python
from pipeline_engine.imports import np, pipeline


@pipeline(
    name="prepare_velocity_signal",
    description="Prepare one velocity signal for later pipelines.",
    input_slot="both",
    dag_produces=["velocity_signal"],
)
def run(ctx) -> None:
    ctx.require_inputs("hd", "dv")

    signal = np.asarray([1.0, 2.0, 3.0], dtype=np.float32)
    ctx.state.set("velocity_signal", signal)
```

Consumer:

```python
from pipeline_engine.imports import np, pipeline


@pipeline(
    name="summarize_velocity_signal",
    description="Summarize the velocity signal prepared upstream.",
    input_slot="both",
    dag_requires=["velocity_signal"],
)
def run(ctx) -> None:
    signal = ctx.state.get("velocity_signal")
    if signal is None:
        raise RuntimeError("Expected prepare_velocity_signal to set velocity_signal.")

    ctx.output.h5.write(
        "analysis/velocity_signal/mean",
        np.float32(np.nanmean(signal)),
        unit="mm/s",
    )
```

If the user selects `summarize_velocity_signal`, the DAG automatically runs
`prepare_velocity_signal` first because of `dag_requires`.

## DAG Rules

Use DAG keys to express pipeline dependencies.

- `dag_produces=["some_key"]`: this pipeline creates a named result.
- `dag_requires=["some_key"]`: this pipeline needs a named result.
- A pipeline also implicitly produces its own pipeline name.
- If no pipeline produces a required key, the key is treated as external input.

Good DAG keys describe data, not implementation details:

```python
dag_produces=["velocity_per_beat"]
dag_requires=["velocity_per_beat"]
```

Avoid vague keys:

```python
dag_produces=["step1"]
dag_requires=["thing"]
```

## Optional Python Dependencies

If your pipeline needs an import that is not always installed, declare it:

```python
@pipeline(
    name="my_pandas_pipeline",
    description="Example pipeline with an optional dependency.",
    requires=["pandas"],
    input_slot="both",
)
def run(ctx) -> None:
    import pandas as pd

    ...
```

If `pandas` is missing, the pipeline appears as unavailable instead of crashing
the whole catalog import.

Use the importable module name in `requires`. For example, use `"skimage"` for
scikit-image.

If the dependency should be part of normal project installation, also add it to
`pyproject.toml`.

## When To Use A Package Pipeline

A single file is best for small pipelines.

Use a package when the pipeline grows and needs helper modules:

```text
src/pipelines/my_big_pipeline/
    __init__.py
    runner.py
    models.py
```

`__init__.py` registers the pipeline:

```python
from pipeline_engine.imports import ProcessResult, pipeline

from .runner import run_my_big_pipeline


@pipeline(
    name="my_big_pipeline",
    description="Compute a larger metric with helper modules.",
    input_slot="both",
    dag_produces=["my_big_pipeline"],
)
def run(ctx) -> ProcessResult:
    metrics, attrs = run_my_big_pipeline(ctx)
    return ProcessResult(metrics=metrics, attrs=attrs)
```

`runner.py` contains the implementation:

```python
from pipeline_engine.imports import np


def run_my_big_pipeline(ctx) -> tuple[dict[str, object], dict[str, object]]:
    ctx.require_inputs("hd", "dv")

    moment0 = ctx.inputs.hd.array("moment0", dtype=np.float32)
    artery_mask = ctx.inputs.dv.as_dopplerview().retinal_artery_mask()

    artery_mean = np.float32(np.nanmean(moment0[:, artery_mask]))

    metrics = {
        "analysis/my_big_pipeline/artery_mean_moment0": (
            artery_mean,
            {"unit": "a.u."},
        )
    }
    attrs = {
        "my_big_pipeline_hd_path": "moment0",
        "my_big_pipeline_dv_path": "segmentation/Retina/artery_mask",
    }
    return metrics, attrs
```

Register the package exactly like a single-file pipeline:

```python
PIPELINE_MODULES = (
    ...
    "my_big_pipeline",
)
```

## Running Your Pipeline

First check that the catalog imports:

```powershell
python -c "from pipelines import load_pipeline_catalog; a, m = load_pipeline_catalog(); print('available:', [p.name for p in a]); print('missing:', [p.name for p in m])"
```

To run from the pipeline CLI, create a text file with the target pipeline names:

```text
my_metric
```

Then run:

```powershell
eyeflow-cli --data path\to\input.holo --pipelines path\to\pipelines.txt --output path\to\outputs
```

You can also pass a folder containing `.holo` files:

```powershell
eyeflow-cli --data path\to\folder --pipelines path\to\pipelines.txt --output path\to\outputs
```

To run from the GUI:

```powershell
eyeflow
```

Open the Pipeline Library tab, select your pipeline target, choose an input, and
run it.

## Common Mistakes

- Forgetting to add the module to `PIPELINE_MODULES`.
- Reading from `ctx.inputs.hd` or `ctx.inputs.dv` without calling
  `ctx.require_inputs(...)`.
- Writing outputs without units when the value has a physical meaning.
- Using `ctx.state` for values that should be saved in the output H5.
- Using `ctx.output.h5.write(...)` for large temporary values that only the
  next pipeline needs.
- Hiding missing data by returning `nan` when the pipeline should fail clearly.
- Creating a DAG key that no producer creates and then expecting the runtime to
  build it automatically.

## Style Guidelines

- Keep one pipeline responsible for one clear analysis step.
- Put domain calculations in small functions that are easy to test.
- Prefer schema logical names such as `"moment0"` and `"retinal_artery_mask"`
  over hard-coded H5 paths.
- Convert numeric outputs to 32-bit NumPy values where practical.
- Keep output paths stable once users or downstream tools depend on them.
- Raise clear errors when required data is missing or malformed.
- Do not run AI segmentation inside EyeFlowPython. Use masks already present in
  the input data.
