# EyeFlow

EyeFlow is the cohort-analysis engine for retinal Doppler holography. It browses EyeFlow .h5 outputs, reads per-segment metrics, applies QC, compares models, and aggregates results at eye/cohort level (including artery–vein summaries) to help design biomarkers. It exports clean CSV reports for stats, figures, and clinical models.

---

## Setup

### Prerequisites

- `uv`.
- Python 3.10 or higher. The repository includes `.python-version` set to Python 3.11 for uv-managed environments.

This project uses `pyproject.toml` plus a committed `uv.lock` for reproducible installs.

### 1. Recommended Setup With uv

```sh
# Install the pinned Python if uv does not already find a compatible one
uv python install

# Create/update .venv and install the base application
uv sync

# Optional: include pipeline-specific dependencies
uv sync --extra pipelines
```

Run commands through uv so you do not need to activate the virtual environment:

```sh
uv run eyeflow
uv run eyeflow-cli --help
```

If you prefer to activate the environment manually:

```sh
# Windows PowerShell
./.venv/Scripts/activate

# macOS/Linux
source .venv/bin/activate
```

### 2. Pip Fallback

```sh
# Base application
pip install -e .

# Optional pipeline dependencies
pip install -e ".[pipelines]"
```

After changing dependencies in `pyproject.toml`, run `uv lock` and commit the updated `uv.lock`.

---

## Usage

Launch the main application to process files interactively:

### GUI

The GUI processes one `.holo` file at a time and lets you select pipeline targets for that input. The DAG resolves the required upstream pipelines and run order. Outputs are written to the default `*_EF` folder derived from the selected `.holo` file.

Use the Pipeline Library tab to select which coded pipeline outputs you want. Selection preferences are saved per user between app launches.

```sh
# Via uv
uv run eyeflow

# Or via the script
uv run python src/eye_flow.py
```

When developing from the repository checkout, install the project in editable mode or run `python src/eye_flow.py`; edit pipeline code under `src/pipelines/` and restart the app to pick up code changes.

### CLI

The CLI is designed for headless processing in environments or clusters.

```sh
# Via uv
uv run eyeflow-cli

# Or via the script
uv run python src/cli.py
```

---

## Pipeline System

Pipelines are the heart of EyeFlow. Pipeline implementations live in `src/pipelines/`; shared execution contracts and DAG code live in `src/pipeline_engine/`. To add a new analysis, create a file in `src/pipelines/` with a class inheriting from `ProcessPipeline`.

To register it to the app, add the decorator `@registerPipeline`, then add an explicit import block for the class in `src/pipelines/__init__.py` so it is appended to the coded pipeline catalog. The app does not scan or import pipeline files dynamically at runtime.

Pipeline execution is DAG-based. The selected pipeline names are treated as targets, every pipeline implicitly produces its own name, and optional `dag_produces` keys can describe intermediate outputs. Use `dag_requires` to declare which produced keys must be available before a pipeline runs. Keys with no producer are treated as external input data.

To see complete examples, check out `src/pipelines/waveform_shape_metrics.py` and `src/pipelines/dual_input_tutorial.py`.

### Simple Pipeline Structure

```python
from pipeline_engine import ProcessPipeline, ProcessResult, registerPipeline

@registerPipeline(
    name="my_analysis",
    description="Calculates a custom clinical metric.",
    required_deps=["torch>=2.2"],
    dag_requires=["velocity_per_beat"],
    dag_produces=["custom_clinical_metric"],
)
class MyAnalysis(ProcessPipeline):
    def run(self, h5file):
        import torch
        # 1. Read data using h5py
        # 2. Perform calculations
        # 3. Return metrics

        metrics={"peak_flow": 12.5}

        # Optional attributes applied to the pipeline group.
        attrs = {
            "pipeline_version": "1.0",
            "author": "StaticExample"
        }

        return ProcessResult(
            metrics=metrics,
            attrs=attrs
        )
```
