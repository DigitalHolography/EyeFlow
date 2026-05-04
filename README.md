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

# Contributor setup: base app + pipeline deps + dev tools
uv sync --extra pipelines --extra dev
```

Run commands through uv so you do not need to activate the virtual environment:

```sh
uv run eyeflow
uv run eyeflow-cli --help
uv run lint-tool
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

# Contributor setup
pip install -e ".[dev,pipelines]"
```

### 3. Development Setup

```sh
# After installing dev dependencies with uv
uv run pre-commit install

# If using an activated pip environment instead
pre-commit install
```

After changing dependencies in `pyproject.toml`, run `uv lock` and commit the updated `uv.lock`.

> [!NOTE]
> The pre-commit is really usefull to run automatic checks before pushing code, reducing chances of ugly code being pushed.
>
> If a pre-commit hook fails, it will try to fix all needed files, **so you will need to add them again before recreating the commit**.

> [!TIP]
> You can run the linter easily, once the `dev` dependencies are installed, with the command:
>
> ```sh
> # To only run the checks
> uv run lint-tool
>
> # To let the linter try to fix as much as possible
> uv run lint-tool --fix
> ```

---

## Usage

Launch the main application to process files interactively:

### GUI

The GUI processes one `.holo` file at a time and lets you run multiple coded pipelines on that input. Outputs are written to the default `*_EF` folder derived from the selected `.holo` file.

Use the Pipeline Library tab to select which coded pipelines run. Selection preferences are saved per user between app launches, including installed builds.

```sh
# Via uv
uv run eyeflow

# Or via the script
uv run python src/eye_flow.py
```

When developing from the repository checkout, install the project in editable mode or run `python src/eye_flow.py`; edit pipeline code under `src/pipelines/` and restart the app to pick up code changes.

Installed builds only run the pipelines explicitly registered in source code and bundled into the application. They do not load editable pipeline files from folders next to `EyeFlow.exe` or from environment-variable overrides; adding or changing a pipeline requires rebuilding and reinstalling the application.

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

Pipelines are the heart of EyeFlow. To add a new analysis, create a file in `src/pipelines/` with a class inheriting from `ProcessPipeline`.

To register it to the app, add the decorator `@registerPipeline`, then add an explicit import block for the class in `src/pipelines/__init__.py` so it is appended to the coded pipeline catalog. The app does not scan or import pipeline files dynamically at runtime.

To see complete examples, check out `src/pipelines/waveform_shape_metrics.py` and `src/pipelines/dual_input_tutorial.py`.

### Simple Pipeline Structure

```python
from pipelines import ProcessPipeline, ProcessResult, registerPipeline

@registerPipeline(
    name="My Analysis",
    description="Calculates a custom clinical metric.",
    required_deps=["torch>=2.2"],
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
