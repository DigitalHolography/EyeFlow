# UI map

The Tkinter UI owns presentation and interaction. Shared validation,
scheduling, science, HDF5 paths, and settings persistence stay in their owning
packages so CLI and GUI use the same behavior.

| Task | Start with |
|---|---|
| App composition and views | `app.py`, `views/` |
| Input selection and folder/drop handling | `controllers/input.py`; discovery in `src/input_output/runs/` |
| Pipeline and option display | `controllers/pipeline_library.py`; declarations and `PipelineDAG` |
| Run/progress/cancellation | `controllers/run.py`, `controllers/progress.py`; `pipeline_engine.run_service` |
| Settings import and view mode | `controllers/settings.py`; `src/app_settings.py` |
| Reusable widgets and services | `widgets.py`, `services.py` |

Controllers pass selected targets/options to `resolve_run_spec` and execution
to `execute_run`. Do not reimplement dependency or destination checks in a
controller. Run science off the Tk event thread. Settings changes use
`AppSettingsStore`. For a presentation change, inspect one view and its
controller; use focused controller tests without creating a real window when
possible.
