# GUI context

## Boundary

The Tkinter UI owns presentation and user interaction. It does not own pipeline
dependency logic, run validation, scientific calculations, external HDF5 path
rules, or settings-file serialization. Those services must remain usable by the
CLI.

## Structure

- `app.py` composes application state, controllers, and views.
- `views/minimal.py` and `views/advanced.py` build the two layouts;
  `views/app.py` contains shared view construction.
- `controllers/input.py` owns selected input interaction.
- `controllers/pipeline_library.py` displays descriptors/options and their DAG
  requirements; it delegates actual resolution to `PipelineDAG`.
- `controllers/run.py` builds a `RunSpec`, starts execution off the UI thread,
  and translates callbacks/results into UI state.
- `controllers/settings.py` imports settings and switches view mode through
  `AppSettingsStore`.
- progress, resources, and view controllers isolate their named concerns.
- `services.py` contains UI-facing services; `widgets.py` contains reusable
  widgets.

## Workflow contract

Both velocity estimators run automatically; do not expose a method selector.
The selected target names and option names are passed to
`pipeline_engine.resolve_run_spec`; the controller must not duplicate its
availability, dependency, or destination validation. Execution delegates to
`execute_run`. UI callbacks may update widgets, but scientific work must not run
on the Tk event thread.

Pipeline library visibility is user selection, not resolved execution order.
Required hidden or visible dependencies come from the plan and should be shown
as required rather than persisted as user-selected targets.

Settings persistence goes through `AppSettingsStore`. An imported settings file
is validated before replacing active settings. Do not write JSON directly from
a controller.

## Change scopes

For a layout-only change, inspect one view and the controller fields it binds;
the calculation and I/O packages are irrelevant. For run workflow changes,
inspect `controllers/run.py`, `pipeline_engine/run_service.py`, and
`test_run_service.py`. For pipeline-library behavior, inspect its controller,
descriptor/DAG interfaces, and `test_pipeline_library_dependencies.py`. For
settings import/mode behavior, inspect the settings controller,
`app_settings.py`, `test_settings_import.py`, and `test_default_settings.py`.

Most GUI behavior is tested at controller/helper boundaries without creating a
real window. Keep new logic separable from widget construction so it remains
testable in the same style.
