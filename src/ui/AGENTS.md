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

Input interactions are covered by `test_input_selection.py` (popup/drop
controllers, mixed directories, native dialog delegation, and output routing)
and `test_stem_list_inputs.py` (quoted lists and HD/DV resolution).

The input row exposes separate native file and folder dialogs. On Windows,
`TkFileDialogService.askdirectories` uses `windows_folder_dialog.py` and the
system `IFileOpenDialog` with `FOS_PICKFOLDERS | FOS_ALLOWMULTISELECT` to support
Shift/Ctrl folder selection. Other platforms retain Tk's single-folder fallback.
Both dialogs append via
`InputController.add_holo_input_paths`; folders use the same recursive expansion
and deduplication as drops. Clear resets the selection. A text list stays a sole
input: selecting one replaces files, while adding files/folders replaces the list.
Directory expansion delegates to `input_output.runs.discovery`, which accepts
matching HD/DV run directories without `.holo` markers, including the selected
directory itself. Keep folder recognition and path identities out of the UI.
The file and folder pickers share `last_input_directory`, persisted through
`AppSettingsStore`. Remember the chosen folder before expanding it into runs;
appending inputs and Clear must not reset the picker's location. Cancel leaves
it unchanged, and an unavailable saved folder falls back to the selected input
directory, then the working directory.
For multiple folders, remember their containing directory so sibling selections
reopen at the same browsing level.
