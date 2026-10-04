# EyeFlow repository context

EyeFlow is a desktop and command-line analysis application for retinal Doppler
holography. It consumes a `.holo` run plus HoloDoppler (HD) and DopplerView
(DV) companion HDF5 files, resolves a dependency-aware set of analysis
pipelines, and writes one structured EyeFlow result directory per run. It does
not perform vessel segmentation.

## Start here

Read only the context needed for the task:

| Task | Read first | Then inspect |
|---|---|---|
| Entry points, end-to-end flow, settings, or subsystem ownership | [Architecture](docs/architecture.md) | `src/launcher.py`, `src/cli.py`, `src/eye_flow.py`, or `src/app_settings.py` as routed there |
| Pipeline registration, scheduling, options, or execution | [Pipeline engine context](src/pipeline_engine/AGENTS.md) | `src/pipelines/__init__.py` and the affected pipeline declaration |
| A scientific product or pipeline output | [Pipeline context](src/pipelines/AGENTS.md) | The pipeline package, then only the calculation modules it calls |
| A numerical primitive, topology, velocity, or waveform algorithm | [Calculation context](src/calculations/AGENTS.md) | The relevant calculation package and mapped tests |
| Input layout, HDF5 schema, writers, or output folders | [Data contracts](docs/data-contracts.md), then [I/O context](src/input_output/AGENTS.md) | The relevant schema or writer module |
| GUI workflow or presentation | [UI context](src/ui/AGENTS.md) | One controller plus its view; do not load scientific implementations unless the backend contract changes |
| Test failure or regression coverage | [Test map](test/AGENTS.md) | The mapped production module and focused test file |
| Velocity-profile quadratic fitting | [Velocity-profile analysis](docs/velocity_profile_analysis.md) | `src/pipelines/velocity_profile_analysis/` |

## Repository map

| Path | Responsibility |
|---|---|
| `src/pipeline_engine/` | Pipeline metadata, DAG resolution, run specifications, contexts, and execution |
| `src/pipelines/` | Pipeline declarations and orchestration; hidden core pipelines publish shared state for visible products |
| `src/calculations/` | Reusable scientific and numerical code, independent of UI and output-folder policy |
| `src/input_output/` | `.holo` companion discovery, typed HDF5 adapters, output schema, writers, reports, and archives |
| `src/ui/` | Tkinter views and controllers; delegates runs to `pipeline_engine.run_service` |
| `src/utils/` | Process-wide logging; keep domain and orchestration behavior out of this small shared layer |
| `src/app_settings.py` / `default_settings.json` | Persisted settings behavior / fresh-install defaults |
| `test/` | Unit, numerical-reference, pipeline-integration, schema, UI-controller, and launcher tests |
| `benchmarks/` | Checked-in performance/parity observations, not an automated test suite |
| `.github/workflows/release.yml` | Tag-triggered Windows release and installer publication; it is not the test workflow |

## Execution in one view

```text
GUI (`eye_flow.main`) ----+
                         +-> resolve_run_spec -> PipelineDAG -> execute_run
CLI (`cli.main`) --------+                         |
                                                   v
                HD + DV typed sources -> ordered pipelines
                    -> shared run state + EyeFlow HDF5/sidecars
```

Both front ends must converge on `pipeline_engine.run_service`; do not create a
second scheduling or execution path in a controller or CLI handler.

## Global invariants

- Array axes are explicit. HD moments and frequency bands are `(frame, y, x)`;
  masks are `(y, x)`. Do not silently transpose scientific arrays outside the
  canonical input-alignment layer.
- HDF5 paths have owners. External paths belong in `src/input_output/schema/`;
  current EyeFlow output paths come from `EyeFlowOutputPaths.active()`.
- `@pipeline`/`@registerPipeline` declarations are the source of pipeline
  dependencies and options. Discovery scans pipeline modules automatically.
- `dag_requires` names produced data, while `requires` names importable Python
  packages. Do not interchange them.
- `ctx.state` is run-local and shared between ordered pipelines; persisted
  results go through `ctx.output.h5`. Never make a later pipeline depend on an
  undeclared state producer.
- `execute_run` replaces an existing `<stem>_EF` directory before each run.
  Treat that directory as entirely EyeFlow-owned.
- `velocity_estimation_method="frequency_bands"` produces a dimensionless
  relative index, not calibrated velocity. It does not disable downstream
  pipelines. Preserve method/quantity/unit provenance throughout outputs and
  presentation.
- Preserve NaN, boundary, axis, unit, and sign behavior in scientific changes;
  these are tested contracts, not incidental implementation details.
- Keep calculation code independent of Tkinter, settings persistence, and
  filesystem layout. Pipelines adapt calculations to runtime state and output.

## Development commands

From the repository root on PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .
pip install pytest
python -m pytest
```

Run focused tests first, for example:

```powershell
python -m pytest test/test_run_service.py
python -m pytest test/test_frequency_band_velocity.py
```

The source layout requires `src` on `sys.path`; an editable install provides
that. If using an existing interpreter without an editable install, set
`PYTHONPATH` to the repository's `src` directory for the command.

Before changing a contract, search both its producer and consumers. Before
finishing, run the focused tests from [the test map](test/AGENTS.md), then the
full suite when the change crosses subsystem boundaries.
