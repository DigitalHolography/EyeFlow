# EyeFlow repository guide

EyeFlow analyzes retinal Doppler holography from a `.holo` run plus HD and DV
HDF5 companions. It does not segment vessels. GUI and CLI both use
`pipeline_engine.run_service`; one run writes one `<stem>_EF` result directory.

## Read only what the task needs

| Task | Start with |
|---|---|
| Entry points, settings, or ownership | [Architecture](docs/architecture.md) |
| Scheduling, options, or execution | [Pipeline engine map](src/pipeline_engine/AGENTS.md) |
| Scientific product or pipeline output | [Pipeline map](src/pipelines/AGENTS.md) |
| Numerical algorithm or geometry | [Calculation map](src/calculations/AGENTS.md) |
| Input layout, HDF5 schema, or writer | [I/O map](src/input_output/AGENTS.md); use [data contracts](docs/data-contracts.md) when changing a public contract |
| GUI workflow | [UI map](src/ui/AGENTS.md) |
| Test or regression | [Test guide](test/AGENTS.md) |
| Quadratic velocity-profile fit | [Profile analysis](docs/velocity_profile_analysis.md) |

Follow the relevant map to the owning module and focused tests. Search both
producers and consumers before changing a contract.

## Repository boundaries

- `src/pipeline_engine/` owns metadata, DAG resolution, run specs, and execution.
- `src/pipelines/` declares dependencies/options, adapts sources and
  calculations, owns run state, and publishes products.
- `src/calculations/` owns reusable science and has no UI, settings, or
  output-directory dependency.
- `src/input_output/` owns HD/DV adapters, schemas, output paths, and writers.
- `src/ui/` owns Tkinter presentation; `src/app_settings.py` owns persisted
  settings. Both front ends must converge on the shared run service.

## Contracts to preserve

- HD image volumes use `(frame, y, x)`; native masks use `(y, x)`. Published
  2-D EyeFlow maps use lower-left `(x, y)` alignment. Do not transpose outside
  the canonical input/output boundaries.
- External HDF5 paths belong in `src/input_output/schema/`;
  `EyeFlowOutputPaths.active()` owns current EyeFlow paths.
- Pipeline declarations own options and dependencies. `requires` names Python
  packages; `dag_requires` names produced data. A `ctx.state` consumer needs
  a declared producer. Durable results go through `ctx.output.h5`.
- `execute_run` replaces the EyeFlow-owned `<stem>_EF` directory for each
  attempt. Validate destinations before execution.
- Selected velocity workflows produce physical `mm/s`. Band-ratio estimation
  converts `HF/LF` to Hz with explicit calibration before the shared velocity
  conversion. Preserve method, calibration, quantity, unit, and timing
  provenance; reject dimensionless legacy velocity metadata.
- Scientific changes must preserve defined NaN, boundary, axis, unit, and sign
  behavior unless the contract itself is intentionally changed.

## Documentation and instruction maintenance

Update `docs/` only for a material change to user workflow, architecture,
scientific meaning, or an input/output contract, or to correct guidance that is
wrong. Do not revise docs for routine refactors or add implementation diaries.
When text is obsolete, duplicated, or no longer useful, remove it rather than
adding another qualification. Keep detailed numerical facts in the owning
code/tests and focused reference pages.

Update any `AGENTS.md` sparingly: only when a durable ownership boundary,
contract, or navigation path changes. Keep it a short map for future work;
do not record every implementation step, incident, or test result. Delete stale
instructions and avoid repeating rules already in this root file.

## Verification

For documentation-only changes, do not run tests. For code changes, start
with the smallest focused test that protects the affected contract; expand only
for a concrete uncovered risk or when requested. See [test guidance](test/AGENTS.md).
The source layout needs an editable install or `PYTHONPATH=src`.
