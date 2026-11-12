# Test guide

Use the smallest test that uniquely protects a scientific, schema, scheduling,
or user-visible contract. Search existing coverage first; extend a focused
case when it already reaches the relevant failure mode. Test the rule at its
owning layer and add integration coverage when the boundary wiring itself can
regress. Avoid full artifact families and tests that merely mirror code.

| Change | Start with |
|---|---|
| Numerical formula, NaN, axes, or units | Owning reference test, such as `test_cardiac_cycle_analysis.py`, `test_velocity_profile_analysis.py`, or `test_spectral_spatial_reduction.py` |
| HDF5 path, source alignment, or serialization | `test_h5_writer.py`, `test_dopplerview_compat.py`, or the producer's output test |
| Pipeline discovery, option, or DAG dependency | `test_pipeline_discovery.py`, `test_pipeline_library_dependencies.py` |
| Run service or velocity workflow selection | `test_run_service.py`, `test_velocity_workflow_selection.py`, `test_frequency_band_velocity.py` |
| Stage graph or process coordination | `test_pipeline_graphs.py`, `test_parallel_coordinator.py` |
| GUI interaction or settings | `test_input_selection.py`, `test_settings_import.py`, or the affected controller test |
| Installer or release logic | `test_bundle_policy.py`, `test_installer_builder.py`, `test_release_workflow.py` |

Documentation or instruction changes need no runtime test. Run the narrowest
pytest node or file first; expand only after a failure or a concrete uncovered
risk. Run the full suite when requested or when focused checks cannot establish
the changed behavior.

Find a specific case with `rg -n 'def test_.*<behavior>' test/`, then run its
pytest node (`python -m pytest test/test_file.py::test_name`). Measure duration
with `--durations=10` when a focused check is slow; classify tests from those
measurements before creating a fast/slow split. If one focused test takes over
a minute, reduce its fixture or isolate the contract before adding similar
coverage.
