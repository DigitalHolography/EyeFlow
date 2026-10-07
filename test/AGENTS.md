# Test navigation

The suite mixes pytest functions and `unittest.TestCase`; run it through pytest.
Most tests construct small synthetic arrays/HDF5 files and patch expensive
boundaries. They are contracts for scientific ordering, axes, units, NaN and
boundary behavior—not only regression smoke tests.

## Production-to-test map

| Production area | Primary tests | Coverage character |
|---|---|---|
| Pipeline metadata, discovery, DAG, UI selection | `test_pipeline_discovery.py`, `test_pipeline_library_dependencies.py`, `test_blood_volume_rate_pipeline.py`, `test_velocity_analysis_pipeline.py` | unit + integration of declared dependency closure |
| Run specs, input expansion, destinations, execution | `test_run_service.py`, `test_stem_list_inputs.py`, `test_pipeline_context.py` | orchestration and filesystem/HDF5 integration |
| Settings/default/import | `test_default_settings.py`, `test_settings_import.py` | persistence normalization and controller behavior |
| Run logging | `test_logger.py` | level normalization, callbacks, and settings-adjacent snapshot paths |
| Launch dispatch | `test_launcher.py` | CLI-versus-GUI routing |
| HD/DV schema and velocity methods | `test_velocity_pipeline.py`, `test_frequency_band_velocity.py`, `test_frequency_band_pipeline_integration.py`, `test_dual_velocity_workflows.py`, `test_dopplerview_compat.py` | strict data contract, numerical behavior, physical-output integration |
| HDF5 and artifact writers | `test_h5_writer.py`, `test_avi_writer.py`, `test_eps_writer.py`, `test_artifact_names.py` | serialization, attributes, type handling, names, collision-safe migration |
| Heartbeat/per-beat/math | `test_heartbeat_analysis.py`, `test_per_beat_runner.py`, `test_calculations_math.py`, `test_periodic_sliding_windows.py` | numerical unit/reference and boundary behavior |
| Topology/optic disc | all `test_topology_*.py`, `test_optic_disc_*.py`, `test_segment_profile_models.py` | geometry, cache identity, transforms, chunks, profiles, memory planning |
| Velocity analysis outputs | `test_velocity_signal_outputs.py`, `test_cross_section_velocity_profiles.py`, `test_velocity_fft_profiles.py`, `test_segment_map_outputs.py`, `test_segment_velocity_map_avi.py`, `test_waveform_scratch_analysis.py` | schema, axes, units, artifacts, scratch behavior |
| Spatial-gradient/lumen | `test_spatial_gradient_pipeline.py`, `test_spatial_gradient_profiles.py`, `test_lumen_size_pngs.py`, `test_gaussian2d_blur.py`, `test_unsharpen.py`, relevant topology chunk tests | exact processing order and ImageJ-style references |
| Blood-volume rate | `test_blood_volume_rate.py`, `test_blood_volume_rate_pipeline.py`, `test_topology_mask_area.py` | formula/reference, schema/artifacts, option dependencies |
| Waveform metrics | `test_waveform_shape_metrics_*.py`, `test_absolute_waveform_metrics.py`, `test_lowrank_quadrant_outputs.py` | metric formulas, segment geometry, quadrants, plots |
| Velocity-profile fitting | `test_velocity_profile_analysis.py` | weighted solver, degeneracy/NaN behavior, waveform option/output paths |
| Displacement | `test_displacement_map_pipeline.py`, `test_displacement_cross_section_outputs.py` | registration, source validation, output schema |
| PDF report | `test_pdf_report_image_lookup.py`, `test_pdf_report_runner_paths.py` | artifact lookup and runner paths |
| Installer/release support | `test_installer_script.py` | installer-script contract; release workflow remains separate |

## Selection guidance

Start with the smallest listed file. Add producer/consumer integration tests
when changing a shared contract:

- an output path or axis: writer/schema test plus producer and every direct
  consumer test;
- topology orientation/interpolation: topology unit tests plus the affected
  velocity, gradient, or displacement output test;
- pipeline dependency/option: DAG/UI dependency tests plus the pipeline runner
  test;
- velocity semantics: band/moment estimator tests plus every output family that
  displays or stores velocity units;
- settings: default, import, CLI/run-spec consumption, and GUI controller only
  where exposed.

Avoid updating expected values until the production behavior is established
from source and the scientific intent is explicit. Existing numerical fixtures
often encode boundary and NaN policy.

## Commands

```powershell
python -m pytest test/test_frequency_band_velocity.py
python -m pytest test/test_run_service.py -k resolve
python -m pytest
```

If the package is not installed editable, set `PYTHONPATH` to `src`. Some
artifact tests skip when their optional rendering dependency is unavailable.
`benchmarks/rtx4090_cross_section.json` is a hardware-specific observation with
parity errors and timings; it is useful when optimizing cross-section code but
is not executed by pytest.
