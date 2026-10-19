# Test navigation

The suite mixes pytest functions and `unittest.TestCase`; run it through pytest.
Most tests construct small synthetic arrays/HDF5 files and patch expensive
boundaries. They are contracts for scientific ordering, axes, units, NaN and
boundary behavior—not only regression smoke tests.

## Production-to-test map

| Production area | Primary tests | Coverage character |
|---|---|---|
| Pipeline metadata, discovery, DAG, UI selection | `test_pipeline_discovery.py`, `test_pipeline_library_dependencies.py`, `test_blood_volume_rate_pipeline.py`, `test_velocity_analysis_pipeline.py` | unit + integration of declared dependency closure |
| Run specs, input expansion, destinations, execution | `test_run_service.py`, `test_pipeline_failure_isolation.py`, `test_stem_list_inputs.py`, `test_pipeline_context.py` | orchestration and filesystem/HDF5 integration |
| GUI and text-list input selection | `test_input_selection.py`, `test_windows_folder_dialog.py`, `test_stem_list_inputs.py` | quoted cross-drive lists, independent companion/output roots, mixed-file/folder drops, native Windows multi-folder selection, COM cleanup, and shared last-folder persistence across selection changes and restarts |
| Marker-free run-folder inputs | `test_run_folder_inputs.py` | self/nested discovery, matching HD/DV names, optional EF/AE, GUI/CLI execution, marker deduplication, HDF5 validation, output routing, and real-source diagnostics |
| Settings/default/import | `test_default_settings.py`, `test_settings_import.py` | persistence normalization and controller behavior |
| Run logging | `test_logger.py` | level normalization, callbacks, and settings-adjacent snapshot paths |
| Launch dispatch | `test_launcher.py` | CLI-versus-GUI routing |
| HD/DV schema and velocity methods | `test_velocity_pipeline.py`, `test_frequency_band_velocity.py`, `test_frequency_band_pipeline_integration.py`, `test_dual_velocity_workflows.py`, `test_dopplerview_compat.py` | strict data contract, numerical behavior, physical-output integration |
| HD spectograms input | `test_spectograms_input.py` | exact root path, lazy 4-D numeric source, optional workflow input, and preservation of raw exported axes |
| Gap-spectrum preprocessing | `test_gap_spectrogram_preprocessing.py`, `test_topology_capillary_pipeline.py`, `test_dual_velocity_workflows.py` | explicit Y flip and spatial axis swap, spatial median band power, invalid/empty support, nearest-neighbor mask interpolation to HD spectral dimensions, and once-per-run shared HDF5 output |
| Beat-resolved gap spectra | `test_gap_spectrogram_per_beat.py`, `test_dual_velocity_workflows.py` | moments-preferred boundaries with band-ratio fallback, linear phase alignment, NaN gaps, early/late phase medians of MedianBeat (including shared one-third boundary and per-frequency validity), optional skips, dependency scheduling, and shared Spectrograms hierarchy |
| Gap-spectrum and phase metrics | `test_gap_spectrum_metrics.py`, `test_gap_phase_metrics.py`, `test_gap_spectrogram_per_beat.py`, `test_dual_velocity_workflows.py` | pooled and beat-only medians, fixed-peak areas in kHz after the beat median, trapezoidal A1 and common peak-index-distance widths, physical scaling/units and calibration skips, R/Z ratios, closed early/late windows and signed deltas, shared Metrics schema and removed heights, domain-relative peaks/valleys including odd/small domains, quarter-domain neighborhood and overrides, NaN/zero policies, diagnostics, axes, and provenance |
| Low-frequency spectral pipeline | `test_low_freq_spectra_pipeline.py`, `test_low_freq_spectra_artifacts.py`, `test_gap_spectrum_std.py`, `test_dual_velocity_workflows.py` | hidden scheduling and shared execution, timing preference/fallback, per-beat phase medians followed by population STD, valid beat counts, HD sampling-frequency publication, signed full-domain kHz axes, actual positive-half coordinates and shared stored-peak crop, per-curve normalization by the same central-quarter Hv used by A0, styled envelopes, Z1/Z2 phase plots, PNG/EPS routing, read-only figure export, missing calibration/source handling, and failure isolation |
| Spatial and filtered spectra | `test_velocity_spectrograms.py`, `test_gap_spectrogram_signal_processing.py`, `test_dual_velocity_workflows.py` | single-pass full-frame mean/median and gap reduction, shared ProcessingSpectrograms hierarchy, physical FFT cutoff/DC, chronological centered SVD mode removal, NaN support, and moments-failure isolation |
| HDF5 and artifact writers | `test_h5_writer.py`, `test_avi_writer.py`, `test_eps_writer.py`, `test_artifact_names.py` | serialization, attributes, type handling, names, collision-safe migration |
| Cardiac-cycle/per-beat/math | `test_cardiac_cycle_analysis.py`, `test_per_beat_runner.py`, `test_calculations_math.py`, `test_periodic_sliding_windows.py` | numerical unit/reference and boundary behavior |
| Native topology/optic disc and vessel-segment analysis | all `test_topology_*.py`, `test_optic_disc_*.py`, `test_segment_profile_models.py`, `test_calculation_package_boundaries.py` | native geometry; sampling cache identity, transforms, chunks; profile measurement and memory planning; import boundaries |
| Capillary analysis zones | `test_topology_capillary.py`, `test_topology_capillary_pipeline.py`, `test_gap_masks.py`, `test_optic_disc_geometry.py`, `test_dopplerview_compat.py` | legacy cutoff/exclusions, tiled temporal median, moment0ff enhancement/Otsu/smoothing, bright-vessel square dilation before gap inversion, vessel OR before inversion (gap intersection) with final outer/optic-disc constraints and fallback, spectral consumer, DV alignment, topology ownership, missing-input warnings, physical disc fallback, and HDF5 coordinates |
| Velocity analysis outputs | `test_velocity_signal_outputs.py`, `test_cross_section_velocity_profiles.py`, `test_velocity_fft_profiles.py`, `test_segment_map_outputs.py`, `test_segment_map_video.py`, `test_waveform_scratch_analysis.py` | schema, axes, units, artifacts, scratch behavior |
| Spatial-gradient/lumen | `test_spatial_gradient_pipeline.py`, `test_spatial_gradient_profiles.py`, `test_lumen_size_pngs.py`, `test_gaussian2d_blur.py`, `test_unsharpen.py`, relevant topology chunk tests | exact processing order and ImageJ-style references |
| Blood-volume rate | `test_blood_volume_rate.py`, `test_blood_volume_rate_pipeline.py`, `test_topology_mask_area.py` | formula/reference, schema/artifacts, option dependencies |
| Waveform metrics | `test_waveform_shape_metrics_*.py`, `test_absolute_waveform_metrics.py`, `test_lowrank_quadrant_outputs.py` | metric formulas, segment geometry, quadrants, plots |
| Velocity-profile fitting | `test_velocity_profile_analysis.py` | weighted solver, degeneracy/NaN behavior, waveform option/output paths |
| Displacement | `test_displacement_map_pipeline.py` | registration, source validation, temporary field state, and MP4 artifact routing |
| PDF report | `test_pdf_report_image_lookup.py`, `test_pdf_report_runner_paths.py` | artifact lookup and runner paths |
| Installer/release support | `test_installer_script.py` | installer-script contract; release workflow remains separate |

## Selection guidance

Use the map to locate the smallest relevant test cases, not as a mandatory
list of files to run. Follow the proportional-verification policy in the root
`AGENTS.md`: simple output-path/name changes need only a focused contract test,
and documentation-only changes need no test run. Do not automatically run slow
end-to-end workflows or the full suite because a shared schema is touched.
Add producer/consumer integration coverage when a concrete behavior risk
requires it:

- an output path: a focused producer/writer contract test and a reference
  search for consumers; add consumer tests only if consumer behavior changes;
- an array axis: focused producer and affected consumer tests;
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
