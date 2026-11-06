# Test navigation

Run pytest from the repository root with the editable-install interpreter:
`python -m pytest test/test_run_service.py`. Without an editable install, put
`src` on `PYTHONPATH`. The suite mixes pytest and unittest; optional rendering
or registration dependencies may cause skips.

## Production-to-test map

| Area | Focused tests |
|---|---|
| Catalog, dependencies and options | `test_pipeline_discovery.py`, `test_pipeline_library_dependencies.py`, affected `test_*_pipeline.py` |
| Run specs, dispatch and failure isolation | `test_run_service.py`, `test_pipeline_failure_isolation.py`, `test_pipeline_context.py` |
| Local stage planning/execution | `test_stage_graph_ownership.py`, `test_pipeline_graphs.py` |
| Product exchange and worker lifecycle | `test_product_exchange.py`, `test_product_exchange_processes.py`, `product_exchange_fixtures.py` |
| Replay, output and admission | `test_input_replay.py`, `test_output_model.py`, `test_process_admission.py`, `test_parallel_coordinator.py`, `test_resource_budget.py` |
| Input selection and run layouts | `test_input_selection.py`, `test_stem_list_inputs.py`, `test_run_folder_inputs.py`, `test_windows_folder_dialog.py` |
| Settings, logging and launch | `test_default_settings.py`, `test_settings_import.py`, `test_logger.py`, `test_launcher.py` |
| HD/DV inputs and physical velocity | `test_velocity_pipeline.py`, `test_frequency_band_velocity.py`, `test_frequency_band_pipeline_integration.py`, `test_dual_velocity_workflows.py`, `test_dopplerview_compat.py`, `test_spectograms_input.py` |
| Spectral pipeline and artifacts | `test_low_freq_spectra_pipeline.py`, `test_low_freq_spectra_artifacts.py`, `test_spectrogram_frequency.py`, `test_velocity_spectrograms.py` |
| Spectral calculations | `test_gap_spectrogram_*.py`, `test_gap_spectrum_*.py`, `test_gap_phase_metrics.py`, `test_gap_flood_metrics.py`, `test_phase_windows.py` |
| Math and cardiac timing | `test_calculations_math.py`, `test_cardiac_cycle_analysis.py`, `test_per_beat_runner.py`, `test_periodic_sliding_windows.py` |
| Topology, masks and sampling | `test_topology_*.py`, `test_optic_disc_*.py`, `test_gap_masks.py`, `test_segment_profile_models.py`, `test_calculation_package_boundaries.py` |
| Velocity outputs | `test_velocity_signal_outputs.py`, `test_cross_section_velocity_profiles.py`, `test_velocity_fft_profiles.py`, `test_velocity_profile_analysis.py`, `test_segment_map_outputs.py`, `test_segment_map_video.py`, `test_waveform_scratch_analysis.py` |
| Spatial gradients and lumen geometry | `test_spatial_gradient_pipeline.py`, `test_spatial_gradient_profiles.py`, `test_lumen_size_pngs.py`, `test_gaussian2d_blur.py`, `test_unsharpen.py` |
| Waveform metrics and volume rate | `test_waveform_shape_metrics_*.py`, `test_absolute_waveform_metrics.py`, `test_lowrank_quadrant_outputs.py`, `test_blood_volume_rate.py`, `test_blood_volume_rate_pipeline.py` |
| Displacement | `test_displacement_map_pipeline.py`, `test_displacement_profile_outputs.py` |
| Writers, reports and artifact names | `test_h5_writer.py`, `test_avi_writer.py`, `test_artifact_names.py`, `test_pdf_report_runner_paths.py` |
| Installer/release support | `test_installer_script.py`, `test_installer_builder.py` |

## Selection and maintenance

Start with the mapped test and its producer/consumer. Use focused checks for
path, naming or presentation changes; add integration coverage for actual
producer/consumer risks. Run broader tests only when requested or a concrete
risk cannot be checked narrowly. Do not repeat passing checks without a relevant
change or unresolved concern. Documentation-only edits need no tests.

Tests verify established contracts; they do not define intended behavior. Keep
scientific axes, units, NaNs, sign and boundary references intact. Prefer updating
an existing regression to adding another fixture/file, and remove obsolete or
redundant cases when changing a contract.
