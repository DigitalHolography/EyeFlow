# Waveform Shape Metrics Matlab Trace Review

This document maps the recent Python changes for beat division and velocity-signal smoothing back to the old Matlab EyeFlow implementation.

## 1. Beat Division

### Matlab Source Of Truth

File: `../eyeflow/EyeFlow/BloodFlowVelocity/find_systole_index.m`

Function: `find_systole_index`

Context: This function detects systolic acceleration peak indexes. Those indexes become the beat boundaries used by per-beat velocity calculations.

Relevant lines:

- `find_systole_index.m:18`
  - Defines the default systole low-pass cutoff: `options.lowpass_freq = 15`.
- `find_systole_index.m:23-24`
  - Computes the temporal sampling used by the detector:
    - `fs = ToolBox.fs * 1000 / ToolBox.stride`
    - `dt = 1 / fs`
- `find_systole_index.m:28-30`
  - Applies a 4th-order Butterworth low-pass filter to the artery pulse:
    - `[b, a] = butter(4, options.lowpass_freq / (fs / 2), 'low')`
    - `pulse_artery_filtered = filtfilt(b, a, pulse_artery)`
- `find_systole_index.m:36-37`
  - Computes the derivative used for systole detection:
    - `diff_artery_signal = gradient(pulse_artery_filtered)`
- `find_systole_index.m:43-49`
  - Detects systolic peaks:
    - `min_duration = 0.5`
    - `min_peak_height = prctile(diff_artery_signal, 95)`
    - `min_peak_distance = floor(min_duration / dt)`
    - `findpeaks(..., 'MinPeakHeight', ..., 'MinPeakDistance', ...)`
- `find_systole_index.m:51-52`
  - Removes peaks closer than 10 frames:
    - `sys_idx_list = validate_peaks(sys_idx_list, 10)`
- `find_systole_index.m:150-164`
  - Implementation of `validate_peaks`, which removes the second peak in any too-close pair.

File: `../eyeflow/EyeFlow/BloodFlowVelocity/perBeatSignalAnalysis.m`

Function: `perBeatSignalAnalysis`

Context: This confirms how detected systolic indexes become beat count.

Relevant lines:

- `perBeatSignalAnalysis.m:8`
  - Beat count is `numberOfBeats = numel(sysIdxList) - 1`.
- `perBeatSignalAnalysis.m:9`
  - Per-beat sample count is `N_fft = 2 ^ nextpow2(max(diff(sysIdxList)))`.
- `perBeatSignalAnalysis.m:14-18`
  - Each beat is sliced from one systolic index to the next, then Fourier-interpolated.

File: `../eyeflow/EyeFlow/BloodFlowVelocity/perBeatAnalysis.m`

Function: `perBeatAnalysis`

Context: This shows where the cached systolic indexes are consumed by per-beat analysis.

Relevant lines:

- `perBeatAnalysis.m:6`
  - Reads `BandLimitedSignalHarmonicCount` from `params.json.PulseAnalysis`.
- `perBeatAnalysis.m:7`
  - Reads `sys_idx_list` from `ToolBox.Cache.sysIdxList`.
- `perBeatAnalysis.m:11-12`
  - Writes beat periods from `diff(sys_idx_list)`.
- `perBeatAnalysis.m:14-15`
  - Runs per-beat analysis for vein and artery using the same `sys_idx_list`.

### Python Implementation

File: `src/calculations/blood_flow_velocity/find_systole_index.py`

Function: `find_systole_index`

Context: New Python equivalent of Matlab `BloodFlowVelocity/find_systole_index.m`.

Relevant lines:

- `find_systole_index.py:19-26`
  - Public systole-detection function.
  - Defaults mirror Matlab:
    - `lowpass_freq_hz = 15.0`
    - `min_duration_seconds = 0.5`
    - `validation_distance = 10`
- `find_systole_index.py:29-32`
  - Filters the artery signal, computes derivative, computes 95th-percentile threshold, and computes peak distance.
- `find_systole_index.py:33-38`
  - Calls SciPy `find_peaks` with height and distance, then validates too-close peaks.
- `find_systole_index.py:50-59`
  - Implements Matlab's 4th-order Butterworth + `filtfilt`.
- `find_systole_index.py:71-74`
  - Implements Matlab's `floor(0.5 / dt)` minimum peak distance.
- `find_systole_index.py:77-84`
  - Implements the same keep-first/remove-next close-peak validation pattern.

File: `src/calculations/steps/arterial_waveform_analysis.py`

Class/function: `ArterialWaveformAnalysisStep.run`

Context: This step calls the new Matlab-aligned systole detector and stores the detected beat indexes in the pipeline cache.

Relevant lines:

- `arterial_waveform_analysis.py:52-54`
  - Calls `find_systole_index(sig, dt_seconds=stride / fs)`.
  - Stores:
    - `peaks = detection.systole_indexes`
    - `sig_filtered = detection.artery_signal_filtered`
- `arterial_waveform_analysis.py:58-60`
  - Writes:
    - `retinal_artery_velocity_signal_filtered_perbeat`
    - `retinal_artery_velocity_signal_filtered`
    - `beat_indices`
- `arterial_waveform_analysis.py:61-66`
  - Computes `time_per_beat` from `np.diff(peaks) * stride / fs`.
- `arterial_waveform_analysis.py:67-68`
  - Stores debug values for the peak detector.

Class/function: `ArterialWaveformAnalysisStep.slice_interp_beats`

Context: This was corrected because the beat dimension must be `len(peaks) - 1`, matching Matlab `numel(sysIdxList) - 1`.

Relevant lines:

- `arterial_waveform_analysis.py:28-33`
  - Uses `nbeat = max(0, len(peaks) - 1)`.
- `arterial_waveform_analysis.py:35-42`
  - Fills exactly one row per interval between consecutive systolic indexes.

File: `src/calculations/blood_flow_velocity/per_beat_signal.py`

Function: `_analyze_cycles`

Context: Existing Python per-beat implementation already matched the Matlab beat-count rule.

Relevant lines:

- `per_beat_signal.py:51`
  - Beat count is `cycle_boundaries.size - 1`.
- `per_beat_signal.py:52`
  - Uses max period length to compute power-of-two FFT length.
- `per_beat_signal.py:55-66`
  - Slices each beat boundary interval, interpolates to `n_fft + 1`, drops last sample, FFTs, then band-limits.

## 2. Velocity Signal Smoothing

### Matlab Source Of Truth

File: `../eyeflow/EyeFlow/BloodFlowVelocity/pulseAnalysis.m`

Function: `pulseAnalysis`

Context: Old EyeFlow smooths artery and vein velocity signals before caching them and before systole detection/per-beat analysis.

Relevant lines:

- `pulseAnalysis.m:242-248`
  - Conditional smoothing block:
    - `if filterSignals`
    - `[b, a] = butter(4, 15 / (fs / 2), 'low')`
    - `v_artery_signal = filtfilt(b, a, v_artery_signal)`
    - `v_vein_signal = filtfilt(b, a, v_vein_signal)`
- `pulseAnalysis.m:258-263`
  - Writes and caches the filtered signals:
    - `/Artery/Velocity/VelocitySignal`
    - `/Vein/Velocity/VelocitySignal`
    - `ToolBox.Cache.ArterialVelocity`
    - `ToolBox.Cache.VeinVelocity`
- `pulseAnalysis.m:271-272`
  - Runs systole detection on the filtered velocity signals:
    - `find_systole_index(v_artery_signal, 'pulseVein', v_vein_signal)`
- `pulseAnalysis.m:276`
  - Writes systolic acceleration peak indexes to `/Artery/Velocity/SystolicAccelerationPeakIndexes`.

File: `../eyeflow/EyeFlow/Parameters/DefaultEyeFlowParams.json`

Context: Confirms old default settings for smoothing and harmonic count.

Relevant lines:

- `DefaultEyeFlowParams.json:69`
  - `PulseAnalysis` section begins.
- `DefaultEyeFlowParams.json:78`
  - `"FilterSignals": true`
- `DefaultEyeFlowParams.json:81`
  - `"BandLimitedSignalHarmonicCount": 13`

### Python Implementation

File: `src/calculations/math/filtering.py`

Function: `butter_lowpass_filtfilt`

Context: Shared Python equivalent of Matlab 4th-order Butterworth + `filtfilt` low-pass smoothing.

Relevant lines:

- `filtering.py:8-14`
  - Public smoothing helper.
- `filtering.py:15-20`
  - Builds a 4th-order Butterworth low-pass filter and applies `filtfilt`.
  - Includes a short-signal guard so synthetic or tiny debug signals do not crash.
- `filtering.py:23-29`
  - Converts `dt_seconds` to normalized cutoff via Nyquist frequency, matching Matlab's `15 / (fs / 2)`.

File: `src/calculations/steps/vessel_velocity_estimator.py`

Class/function: `VesselVelocityEstimatorStep.run`

Context: The velocity estimator now filters artery and vein velocity signals before storing them in the analysis cache.

Relevant lines:

- `vessel_velocity_estimator.py:110-112`
  - Computes artery and vein velocity signals from the velocity map and masks.
- `vessel_velocity_estimator.py:113-127`
  - Applies `butter_lowpass_filtfilt` to both signals when `FilterSignals` is enabled.
- `vessel_velocity_estimator.py:127-129`
  - Stores the filtered signals as:
    - `retinal_artery_velocity_signal`
    - `retinal_vein_velocity_signal`
- `vessel_velocity_estimator.py:159-166`
  - Reads smoothing settings from `ctx.dopplerview_config["PulseAnalysis"]`.
- `vessel_velocity_estimator.py:169-174`
  - Computes `dt_seconds = batch_stride / sampling_freq`.

File: `src/pipelines/waveform_shape_metrics/constants.py`

Context: Stores current legacy defaults used by the waveform-shape sandbox.

Relevant lines:

- `constants.py:1`
  - `LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT = 13`
- `constants.py:3`
  - `LEGACY_FILTER_VELOCITY_SIGNALS = True`
- `constants.py:4`
  - `LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ = 15.0`

File: `src/pipelines/waveform_shape_metrics/dopplerview.py`

Function: `run_dopplerview_analysis`

Context: Passes Matlab-compatible smoothing settings into the DopplerView replay step.

Relevant lines:

- `dopplerview.py:36-39`
  - Supplies:
    - `FilterSignals = LEGACY_FILTER_VELOCITY_SIGNALS`
    - `LowpassFreqHz = LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ`

## 3. What Was Fixed And Why

### Beat Division Fix

The old incorrect Python behavior could produce the correct sample dimension but the wrong beat dimension because systole detection was not Matlab-equivalent.

The critical Matlab rule is:

- `numberOfBeats = numel(sysIdxList) - 1`

Therefore:

- Old Matlab output with 6 beats implies 7 systolic boundary indexes.
- Python output with 8 beats implied 9 systolic boundary indexes.

The Python detector was changed to use the same filter order, cutoff, peak threshold, peak distance, and close-peak validation used by `find_systole_index.m`.

### Smoothing Fix

The old EyeFlow `pulseAnalysis.m` filtered artery and vein velocity signals before caching and before systole/per-beat analysis.

The Python implementation now applies the same 4th-order 15 Hz Butterworth `filtfilt` smoothing to `retinal_artery_velocity_signal` and `retinal_vein_velocity_signal` before downstream beat detection and per-beat calculations.

## 4. Remaining Review Notes

- The smoothing change reproduces the Matlab signal-filtering block, but it does not make the full Python `VesselVelocityEstimatorStep` equivalent to old `pulseAnalysis.m`.
- If values still differ after beat count and smoothness align, the next review target is the upstream velocity-map/background computation:
  - Matlab: `../eyeflow/EyeFlow/BloodFlowVelocity/pulseAnalysis.m`
  - Python: `src/calculations/steps/vessel_velocity_estimator.py`
- The runtime H5 debug attributes to inspect after a run are:
  - `debug_beat_count`
  - `debug_beat_indices`
  - `debug_beat_period_idx`
  - `debug_beat_detection_min_peak_distance`
  - `debug_beat_detection_min_peak_height`
  - `filter_velocity_signals`
  - `velocity_signal_lowpass_hz`
