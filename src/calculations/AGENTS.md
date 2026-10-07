# Scientific calculations context

## Boundary

`calculations` contains reusable numerical and domain code. It may depend on
NumPy/SciPy/scikit-image and the optional compute backend, but it must not know
about Tkinter, pipeline selection, persisted settings, `.holo` layout, or
output-directory policy. Pipeline packages translate typed source data into
calculation inputs and pack results into the EyeFlow schema.

## Context scopes

| Scientific task | Start with | Important consumers/tests |
|---|---|---|
| Local retinal velocity | `retinal_velocity/vessel_velocity_estimator.py`, `_masks.py` | `heartbeat_core`, `waveform_velocity_core`; `test_frequency_band_velocity.py`, `test_heartbeat_core.py` |
| Heartbeat and waveform boundaries | `blood_flow_velocity/signal_analysis/heartbeat/` | heartbeat core and all per-beat products; `test_heartbeat_analysis.py`, `test_per_beat_runner.py` |
| Per-beat resampling | `blood_flow_velocity/signal_analysis/per_beat/` | waveform core/metrics; `test_per_beat_runner.py` |
| Waveform morphology/metrics | `blood_flow_velocity/signal_analysis/waveform/` | shape/absolute/low-rank pipelines; their focused tests |
| Branch, annulus, segment, or cross-section behavior | `topology/` | topology core, velocity/spatial-gradient/displacement profiles; all `test_topology_*`, profile, and segment-output tests |
| Spatial gradient/filter behavior | `math/spatial_gradient.py`, `math/image.py`, periodic windows | spatial-gradient pipeline; `test_spatial_gradient_*`, `test_gaussian2d_blur.py`, `test_unsharpen.py`, `test_periodic_sliding_windows.py` |
| Blood-volume-rate formulas | `blood_volume_rate.py` | `pipelines/blood_volume_rate/`; `test_blood_volume_rate.py`, pipeline dependency test |
| Shared statistics/Fourier utilities | `math/` | search call sites; `test_calculations_math.py` and affected pipeline tests |
| CPU/GPU behavior | `compute_backend.py`, topology transforms/workflow | topology transform/chunk tests and benchmark snapshot |

## Velocity calculation contract

The estimator processes `(frame, y, x)` volumes in fixed-size temporal chunks.
It builds one dilated vessel-background mask, inpaints each chunk, computes a
signed difference from local background, and returns maps plus artery/vein
signals. The active source differs by method:

- moments: `sqrt(moment2 / spatial_mean(moment0))`, guarded where the mean is
  zero, then physical scaling after background subtraction;
- bands: `HF / LF`, with exact-zero LF mapped to zero, then the same background
  path after Hz-per-ratio calibration and with physical mm/s scaling.

Inputs must be finite and non-negative in band mode. Do not introduce epsilon
bias, infinity, or a silent fallback. The estimator cache key includes the
method and only the active source identities so heartbeat and waveform cores can
reuse exactly matching work.

## Topology and cross-sections

`topology` owns optic-disc geometry, annulus construction, branch identity,
segment masks/centers, orientation, transforms, profiles, native mask areas,
cache keys, and bounded chunk preparation. The optic-disc center defines radial
direction. Its mask is removed from vessel support before topology products.

Velocity/displacement segments normally use fused resize-and-rotate sampling.
Spatial gradients deliberately use staged processing:

```text
interpolate -> centered temporal mean -> Sobel magnitude
 -> ImageJ Gaussian -> ImageJ unsharp -> centered temporal mean -> rotate
```

The gradient path carries periodic temporal halo data across chunks so results
match whole-recording filtering. The chunk planner applies one scratch-memory
budget across workers; retained result arrays are not counted. Changes to
orientation, interpolation, mask order, halos, or branch identity affect
multiple consumers and require topology, velocity-profile, spatial-gradient,
displacement, and segment-output tests as appropriate.

The current production chain uses two centered 7-frame means, ImageJ Gaussian
radius 6, and ImageJ unsharp radius 8 with weight 0.6; the staged chunk halo is
six frames on each side. `math/spatial_gradient.py` and the reference tests are
the source of truth for these constants and border/NaN behavior.

## Numerical invariants

- Treat axis order, sign, units, index base, and NaN behavior as API.
- Do not replace NaN-aware reducers with ordinary reducers or vice versa without
  tracing the complete output contract.
- Preserve float32 persisted/result arrays unless a solver explicitly uses
  float64 internally for stability.
- Boundary behavior is intentional: periodic waveform windows, shrinking or
  halo-backed filters, and spatial border zeroing are distinct contracts.
- Masks are boolean support, not weights. Functions should not mutate caller
  arrays unless documented.
- Validate physical calibration before applying it. `PixelPitch.isotropic_*`
  intentionally rejects materially anisotropic sampling.

For detailed weighted quadratic profile behavior, read
[velocity-profile analysis](../../docs/velocity_profile_analysis.md); its solver
lives with the pipeline because it is specific to that output product.
