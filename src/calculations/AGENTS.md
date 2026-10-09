# Scientific calculations context

## Boundary

`calculations` contains numerical and transform logic reusable by more than one
pipeline. It may depend on NumPy/SciPy/scikit-image and the optional compute
backend, but it must not know about Tkinter, pipeline selection, persisted
settings, `.holo` layout, or output-directory policy. Logic whose purpose is
one product belongs in that pipeline package; pipelines translate typed source
data into reusable calculation inputs and pack results into the EyeFlow schema.

## Context scopes

| Scientific task | Start with | Important consumers/tests |
|---|---|---|
| Cardiac-cycle detection | `blood_flow_velocity/signal_analysis/cardiac_cycle/` | `pipelines/velocity/cardiac_cycle.py`; `test_cardiac_cycle_analysis.py`, `test_velocity_pipeline.py` |
| Per-beat resampling | `blood_flow_velocity/signal_analysis/per_beat/` | velocity analysis and waveform metrics; `test_per_beat_runner.py` |
| Waveform morphology/metrics | `blood_flow_velocity/signal_analysis/waveform/` | shape/absolute/low-rank pipelines; their focused tests |
| Branch, annulus, segment, or cross-section behavior | `topology/` | topology core and velocity/spatial-gradient profiles; all `test_topology_*`, profile, and segment-output tests |
| Spatial gradient/filter behavior | `math/spatial_gradient.py`, `math/image.py`, periodic windows | spatial-gradient pipeline; `test_spatial_gradient_*`, `test_gaussian2d_blur.py`, `test_unsharpen.py`, `test_periodic_sliding_windows.py` |
| Blood-volume-rate formulas | `blood_volume_rate.py` | `pipelines/blood_volume_rate/`; `test_blood_volume_rate.py`, pipeline dependency test |
| Shared statistics/Fourier utilities | `math/` | search call sites; `test_calculations_math.py` and affected pipeline tests |
| CPU/GPU behavior | `compute_backend.py`, topology transforms/workflow | topology transform/chunk tests and benchmark snapshot |

## Velocity calculation ownership

Retinal velocity estimation is the literal purpose of the `velocity` pipeline,
so its product-specific estimator lives in
`pipelines/velocity/estimation.py`, with mask preparation in
`pipelines/velocity/masks.py`. Read the pipeline context and data contracts for
its method, calibration, and output invariants. Reusable cardiac-cycle,
per-beat, waveform, topology, and math transforms remain in `calculations`.

## Topology and cross-sections

`topology` owns optic-disc geometry, annulus construction, branch identity,
segment masks/centers, orientation, transforms, profiles, native mask areas,
cache keys, and bounded chunk preparation. The optic-disc center defines radial
direction. Its mask is removed from vessel support before topology products.

Velocity segments use fused resize-and-rotate sampling. The dormant
`pipelines/displacement_map/segments.py` helper can use the same topology
transforms, but the current displacement runner does not call it. Spatial
gradients deliberately use staged processing:

```text
interpolate -> centered temporal mean -> Sobel magnitude
 -> ImageJ Gaussian -> ImageJ unsharp -> centered temporal mean -> rotate
```

The gradient path carries periodic temporal halo data across chunks so results
match whole-recording filtering. The chunk planner applies one scratch-memory
budget across workers; retained result arrays are not counted. Changes to
orientation, interpolation, mask order, halos, or branch identity affect
multiple consumers and require topology, velocity-profile, spatial-gradient,
and segment-output tests as appropriate. Add displacement coverage if its
segment helper becomes part of the production runner.

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
