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
| Native branch, annulus, segment geometry | `topology/` | topology core; segment, orientation, mask-area, quadrant, and optic-disc tests |
| Segment sampling, measurement, or cross-section behavior | `vessel_segments/` | topology core and velocity/spatial-gradient profiles; cache, workflow, transform, chunk, profile, and segment-output tests |
| Spatial gradient/filter behavior | `math/spatial_gradient.py`, `math/image.py`, periodic windows | spatial-gradient pipeline; `test_spatial_gradient_*`, `test_gaussian2d_blur.py`, `test_unsharpen.py`, `test_periodic_sliding_windows.py` |
| Blood-volume-rate formulas | `blood_volume_rate.py` | `pipelines/blood_volume_rate/`; `test_blood_volume_rate.py`, pipeline dependency test |
| Shared statistics/Fourier utilities | `math/` | search call sites; `test_calculations_math.py` and affected pipeline tests |
| CPU/GPU behavior | `compute_backend.py`, `vessel_segments/sampling/transforms.py` and `streaming.py` | topology transform/chunk tests and benchmark snapshot |

## Velocity calculation ownership

Retinal velocity estimation is the literal purpose of the `velocity` pipeline,
so its product-specific estimator lives in
`pipelines/velocity/estimation.py`, with mask preparation in
`pipelines/velocity/masks.py`. Read the pipeline context and data contracts for
its method, calibration, and output invariants. Reusable cardiac-cycle,
per-beat, waveform, topology, and math transforms remain in `calculations`.

## Native topology and sampled vessel segments

`topology` owns only source-coordinate geometry: optic-disc geometry, annuli,
branch identity, segment masks/centers/windows, centerline orientation, native
mask areas, and anatomical quadrants. The optic-disc center defines radial
direction. Its mask is removed from vessel support before topology products.
It must not import `vessel_segments`, measurement execution, or the compute backend.

`vessel_segments` owns analysis of maps sampled using that geometry:

| Package/module | Responsibility |
|---|---|
| `sampling/models.py` | `SegmentSamplingPlan`, sampled scalar/vector segments, and chunks |
| `sampling/preparation.py` | Spatial plans, shared window sizes, reference-based orientation fallback |
| `sampling/extraction.py` | Padded scalar or component-valued map patches |
| `sampling/transforms.py` | CPU/GPU interpolation, rotation, and mask transforms |
| `sampling/streaming.py` | Ordered temporal chunks, halos, and scratch-memory planning |
| `sampling/cache.py` | Pure plan-cache identity; run-state ownership stays in `pipelines/topology_core/cache.py` |
| `measurement/` | Scalar-map measurement settings, streaming accumulation, and `SegmentMeasurements` |
| `profiles/reductions.py`, `per_beat.py` | Array-based spatial reductions and temporal resampling |
| `profiles/fft.py` | Transverse FFT magnitude profiles and `SegmentFftAccumulator` |
| `profiles/fits/quadratic.py` | Weighted quadratic fits and downward-opening geometry |

Dependency direction is measurement -> sampling + profile reductions, and
sampling -> native topology. Profile reductions and fits do not consume
topology objects; streamed FFT support may use sampling mask transforms.
Keep package initializers narrow, particularly so importing the quadratic
fitter does not load spatial execution or optional GPU machinery.

A sampling plan may use a reference map to resolve indeterminate orientation,
but contains spatial data only, not temporal execution settings. Cache identity
describes the segmentation-derived plan before that refinement; changing cache
semantics is a separate task. Sampling preserves scalar/vector axes, whereas
the measurement runner currently accepts scalar maps only.

Canonical callers use `SegmentSamplingPlan`, `SegmentMeasurements`,
`SegmentMeasurementSettings`, and `fit_quadratic_profiles`. Transitional aliases
live in their new owning packages, never in `topology`. The existing result
field `SegmentMeasurements.topology` holds a sampling plan; `.topology.native`
provides native geometry. Runtime state keys and output schema remain unchanged.

Velocity segments use fused resize-and-rotate sampling. The dormant
`pipelines/displacement_map/segments.py` helper can use the same segment sampling
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
[velocity-profile analysis](../../docs/velocity_profile_analysis.md). The generic
solver lives in `vessel_segments/profiles/fits/quadratic.py`; the velocity-analysis
pipeline owns source-path selection, provenance, and output packing.
