# Annular vessel sampling

EyeFlow configures annular sampling with one positive integer:

```json
{"number_of_radii_in_half_diagonal": 42}
```

The count includes intermediate annuli. For an image `(ny, nx)`, let
`R = max(hypot((ny - 1)/2, (nx - 1)/2), 1)` be its half-diagonal in native
pixels. For count `N`, stride is **R/N** and radial measurement thickness is
**2R/N**. The width-to-stride ratio is fixed at two, giving 50% overlap. For
512 x 512 with N=32, stride is approximately **11.29 pixels** and thickness
**22.58 pixels**. Each odd zero-based annulus is exactly halfway between its
even-indexed neighbors. Size and stride are resolved geometry, not application
preferences. GUI imports and the CLI validate the same count through the shared
run service before replacing an output directory.

`pipelines/topology_core/runner.py` resolves one geometry for prepared
artery/vein topology using the count captured in the run context. Calculation
code receives the geometry, independently of
settings persistence. Velocity, spatial gradients, BVR and displacement profiles
reuse that topology rather than deriving their own annuli.

## Coordinates and terminal sections

Internal radial coordinates are fractions of R. The optic-disc center is
unchanged. The initial boundary is
`r0 = ceil(min(disc_width, disc_height)/2)`. For center `(cx, cy)`, the last
possible start is bounded by the farthest image corner at radius:

```text
farthest = hypot(max(abs(cx), abs(nx-1-cx)), max(abs(cy), abs(ny-1-cy)))
```

Annulus `i` has full-width bounds:

```text
inner[i] = r0 + i * stride
outer[i] = inner[i] + size
membership: inner < distance_from_disc_center <= outer
```

Every start before `farthest` is retained, including tiny terminal sections.
Outer boundaries are never truncated to R or to `farthest`; only pixels outside
the rectangular image are unavailable. Branch identification also covers the
farthest corner, rather than cutting off off-center fields at R. Its support
radius is rounded outward to accommodate float32 rasterization. A float64
tolerance of `32 * eps * max(farthest-r0, size, stride, 1)` excludes numerically
zero terminal extents. If the optic-disc exclusion covers the radial extent,
there are no annuli.

N controls sampling density across a half-diagonal; it is not a promise of
exactly N measured annuli. The actual count depends on the optic-disc radius
and its position relative to the corners.

## Saved configurations and explicit calculation APIs

Application settings read `number_of_radii_in_half_diagonal` directly, with
default 42 when absent. The same scalar passes through the run specification
and worker context. JSON loading and saving do not migrate or rewrite annulus
preferences; old count names and size/stride fields are not used.

The numerical `annulus_geometry()` API still accepts explicit native-pixel
`annulus_size`/`annulus_stride`, including fractional sizes, non-overlapping
or gapped bands, and the legacy `number_of_radii_in_fov` calculation argument.
Those explicit calculation paths retain their existing clipping and terminal
tail rules. They are not exposed as application settings.

Counts embedded in source HD/DV HDF5 or sidecar metadata are acquisition
attributes, not an intentional EyeFlow sampling choice, and are ignored.

## Measurements and products

`delta_radius` is `outer - inner` in native pixels, so count-based sampling
retains the full nominal thickness even in outer annuli. It is never the stride.
Mask-derived lumen diameter remains native mask area
divided by this thickness. Physical diameter uses pixel pitch, and both mask
and gradient BVR retain their established physical formulas and mm^3/s units.
Total mask BVR still sums branches and takes a median over radii after its
existing temporal filter. No overlap correction or extra statistical weights
are introduced, and lumen histograms retain count semantics.

The larger radial axis naturally affects segment masks/centers, empirical
window-size percentiles, lumen samples and QC, gradient edge profiles, velocity
signals/maps/profiles/fits/FFT, waveform metrics, quadrant/radius aggregations,
BVR and optional displacement profiles. These observations overlap and are
correlated; more samples do not imply more independent observations. Aggregate
values and window sizes may change naturally. Global calculations independent
of annuli and full-field displacement registration, fields, magnitude maps and
videos are unchanged.

The annulus diagnostic images cycle outline shades over the overlapping stride
phases (`ceil(size / stride)`, bounded by the available annuli). At the default
50% overlap, annulus 1 is white (255), annulus 2 gray (128), annulus 3 white,
and so on. Denser overlap adds intermediate gray levels; adjacent/gapped
sampling has a single white phase. Both boundaries use the same annulus shade.
Shared boundaries and raster collisions retain the earliest annulus's
shade. The innermost optic-disc boundary is white. These outlines are used in
diagnostic images; annulus geometry and outline intensity are not written as a
separate `/Segmentation/Annuli/` HDF5 group. Branch/segment label maps retain
their existing numeric label codes; the grayscale image is a separate
visualization, not a measurement mask. Velocity provenance describes the
resolved annuli.

Sampling density and the farthest corner determine annular work and radial array sizes.
Existing CPU/GPU temporal chunk budgets, sparse/compressed HDF5 storage and
artifact resource management remain in use; no new execution path or tuning
constraints are introduced. Full-resolution acquired-data benchmarking is
separate from synthetic regression validation.
