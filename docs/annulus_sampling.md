# Annular vessel sampling

EyeFlow's annulus size defines radial measurement support; its stride defines
the displacement between starting boundaries. Settings are native pixels and
can be fractional. In `default_settings.json`, both values are `null`, meaning:

```json
{"annulus_size": null, "annulus_stride": null}
```

For an image `(ny, nx)`, size defaults to `max(nx, ny) / 25`, exactly the old
default thickness. A missing/null stride is half the resolved size, including
when a user selects a different size. For 512 x 512, size is **20.48 pixels**
and stride is **10.24 pixels**. Both explicit values must be finite and positive.
Equal size/stride requests adjacent sections; stride greater than size permits
gaps. GUI configuration imports and the CLI use the same saved settings and
run-service validation, before replacing an output directory.

`pipelines/topology_core/annulus_sampling.py` resolves one geometry for prepared
artery/vein topology. Calculation code receives the geometry, independently of
settings persistence. Velocity, spatial gradients, BVR and displacement profiles
reuse that topology rather than deriving their own annuli.

## Coordinates and terminal sections

The native-pixel radius scale remains
`max(hypot((ny - 1)/2, (nx - 1)/2), 1)`. Internal radial coordinates are fractions
of that scale. The optic-disc center is unchanged. The initial boundary is
`min(ceil(min(disc_width, disc_height)/2), radius_scale)` and the outer limit is
`radius_scale`. Annulus `i` has bounds:

```text
inner[i] = initial + i * stride
outer[i] = min(radius_scale, inner[i] + size)
membership: inner < distance_from_disc_center <= outer
```

Count is derived from the starts before the outer limit. Partial final annuli
remain allowed. Only an additional overlapping terminal slice thinner than
`min(1 native pixel, 0.05 * size, stride)` is omitted. The cap avoids subpixel
slivers at ordinary sizes; the relative threshold preserves intentional fine
sampling. At 50% overlap, an even zero-based index corresponds to a legacy
section and retains its clipped tail even below this threshold. Non-overlapping
and gapped configurations also retain their final partial section. A float64
tolerance of `32 * eps * max(extent, size, stride, 1)` excludes numerical
zero-width candidates; an empty radial extent yields no annuli.

All legacy default sections remain at odd one-based indices. An extent of
`N * size` has `2N-1` complete overlapping sections, followed by a useful clipped
half-width section. Clipping is also applied when drawing annulus boundaries.

## Explicit legacy configurations

An older **EyeFlow settings file** containing only `number_of_radii_in_FOV`
(or an accepted historical spelling) resolves to `max(nx, ny)/count` for both
size and stride. This preserves intentional historical non-overlapping sampling
without misinterpreting a count as pixels. This compatibility input remains
dynamic across image sizes and is absent from new default configurations.
Any explicit non-null new size or stride takes precedence over all legacy
count fields; the other dimension uses its new default. To opt an old file
into overlapping defaults, remove its legacy count field.

Counts embedded in source HD/DV HDF5 or sidecar metadata are acquisition
attributes, not an intentional EyeFlow sampling choice, and are ignored.

## Measurements and products

`delta_radius` is `outer - inner` in native pixels, including clipped widths;
it is never the stride. Mask-derived lumen diameter remains native mask area
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

Resolved geometry is published once under `/Segmentation/Annuli/`:

| Dataset suffix | Meaning |
| --- | --- |
| `Size/value` | Nominal thickness in native pixels |
| `Stride/value` | Starting-boundary displacement in native pixels |
| `Count/value` | Derived annulus count |
| `BoundsPixels/value` | `(radius, boundary)` actual native-pixel bounds |
| `BoundsFraction/value` | Matching fractional half-diagonal bounds |
| `RadiusScalePixels/value` | Conversion scale in pixels |
| `OutlineIntensity/value` | Shared uint8 monochrome outline image, lower-left `(x, y)` |

The annulus diagnostic images cycle outline shades over the overlapping stride
phases (`ceil(size / stride)`, bounded by the available annuli). At the default
50% overlap, annulus 1 is white (255), annulus 2 gray (128), annulus 3 white,
and so on. Denser overlap adds intermediate gray levels; adjacent/gapped
sampling has a single white phase. Both boundaries use the same annulus shade.
Shared or clipped boundaries and raster collisions retain the earliest annulus's
shade. The innermost optic-disc boundary is white. `OutlineIntensity` stores the
same overlay on a black (0) background with palette and boundary-policy
attributes. Branch/segment label maps retain their existing numeric label codes;
the grayscale image is a separate visualization, not a measurement mask.

Boundary columns are inner then outer; indices are zero-based. Geometry is
serialized with the existing float32 numeric policy. `DeltaRadius` retains its
existing NaN convention for sections without any vessel support. The explicit
boundary datasets remain available even for empty vessel sections. Velocity
provenance also describes the actual resolved annuli.

Default sampling roughly doubles annular work and retained radial arrays.
Existing CPU/GPU temporal chunk budgets, sparse/compressed HDF5 storage and
artifact resource management remain in use; no new execution path or tuning
constraints are introduced. Full-resolution acquired-data benchmarking is
separate from synthetic regression validation.
