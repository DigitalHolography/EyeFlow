# Spectral pipeline performance

This is a dated local measurement, not a run-time guarantee. The fixture was a
synthetic, uncompressed `(640, 64, 128, 128)` float32 `/spectrograms` dataset
(2.5 GiB), with a synthetic 75% gap mask and eight beats. Measurements were
made on Windows with 16 logical CPUs on 2026-10-28. They exclude upstream
velocity/topology work and do not measure NAS or a clinical run.

| Version measured on the fixture | Spatial pass | Spectral stages including figures |
|---|---:|---:|
| Original serial reducer | 34.33 s | 48.47 s |
| Bounded native-coordinate reducer | 14.91 s | 27.87 s |
| Eight-worker reduction | 3.38 s | 16.50 s |

The implemented reducer uses a 64 MiB read buffer, native HD spatial order,
exact spatial medians, float64 mean accumulation, and at most eight available
CPU workers for sufficiently large blocks. All 149 persisted datasets in the
measured fixture matched the preceding implementation exactly. In the
measured eight-worker run, figure export took 11.32 s and dominated local
execution. Current stage-graph scheduling and shared background rendering
can overlap exports with later work; the table predates those execution paths.

The earlier full-size serial pass made 640 HDF5 reads. Bounded time slabs
reduced this to 40 calls for this fixture. Call count does not imply a NAS
transfer count or speedup. A separate isolated export comparison found
pixel-identical PNGs and matching EPS content after date/title metadata was
excluded, but did not measure a full run. Measure representative source files
through final HDF5 publication before changing resource limits or claiming
end-to-end gains.

The checked-in `benchmarks/profile_low_freq_spectra.py` can profile the
spectral stages; its `--help` gives current options. Use a fresh output
directory; full-size synthetic input generation needs at least 2.5 GiB of
space. The current implementation and scientific rules are in
`src/pipelines/low_freq_spectra/` and
[data contracts](data-contracts.md).
