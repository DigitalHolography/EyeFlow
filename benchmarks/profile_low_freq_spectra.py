"""Profile spectral stages without replacing an existing EyeFlow result.

Input HDF5 files are opened read-only. All derived files go to --output.
This is an isolated spectral benchmark: supplied masks/timing replace upstream
velocity/topology work. Synthetic CPU-only runs do not measure file I/O;
synthetic HDF5 runs measure the generated storage layout. Neither establishes
the performance of an actual HD export. Run in an editable environment:

python benchmarks/profile_low_freq_spectra.py --hd scan_HD_output.h5 --frames 32
python benchmarks/profile_low_freq_spectra.py --synthetic-shape 64 64 256 256
python benchmarks/profile_low_freq_spectra.py --generate-synthetic-hd 640 64 128 128 --output .tmp/spectral_profile_new
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import ExitStack
import cProfile
import csv
from functools import wraps
import io
import json
from pathlib import Path
import platform
import pstats
import shutil
from time import perf_counter
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

from input_output.h5_access import PipelineH5Output, RawH5SourceReader
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager
from input_output.schema import EyeFlowOutputPaths, HolodopplerSource
from pipelines.low_freq_spectra import artifacts, gap_analysis, qc, runner, spectrograms
from pipelines.low_freq_spectra.frequency import write_frequency_bin_table
from pipelines.low_freq_spectra.preprocessing import prepare_gap_spectrograms
from pipelines.vessel_inputs import align_spectrogram_frame


class Timings:
    def __init__(self):
        self.seconds = defaultdict(float)
        self.calls = defaultdict(int)

    def wrap(self, name, function):
        @wraps(function)
        def measured(*args, **kwargs):
            start = perf_counter()
            try:
                return function(*args, **kwargs)
            finally:
                self.seconds[name] += perf_counter() - start
                self.calls[name] += 1
        return measured

    def run(self, name, function, *args, **kwargs):
        return self.wrap(name, function)(*args, **kwargs)

    def rows(self):
        return [dict(stage=name, seconds=self.seconds[name], calls=self.calls[name])
                for name in sorted(self.seconds, key=self.seconds.get, reverse=True)]


class Source:
    """Count actual source accesses; keep at most one synthetic frame in memory."""
    def __init__(self, dataset=None, *, shape=None, frame_limit=None):
        self.dataset = dataset
        self.shape = tuple(dataset.shape if dataset is not None else shape)
        if frame_limit:
            self.shape = (min(frame_limit, self.shape[0]), *self.shape[1:])
        self.ndim, self.dtype = 4, np.dtype(np.float32) if dataset is None else dataset.dtype
        self.attrs = {"unit": "a.u."} if dataset is None else dict(dataset.attrs)
        self.seconds, self.reads, self.logical_bytes = 0., 0, 0
        if dataset is None:
            generator = np.random.default_rng(42)
            self.frame = generator.random(self.shape[1:], dtype=np.float32) + .1
            # Fixed low-frequency peaks, with spatially variable real power.
            weights = 1 + 5*np.exp(-((np.arange(self.shape[1])-3)/3)**2)
            weights += 4*np.exp(-((np.arange(self.shape[1])-(self.shape[1]-4))/3)**2)
            self.frame *= weights[:, None, None].astype(np.float32)

    def __getitem__(self, index):
        start = perf_counter()
        data = self.dataset[index] if self.dataset is not None else self.frame * np.float32(1+.1*np.sin(index))
        if data.dtype != self.dtype:
            raise TypeError("Benchmark source returned a different dtype than it advertised.")
        self.seconds += perf_counter() - start
        self.reads += 1
        self.logical_bytes += data.nbytes
        return data


def save_profile(profile, output, name):
    profile.dump_stats(str(output / f"{name}.prof"))
    stream = io.StringIO()
    pstats.Stats(profile, stream=stream).strip_dirs().sort_stats("cumulative").print_stats(45)
    (output / f"{name}.txt").write_text(stream.getvalue(), encoding="utf-8")


def spatial_breakdown(source, zones, count):
    """Time the same per-frame operations separately; outside the main run."""
    timer = Timings()
    prepared = prepare_gap_spectrograms(source, zones)
    for index in range(min(count, source.shape[0])):
        raw = timer.run("source_access", source.__getitem__, index)
        frame = timer.run("align_and_cast_float64", lambda: np.asarray(align_spectrogram_frame(raw), dtype=np.float64))
        full = timer.run("reshape_spatial_axes", frame.reshape, frame.shape[0], -1)
        valid = timer.run("full_finite_mask", lambda: np.all(np.isfinite(full), axis=1))
        mean_samples = timer.run("full_mean_boolean_copy", lambda: full[valid])
        timer.run("full_mean", np.mean, mean_samples, axis=1)
        median_samples = timer.run("full_median_boolean_copy", lambda: full[valid])
        timer.run("full_median", np.median, median_samples, axis=1)
        if prepared is not None:
            samples = timer.run("gap_boolean_gather", lambda: frame[:, prepared[0]])
            gap_valid = timer.run("gap_finite_mask", lambda: np.all(np.isfinite(samples), axis=1))
            median_samples = timer.run("gap_median_boolean_copy", lambda: samples[gap_valid])
            timer.run("gap_median", np.median, median_samples, axis=1)
    return timer.rows()


def hd_metadata(path):
    with h5py.File(path, "r") as h5:
        hd = HolodopplerSource(RawH5SourceReader(h5file=h5, label="HD"))
        ctx = SimpleNamespace(require_inputs=lambda *args: None,
                              inputs=SimpleNamespace(hd=SimpleNamespace(as_holodoppler=lambda: hd)))
        source = hd.optional_spectograms_dataset()
        if source is None:
            start = perf_counter()
            runner.run_low_freq_spectra(ctx)
            return {"path": str(path), "spectrograms_available": False,
                    "top_level_keys": list(h5), "skipped_pipeline_seconds": perf_counter()-start}
        return {"path": str(path), "spectrograms_available": True, "shape": source.shape,
                "dtype": str(source.dtype), "chunks": source.chunks, "compression": source.compression,
                "compression_opts": source.compression_opts, "logical_bytes": source.nbytes,
                "stored_bytes": source.id.get_storage_size(),
                "chunk_cache": source.id.get_access_plist().get_chunk_cache()}


def compare_gap_gather(source, zones):
    """Microbenchmark a candidate using the already-created flat frame.

    Does not replace production code. Exact equality establishes that the
    two indexing expressions select identical values in identical order.
    """
    mask = prepare_gap_spectrograms(source, zones)[0]
    frame = np.asarray(align_spectrogram_frame(source[0]), dtype=np.float64)
    full = frame.reshape(frame.shape[0], -1)
    flat_mask = mask.ravel()
    np.testing.assert_array_equal(frame[:, mask], full[:, flat_mask])
    samples = {"current_3d_boolean_gather": [], "candidate_existing_flat_frame_gather": []}
    for repetition in range(6):
        order = list(samples) if repetition % 2 == 0 else list(reversed(samples))
        for name in order:
            start = perf_counter()
            selected = frame[:, mask] if name.startswith("current") else full[:, flat_mask]
            samples[name].append(perf_counter()-start)
            del selected
    return {"measurement": "alternating_order_warm_memory_microbenchmark_not_end_to_end_speedup",
            "arrays_exactly_equal": True,
            **{name: float(np.median(seconds)) for name, seconds in samples.items()}}


def generate_synthetic_hd(path, shape):
    """Create a fully allocated contiguous float32 HD fixture, one frame at a time."""
    shape = tuple(shape)
    required = int(np.prod(shape))*np.dtype(np.float32).itemsize
    if shutil.disk_usage(path.parent).free < required + 1024**3:
        raise OSError("Not enough disk space for the synthetic source and a 1 GiB reserve.")
    source = Source(shape=shape)
    start = perf_counter()
    with h5py.File(path, "x") as h5:
        h5.create_dataset("HD_parameters", data=json.dumps({
            "sampling_freq": 37037., "batch_size": shape[1]*4, "batch_stride": 256,
            "spectrogram_frequency_binning": 4,
        }))
        dataset = h5.create_dataset("spectrograms", shape=shape, dtype=np.float32)
        dataset.attrs["unit"] = "a.u."
        for frame in range(shape[0]):
            dataset[frame] = source[frame]
        h5.flush()
    return {"seconds": perf_counter()-start, "logical_bytes": required,
            "storage": "contiguous_uncompressed_fully_allocated",
            "cache_note": "Generated immediately before the run; OS cache may be warm. Not a cold-disk benchmark."}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--hd", type=Path)
    mode.add_argument("--synthetic-shape", type=int, nargs=4, metavar=("TIME", "FREQ", "Y", "X"))
    mode.add_argument("--generate-synthetic-hd", type=int, nargs=4, metavar=("TIME", "FREQ", "Y", "X"))
    parser.add_argument("--frames", type=int, help="Use only the first N source frames")
    parser.add_argument("--output", type=Path, default=Path(".tmp/spectral_performance"))
    parser.add_argument("--spatial-detail-frames", type=int, default=4)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summary = {"mode": "real_HD_isolated_pipeline" if args.hd else "synthetic_CPU_only",
               "python": platform.python_version(), "numpy": np.__version__, "h5py": h5py.__version__,
               "timing_note": "Single wall-time run. Nested detailed timers are inclusive and must not be summed with top-level stages.",
               "upstream_note": "All-pixel mask with 25% regularly excluded; synthetic evenly spaced beat boundaries. Not upstream velocity/topology execution."}
    if args.generate_synthetic_hd:
        args.hd = args.output / "synthetic_HD.h5"
        summary["mode"] = "synthetic_hdf5_isolated_pipeline"
        print(f"Generating synthetic HDF5: {tuple(args.generate_synthetic_hd)}", flush=True)
        summary["source_generation"] = generate_synthetic_hd(args.hd, args.generate_synthetic_hd)
    if args.hd:
        summary["hd"] = hd_metadata(args.hd)
        if not summary["hd"]["spectrograms_available"]:
            (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
            print(json.dumps(summary, indent=2), flush=True)
            return
    timer = Timings()
    with ExitStack() as stack:
        dataset, hd = None, None
        if args.hd:
            h5 = stack.enter_context(h5py.File(args.hd, "r"))
            dataset = h5["spectrograms"]
            hd = HolodopplerSource(RawH5SourceReader(h5file=h5, label="HD"))
        source = Source(dataset, shape=args.synthetic_shape, frame_limit=args.frames)
        summary["shape"] = source.shape
        if source.shape[0] < 4 or source.shape[1] < 2:
            raise ValueError("Use at least four frames and two frequency bins.")
        zones = np.ones(source.shape[-2:], dtype=bool)
        zones[::2, ::2] = False
        manager = OutputManager(HoloRunLayout.from_holo(args.output / "profile.holo"))
        work = stack.enter_context(h5py.File(args.output / "profile.h5", "w"))
        output = SimpleNamespace(h5=PipelineH5Output(work), manager=manager, available=True)
        ctx = SimpleNamespace(output=output, attrs={})
        schema = EyeFlowOutputPaths.active()
        fs = hd.sampling_frequency_hz() if hd else 37037.
        output.h5.write(schema.sampling_frequency_hz, fs, unit="Hz")
        output.h5.write(schema.segmentation.gaps, np.flip(zones, axis=0).T)
        if hd is None:
            hd = SimpleNamespace(spectrogram_frequency_parameters=lambda count: {
                "fft_bin_count": 4*count, "frequency_binning_factor": 4,
                "frequency_selection": "full_FFT_domain"})
        write_frequency_bin_table(output, source.shape[1], hd=hd)
        boundaries = np.unique(np.linspace(0, source.shape[0]-1, min(9, source.shape[0]), dtype=int))
        timing = SimpleNamespace(
            cycle_boundary_indexes=boundaries, cardiac_cycle_source="artery",
            cycle_durations_seconds=np.diff(boundaries)*.02,
            continuous=lambda *args, **kwargs: np.zeros(source.shape[0]),
        )
        print(f"Profiling {summary['mode']}, shape={source.shape}", flush=True)
        for name in ("gap_spectrograms_per_beat", "nanmedian", "beat_spectrum_std", "pack_spectrum_metrics",
                     "pack_vertical_metrics", "pack_distribution_metrics", "band_limit_gap_spectrograms",
                     "svd_filter_gap_spectrograms"):
            stack.enter_context(patch.object(gap_analysis, name, timer.wrap("detail/"+name, getattr(gap_analysis, name))))
        start = perf_counter()
        outputs = timer.run("spatial_preprocessing", spectrograms.preprocess_spectrograms, source, zones)
        print(f"Spatial preprocessing: {timer.seconds['spatial_preprocessing']:.3f}s; source access: {source.seconds:.3f}s", flush=True)
        timer.run("write_spatial_outputs", output.h5.write_many, outputs)
        products = timer.run("gap_analysis", gap_analysis.run_gap_analysis, ctx, timing)
        timer.run("write_gap_outputs", output.h5.write_many, products)
        gap = output.h5.get(schema.spectrograms_masked)
        quality = timer.run("qc_calculation", qc.pack_gap_qc, gap[()], source_attrs=dict(gap.attrs), gap_outputs=products)
        timer.run("write_qc_outputs", output.h5.write_many, quality)
        if not args.no_figures:
            from matplotlib.figure import Figure
            original_pair, original_save = artifacts._save_figure_pair, Figure.savefig

            def pair(manager, name, make_figure):
                return timer.run("artifact/"+name, original_pair, manager, name, make_figure)

            def save(fig, filename, *args, **kwargs):
                suffix = kwargs.get("format") or Path(str(filename)).suffix.lstrip(".")
                return timer.run("save_format/"+suffix, original_save, fig, filename, *args, **kwargs)

            stack.enter_context(patch.object(artifacts, "_save_figure_pair", pair))
            stack.enter_context(patch.object(Figure, "savefig", save))
            paths = timer.run("figure_export", artifacts.export_spectrogram_figures, output, dt_seconds=.02)
            summary["artifact_count"] = len(paths)
        timer.run("final_h5_flush", work.flush)
        summary.update(timed_pipeline_seconds=perf_counter()-start,
                       source_access_seconds=source.seconds, source_read_calls=source.reads,
                       source_logical_bytes=source.logical_bytes, timings=timer.rows())
        # Separate profile passes; their cost does not enter pipeline timings.
        detail_source = Source(dataset, shape=args.synthetic_shape, frame_limit=args.frames)
        summary["spatial_operation_timings"] = spatial_breakdown(detail_source, zones, args.spatial_detail_frames)
        summary["gap_gather_candidate"] = compare_gap_gather(detail_source, zones)
        profiler = cProfile.Profile()
        profiler.runcall(gap_analysis.run_gap_analysis, ctx, timing)
        save_profile(profiler, args.output, "gap_analysis")
        profiler = cProfile.Profile()
        profiler.runcall(qc.pack_gap_qc, gap[()], source_attrs=dict(gap.attrs), gap_outputs=products)
        save_profile(profiler, args.output, "qc")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (args.output / "timings.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=("stage", "seconds", "calls"))
        writer.writeheader()
        writer.writerows(summary["timings"])
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
