"""Frequency-resolved Doppler CVI analysis and binary occupancy products.

This product module follows gap_analysis.py: it owns input requirements,
calculations and durable result definitions. Numerical helpers operate on
arrays; the stage runner handles scheduling and publication.
"""

from __future__ import annotations

import warnings
from numbers import Real

import numpy as np
import skimage
from skimage.filters import threshold_otsu

from calculations.math.cycle_boundaries import normalize_cycle_boundaries
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.eyeflow_output import VELOCITY_WORKFLOW_ROOTS, processing_path
from input_output.schema.holodoppler import HD_SPECTOGRAMS_PATH
from utils.logger import Logger

from .frequency import read_frequency_bin_table
from .per_beat import DEFAULT_PHASE_SAMPLE_COUNT, align_power


def run_cvi_analysis(ctx, timing, *, timing_method, **parameters) -> dict[str, object]:
    """Analyze signed source bins using detected beats and describe durable CVI results.

    The stage runner publishes the returned datasets; this analysis owns their
    timing, calibration, scientific definitions and provenance.
    """
    source = ctx.inputs.hd.as_holodoppler().optional_spectograms_dataset()
    if source is None or timing is None or timing.cardiac_cycle_source == "none" or timing.cycle_boundary_indexes.size < 2:
        Logger.log_warning("Skipping Doppler CVI proxy: spectrograms or detected complete beats are unavailable.")
        return {}
    schema = EyeFlowOutputPaths.active()
    frequencies = read_frequency_bin_table(ctx.output, source.shape[1])
    if frequencies is None:
        Logger.log_warning("Skipping Doppler CVI proxy: calibrated signed FrequencyBinTable is unavailable.")
        return {}
    if source.shape[0] != timing.continuous("artery", raw=True).size:
        raise ValueError("CVI spectra and selected cardiac timing must share the same frame count.")
    result = calculate_cvi(
        source, timing.cycle_boundary_indexes, frequencies[()],
        **parameters,
    )
    power_unit = getattr(source, "attrs", {}).get("unit", "a.u.")
    root = schema.spectrograms_cvi_root
    attrs = dict(
        quantity="frequency_resolved_Doppler_CVI_proxy", label="Doppler CVI proxy",
        specificity="functional_power_proxy_not_anatomical_OCT_CVI_or_proven_choroidal_specificity",
        source_path="/" + HD_SPECTOGRAMS_PATH,
        source_power_unit=power_unit, source_dtype=str(source.dtype),
        unit_source="HD_dataset" if "unit" in getattr(source, "attrs", {}) else "uncalibrated_HD_power",
        frequency_source_path="/" + schema.spectrograms_frequency_bin_table,
        frequency_bin_order="HD_export_order_signed_bins_kept_separate",
        cycle_boundary_workflow=timing_method,
        cycle_boundary_vessel_source=timing.cardiac_cycle_source,
        cycle_boundary_source_path="/" + processing_path(
            schema.cardiac_cycle.systolic_peak_frame_indices, VELOCITY_WORKFLOW_ROOTS[timing_method]),
        cycle_boundary_indices=timing.cycle_boundary_indexes, cycle_boundary_index_base=0,
        beat_duration_seconds=timing.cycle_durations_seconds,
        phase_endpoint="excluded", incomplete_edge_intervals="omitted",
        interpolation="spatial_linear_power_before_threshold_no_NaN_bridging_exact_anchors_preserved",
        threshold_reference="pixelwise_median_of_valid_beats_then_pool_phase_and_disc_pixels",
        threshold_method="skimage.filters.threshold_otsu", skimage_version=skimage.__version__,
        numpy_version=np.__version__, histogram_bins=result["histogram_bins"],
        histogram_range="finite_reference_min_max_per_frequency", histogram_weighting="equal_per_valid_phase_pixel",
        comparison="power_strictly_greater_than_threshold", power_preprocessing="none_linear_exported_power",
        negative_value_policy="preserve_finite_negative_power",
        nonfinite_value_policy="any_nonfinite_disc_pixel_invalidates_entire_beat_image_before_reference_median",
        aggregation="median_of_valid_individual_beat_fractions", denominator="fixed_geometric_disc_pixel_count",
        spatial_support="all_pixels_inside_image_centered_disc_including_vessels_and_optic_disc",
        disc_center_xy=[(source.shape[3] - 1) / 2, (source.shape[2] - 1) / 2],
        disc_radius_pixels=0.4 * source.shape[3], disc_boundary="distance_squared_le_radius_squared",
        spatial_coordinates="HD_export_yx_no_transform", working_memory_mb=result["working_memory_mb"],
        tile_pixels=result["tile_pixels"], max_source_passes_per_frequency=4,
        calculation_dtype="float64",
    )

    def payload(value, axes, unit="1", **extra):
        return value, {**attrs, "dimDesc": axes, "unit": unit, **extra}

    return {
        f"{root}/CVI": payload(result["cvi"], ["phase", "freq"]),
        f"{root}/PerBeat": payload(result["per_beat"], ["phase", "beat", "freq"]),
        f"{root}/Thresholds": payload(result["thresholds"], ["freq"], power_unit),
        f"{root}/Phase": payload(result["phase"], ["phase"]),
        f"{root}/Frequency": (frequencies[()], {**attrs, **dict(frequencies.attrs)}),
        f"{root}/DiscMask": payload(result["mask"], ["y", "x"]),
        f"{root}/PixelCount": payload(result["pixel_count"], [], "pixel"),
        f"{root}/ValidBeatCounts": payload(result["valid_counts"], ["phase", "freq"]),
        f"{root}/ThresholdFlags": payload(result["flags"], ["freq"], flag_meanings="0=valid,1=constant_reference,2=no_finite_reference"),
        f"{root}/ThresholdReferenceCounts": payload(result["reference_counts"], ["freq"]),
    }


def run_binary_occupancy_analysis(ctx) -> dict[str, object]:
    """Derive phase/frequency occupancy maps from the completed CVI analysis.

    These spatial products feed both movies and the frequency-averaged image.
    """
    root = EyeFlowOutputPaths.active().spectrograms_cvi_root
    cvi = ctx.output.h5.get(f"{root}/CVI")
    if cvi is None:
        return {}
    peaks = _spectral_peak_bounds(ctx.output, cvi.shape[1])
    if peaks is None:
        Logger.log_warning("Skipping phase occupancy map: valid stored index_max1/index_max2 are unavailable.")
    source = ctx.inputs.hd.as_holodoppler().optional_spectograms_dataset()
    with ctx.output.h5.build_arrays() as allocate:
        result = calculate_binary_occupancy_maps(
            source, cvi.attrs["cycle_boundary_indices"], ctx.output.h5.get(f"{root}/Thresholds")[()],
            np.isfinite(ctx.output.h5.get(f"{root}/PerBeat")[()]), peak_indices=peaks,
            working_memory_mb=float(cvi.attrs["working_memory_mb"]), array_allocator=allocate,
        )
    return _binary_occupancy_datasets(result, dict(cvi.attrs))


def _spectral_peak_bounds(output, bin_count):
    """Read the inclusive frequency-bin bounds selected by gap-spectrum metrics."""
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    peaks = []
    try:
        for name in ("index_max1", "index_max2"):
            node = output.h5.get(f"{root}/{name}")
            value = float(np.asarray(node[()]).reshape(())) if node is not None else np.nan
            if not np.isfinite(value) or value != int(value) or not 0 <= value < bin_count:
                return None
            peaks.append(int(value))
    except (TypeError, ValueError):
        return None
    return tuple(peaks) if peaks[0] <= peaks[1] else None


def _binary_occupancy_datasets(result, cvi_attrs):
    """Describe spatial occupancy maps, their reductions and valid support counts."""
    schema = EyeFlowOutputPaths.active()
    root = schema.spectrograms_cvi_root
    attrs = {**cvi_attrs, "unit": "1", "quantity": "mean_binary_occupancy",
             "definition": "mean of valid beat binary images, then equal-weight mean over bins or phase",
             "beat_aggregation": "arithmetic_mean_over_valid_beats_at_each_phase_and_frequency",
             "aggregation": "mean_valid_beats_then_mean_bins_or_phase",
             "threshold_source_path": f"/{root}/Thresholds",
             "validity_source_path": f"/{root}/PerBeat",
             "nonfinite_value_policy": "exclude_invalid_whole_beat_images; omit_missing_phase_frequency_cells; no_valid_samples_is_NaN",
             "post_aggregation_threshold": "none", "spatial_support": "DiscMask; outside_is_NaN",
             "display_range": [0., 1.], "playback_fps": 10.}
    values = {
        f"{root}/Movies/FrequencyMeanOverBins": (mean_frequency_occupancy(result["frequency_images"]), {
            **attrs, "dimDesc": ["y", "x"],
            "source_path": f"/{root}/Movies/FrequencyMean",
            "definition": "equal-weight mean over all finite signed-frequency movie maps",
            "aggregation_axes": ["freq"], "aggregation": "arithmetic_mean_over_finite_frequency_maps",
            "frequency_selection": "all_exported_signed_bins",
            "source_frequency_bin_count": result["frequency_images"].shape[0]}),
        f"{root}/Movies/FrequencyMean": (result["frequency_images"], {
            **attrs, "dimDesc": ["freq", "y", "x"], "aggregation_axes": ["beat", "phase"],
            "display_frame_order": "ascending_calibrated_signed_frequency; stored_bins_unchanged"}),
        f"{root}/Movies/ValidPhaseCounts": (result["frequency_phase_counts"], {
            "unit": "1", "dimDesc": ["freq"], "definition": "phases containing at least one valid beat image"}),
    }
    if result["phase_images"] is not None:
        values.update({
            f"{root}/Movies/PhaseMean": (result["phase_images"], {
                **attrs, "dimDesc": ["phase", "y", "x"], "aggregation_axes": ["beat", "freq"],
                "frequency_bin_indices": result["selected_bins"], "frequency_interval_endpoints": "both_included",
                "peak_index_source_paths": [f"/{schema.spectrograms_masked_metrics_root}/index_max1",
                                            f"/{schema.spectrograms_masked_metrics_root}/index_max2"]}),
            f"{root}/Movies/ValidFrequencyCounts": (result["phase_bin_counts"], {
                "unit": "1", "dimDesc": ["phase"],
                "definition": "selected bins containing at least one valid beat image"}),
        })
    return values


def disc_mask(height, width):
    y, x = np.ogrid[:height, :width]
    return (x - (width - 1) / 2)**2 + (y - (height - 1) / 2)**2 <= (0.4 * width)**2


def calculate_cvi(source, boundaries, frequencies, *, sample_count=DEFAULT_PHASE_SAMPLE_COUNT,
                  histogram_bins=256, working_memory_mb=64):
    """Return source-order fractions, thresholds, support and validity counts.

    Working arrays scale with phase * beats * tile pixels, rather than image
    size or frequency count. Only reduced scientific outputs are retained.
    """
    if source.ndim != 4 or any(n == 0 for n in source.shape):
        raise ValueError("CVI requires nonempty (time, freq, y, x) linear power.")
    if np.dtype(source.dtype).kind not in "iuf":
        raise TypeError("CVI requires real numeric power.")
    if isinstance(sample_count, bool) or not isinstance(sample_count, int) or sample_count < 2:
        raise ValueError("sample_count must be an integer of at least two.")
    if isinstance(histogram_bins, bool) or not isinstance(histogram_bins, int) or histogram_bins < 2:
        raise ValueError("histogram_bins must be an integer of at least two.")
    if isinstance(working_memory_mb, bool) or not isinstance(working_memory_mb, Real) or not np.isfinite(working_memory_mb) or working_memory_mb <= 0:
        raise ValueError("working_memory_mb must be positive and finite.")
    boundaries = normalize_cycle_boundaries(boundaries, source.shape[0], index_base=0)
    nt, nf, height, width = source.shape
    frequencies = np.asarray(frequencies, dtype=float)
    if frequencies.shape != (nf,) or not np.all(np.isfinite(frequencies)):
        raise ValueError("CVI requires one calibrated signed frequency per source bin.")
    nb = boundaries.size - 1
    support = disc_mask(height, width)
    denominator = int(support.sum())
    # Conservative allowance for source, interpolation, reference and reduction copies.
    tile_pixels = max(1, int(working_memory_mb * 1024**2 // (64 * (nt + sample_count * nb))))
    tile_pixels = min(tile_pixels, width)
    phase = np.arange(sample_count) / sample_count
    per_beat = np.full((sample_count, nb, nf), np.nan)
    valid_counts = np.zeros((sample_count, nf), dtype=np.int32)
    thresholds = np.full(nf, np.nan)
    # 0 valid, 1 constant reference, 2 no finite reference (includes empty disc).
    flags = np.full(nf, 2, dtype=np.uint8)
    reference_counts = np.zeros(nf, dtype=np.int64)

    def tiles(f):
        for y in range(height):
            for x in range(0, width, tile_pixels):
                stop = min(x + tile_pixels, width)
                mask = support[y, x:stop]
                if not mask.any():
                    continue
                values = source[:, f, y:y + 1, x:stop]
                yield y, x, stop, mask, align_power(values.reshape(nt, -1), boundaries, sample_count)

    def reference_chunks(f, valid):
        for _, _, _, mask, aligned in tiles(f):
            if not mask.any():
                continue
            values = aligned[:, :, mask].copy()
            values[~valid] = np.nan
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                reference = np.nanmedian(values, axis=1)
            yield reference[np.isfinite(reference)]

    for f in range(nf):
        valid = np.ones((sample_count, nb), dtype=bool)
        for _, _, _, mask, aligned in tiles(f):
            valid &= np.all(np.isfinite(aligned[:, :, mask]), axis=2)
        if not denominator:
            valid[:] = False
        valid_counts[:, f] = valid.sum(axis=1)
        low, high = np.inf, -np.inf
        for reference in reference_chunks(f, valid):
            if reference.size:
                low = min(low, float(reference.min()))
                high = max(high, float(reference.max()))
                reference_counts[f] += reference.size
        if reference_counts[f]:
            flags[f] = 1 if low == high else 0
            if low == high:
                thresholds[f] = low
            else:
                counts = np.zeros(histogram_bins, dtype=np.int64)
                edges = np.linspace(low, high, histogram_bins + 1)
                for reference in reference_chunks(f, valid):
                    counts += np.histogram(reference, bins=edges)[0]
                centers = edges[:-1] + np.diff(edges) / 2
                thresholds[f] = threshold_otsu(hist=(counts, centers))
        white = np.zeros((sample_count, nb), dtype=np.int64)
        for y, x, stop, mask, aligned in tiles(f):
            white += np.count_nonzero(aligned[:, :, mask] > thresholds[f], axis=2)
        if denominator and np.isfinite(thresholds[f]):
            fractions = white / denominator
            fractions[~valid] = np.nan
            per_beat[:, :, f] = fractions
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        cvi = np.nanmedian(per_beat, axis=1)
    return dict(cvi=cvi, per_beat=per_beat, thresholds=thresholds, flags=flags,
                valid_counts=valid_counts, reference_counts=reference_counts,
                mask=support, pixel_count=denominator, phase=phase, frequencies=frequencies,
                tile_pixels=tile_pixels,
                histogram_bins=histogram_bins, working_memory_mb=working_memory_mb)


def calculate_binary_occupancy_maps(source, boundaries, thresholds, valid_images, *,
                                    peak_indices=None, working_memory_mb=64, array_allocator=None):
    """Average binary images over valid beats, then equally over bins or phase.

    valid_images has (phase, beat, freq) axes and comes from finite PerBeat CVI.
    Whole-image invalidity is therefore identical to the existing CVI contract.
    No threshold is applied after averaging. Returned maps are fractions [0,1],
    with NaN outside the disc and where no valid contributions remain.
    """
    if source.ndim != 4 or any(n == 0 for n in source.shape):
        raise ValueError("Binary occupancy maps require nonempty (time, freq, y, x) power.")
    nt, nf, height, width = source.shape
    boundaries = normalize_cycle_boundaries(boundaries, nt, index_base=0)
    thresholds = np.asarray(thresholds, dtype=float)
    valid = np.asarray(valid_images)
    if (thresholds.shape != (nf,) or valid.ndim != 3 or valid.dtype != np.bool_
            or valid.shape[1:] != (len(boundaries) - 1, nf) or valid.shape[0] < 2):
        raise ValueError("Thresholds and valid-image axes must match phase, complete beats and source bins.")
    if not np.isfinite(working_memory_mb) or working_memory_mb <= 0:
        raise ValueError("working_memory_mb must be positive and finite.")
    phase_count, beat_count, _ = valid.shape
    valid = valid & np.isfinite(thresholds)[None, None, :]
    valid_beat_counts = valid.sum(axis=1)
    valid_phase = valid_beat_counts > 0
    frequency_phase_counts = valid_phase.sum(axis=0)
    selected = np.empty(0, dtype=int)
    if peak_indices is not None:
        pair = np.asarray(peak_indices)
        if (pair.shape != (2,) or pair.dtype.kind not in "iu"
                or not 0 <= pair[0] <= pair[1] < nf):
            raise ValueError("Occupancy-map peak indices must be an inclusive pair of source bin indices.")
        selected = np.arange(pair[0], pair[1] + 1)
    phase_bin_counts = valid_phase[:, selected].sum(axis=1)
    support = disc_mask(height, width)
    if not support.any():
        valid[:] = False
        valid_beat_counts[:] = 0
        valid_phase[:] = False
        frequency_phase_counts[:] = 0
        phase_bin_counts[:] = 0
    allocate = array_allocator or (lambda shape: np.full(shape, np.nan, dtype=np.float32))
    phase_images = allocate((phase_count, height, width)) if selected.size else None
    frequency_images = allocate((nf, height, width))
    # One owned source slab plus interpolants, binary sums and map tiles.
    bytes_per_pixel = nt * nf * np.dtype(source.dtype).itemsize + 64 * phase_count * beat_count
    pixels = max(1, int(working_memory_mb * 2**20 // bytes_per_pixel))
    columns = min(width, pixels)
    rows = max(1, pixels // width) if columns == width else 1
    for y in range(0, height, rows):
        stop_y = min(y + rows, height)
        for x in range(0, width, columns):
            stop_x = min(x + columns, width)
            mask = support[y:stop_y, x:stop_x].reshape(-1)
            if not mask.any():
                continue
            block = source[:, :, y:stop_y, x:stop_x]
            phase_sum = np.zeros((phase_count, mask.size), dtype=np.float64)
            frequency_tile = np.full((nf, mask.size), np.nan, dtype=np.float32)
            for f in range(nf):
                if not valid_phase[:, f].any():
                    continue
                aligned = align_power(block[:, f].reshape(nt, -1), boundaries, phase_count)
                binary = (aligned > thresholds[f]) & valid[:, :, f, None]
                occupancy = np.zeros((phase_count, mask.size), dtype=np.float64)
                np.divide(binary.sum(axis=1), valid_beat_counts[:, f, None],
                          out=occupancy, where=valid_phase[:, f, None])
                frequency_tile[f] = occupancy.sum(axis=0) / frequency_phase_counts[f]
                if selected.size and selected[0] <= f <= selected[-1]:
                    phase_sum += occupancy
            frequency_tile[:, ~mask] = np.nan
            frequency_images[:, y:stop_y, x:stop_x] = frequency_tile.reshape(nf, stop_y - y, stop_x - x)
            if phase_images is not None:
                phase_tile = np.full(phase_sum.shape, np.nan)
                np.divide(phase_sum, phase_bin_counts[:, None], out=phase_tile,
                          where=phase_bin_counts[:, None] > 0)
                phase_tile[:, ~mask] = np.nan
                phase_images[:, y:stop_y, x:stop_x] = phase_tile.reshape(phase_count, stop_y - y, stop_x - x)
            del block
    return dict(phase_images=phase_images, frequency_images=frequency_images,
                phase_bin_counts=phase_bin_counts, frequency_phase_counts=frequency_phase_counts,
                selected_bins=selected)


def mean_frequency_occupancy(images):
    """Average all signed-frequency movie maps, reading one frame at a time."""
    total = np.zeros(images.shape[1:], dtype=np.float64)
    count = np.zeros(images.shape[1:], dtype=np.int32)
    for frequency in range(images.shape[0]):
        frame = np.asarray(images[frequency], dtype=np.float64)
        valid = np.isfinite(frame)
        total += np.where(valid, frame, 0.)
        count += valid
    mean = np.full(total.shape, np.nan)
    np.divide(total, count, out=mean, where=count > 0)
    return mean.astype(np.float32)
