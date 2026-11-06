"""High-level orchestration for retinal displacement-map calculation."""

from __future__ import annotations

from contextlib import nullcontext
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import h5py
import numpy as np

from calculations.math.temporal_median import CenteredMedianBuffer
from calculations.temporal_reference import compute_mean_reference
from input_output.writers.displacement_arrays import OutputCaches
from input_output.readers.frame_sequences import FrameReplay, FrameSequence, load_binary_mask
from input_output.writers.magnitude_avi import write_magnitude_avi
from input_output.writers.sidecars import write_json_file
from input_output.writers.png import write_cv_png_file
from utils.logger import Logger
from runtime_limits import cap_parallel_jobs

try:
    import cv2
except ImportError:
    cv2 = None

from .constants import PDE_REGISTRATION_METHODS
from .parameters import MotionMapConfig, PhotometricConfig
from .preprocessing import (
    photometric_confidence,
    soft_mask,
    structural_registration_signal,
)
from .registration import (
    estimate_registration_field,
    registration_method_class_name,
    require_registration_backend,
)


def create_retinal_motion_map(
    args: MotionMapConfig,
    *,
    analysis_mask_array: np.ndarray | None = None,
    magnitude_video_path: Path | None = None,
    h5_source: h5py.Dataset | None = None,
    frame_replay: FrameReplay | None = None,
    raw_reference: tuple[np.ndarray, int] | None = None,
    render_video: bool = True,
    registration_workers: int = 1,
    native_threads: int | None = None,
    log_label: str = "Displacement",
) -> dict[str, Path]:
    if cv2 is None:
        raise RuntimeError("OpenCV is required to create a displacement map.")
    require_registration_backend(args.registration_method)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequence = frame_replay.source if frame_replay is not None else FrameSequence(
        args.input,
        args.h5_dataset,
        args.h5_frame_axis,
        args.h5_fps,
        args.h5_low_percentile,
        args.h5_high_percentile,
        h5_source=h5_source,
        defer_h5_range=True,
    )
    frame_capacity = sequence.frame_count
    if args.max_frames is not None:
        frame_capacity = min(frame_capacity, args.max_frames)

    if analysis_mask_array is not None:
        supplied_mask = np.squeeze(np.asarray(analysis_mask_array))
        expected_shape = (sequence.height, sequence.width)
        if supplied_mask.shape != expected_shape:
            raise ValueError(
                "analysis mask shape does not match the moment frames: "
                f"{supplied_mask.shape} != {expected_shape}"
            )
        analysis_mask = (supplied_mask != 0).astype(np.uint8) * 255
        if args.invert_analysis_mask:
            analysis_mask = cv2.bitwise_not(analysis_mask)
        analysis_mask = cv2.dilate(analysis_mask, np.ones((3, 3), np.uint8), iterations=20)
    elif args.analysis_mask is not None:
        analysis_mask = load_binary_mask(
            args.analysis_mask,
            (sequence.height, sequence.width),
            args.invert_analysis_mask,
        )
    else:
        if args.invert_analysis_mask:
            raise ValueError("--invert-analysis-mask nécessite --analysis-mask")
        analysis_mask = np.ones((sequence.height, sequence.width), np.uint8) * 255
    analysis_soft = soft_mask(analysis_mask, args.mask_feather_sigma)
    valid_analysis = analysis_mask > 0
    photometric_config = PhotometricConfig(
        mode=args.photometric_normalization,
        illumination_sigma=float(args.illumination_sigma),
        local_contrast_sigma=float(args.local_contrast_sigma),
        clip=float(args.photometric_clip),
        confidence_floor=float(args.photometric_confidence_floor),
        dark_percentile=float(args.photometric_dark_percentile),
    )
    write_cv_png_file(args.output_dir / "analysis_mask.png", analysis_mask)

    frame_context = (
        nullcontext(frame_replay) if frame_replay is not None
        else sequence.cached_for_replay(args.max_frames, args.output_dir)
    )
    with frame_context as replay:
        if raw_reference is None:
            reference_raw, actual_reference_count = compute_mean_reference(
                replay, args.max_frames,
                on_progress=_FrameProgress(f"{log_label}: mean reference", frame_capacity),
            )
        else:
            reference_raw, actual_reference_count = raw_reference
            if reference_raw.shape != (sequence.height, sequence.width) or not 0 < actual_reference_count <= frame_capacity:
                raise ValueError("Shared reference does not match the replay shape/frame limit.")
        frame_capacity = min(frame_capacity, actual_reference_count)
        confidence = photometric_confidence(reference_raw, valid_analysis, photometric_config)
        registration_weight = analysis_soft * confidence
        reference = structural_registration_signal(
            reference_raw, photometric_config, valid_analysis, args.structure_edge_weight
        )

        write_cv_png_file(
            args.output_dir / "reference_mean.png",
            np.round(np.clip(reference_raw, 0.0, 1.0) * 255.0).astype(np.uint8),
        )
        write_cv_png_file(
            args.output_dir / "reference_registration.png",
            np.round(np.clip(reference, 0.0, 1.0) * 255.0).astype(np.uint8),
        )
        write_cv_png_file(
            args.output_dir / "registration_confidence.png",
            np.round(np.clip(confidence, 0.0, 1.0) * 255.0).astype(np.uint8),
        )

        caches = OutputCaches(
            args.output_dir,
            frame_capacity,
            sequence.height,
            sequence.width,
            args.save_field,
        )

        previous_registration_field = np.zeros((sequence.height, sequence.width, 2), dtype=np.float32)
        previous_smoothed_field: np.ndarray | None = None
        median_buffer = CenteredMedianBuffer(args.temporal_median_window)

        def append_temporally_filtered(raw_field: np.ndarray) -> None:

            nonlocal previous_smoothed_field
            field = raw_field.astype(np.float32, copy=True)

            if previous_smoothed_field is not None and args.temporal_alpha < 1.0:
                alpha = np.float32(args.temporal_alpha)
                field = (alpha * field + (np.float32(1.0) - alpha) * previous_smoothed_field).astype(
                    np.float32
                )

            field[analysis_mask == 0] = 0.0
            previous_smoothed_field = field.copy()
            magnitude = np.hypot(field[..., 0], field[..., 1]).astype(np.float32)
            magnitude[analysis_mask == 0] = 0.0
            caches.append(field, magnitude)

        workers = cap_parallel_jobs(max(1, int(registration_workers)))
        # Only this measured backend with zero initialization has independent
        # raw jobs. Previous-field initialization always retains its ordering.
        if args.registration_initialization != "zero" or args.registration_method != "level_set_motion":
            workers = 1

        def register_frame(indexed_frame):
            frame_index, frame = indexed_frame
            gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY).astype(np.float32) / 255.0
            moving = structural_registration_signal(
                gray, photometric_config, valid_analysis, args.structure_edge_weight
            )
            initial_field = (
                previous_registration_field
                if args.registration_initialization == "previous"
                else None
            )
            raw_field, _registration_metric = estimate_registration_field(
                fixed_full=reference,
                moving_full=moving,
                mask_soft_full=registration_weight,
                initial_full=initial_field,
                method=args.registration_method,
                scale=args.scale,
                iterations=args.iterations,
                field_sigma=args.field_sigma,
                update_sigma=args.update_sigma,
                metric_radius=args.registration_metric_radius,
                learning_rate=args.registration_learning_rate,
                native_threads=native_threads,
            )
            return frame_index, raw_field

        def registered_frames():
            frames = iter(enumerate(replay.iter_frames(frame_capacity), start=1))
            if workers == 1:
                for frame in frames:
                    yield register_frame(frame)
                return
            # The owning thread reads replay frames, consumes results in input
            # order and applies all temporal filtering/output writes. Jobs own
            # separate frames and SimpleITK filters, with bounded in-flight data.
            with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="displacement-frame") as executor:
                pending = deque()
                for _ in range(workers):
                    frame = next(frames, None)
                    if frame is not None:
                        pending.append(executor.submit(register_frame, frame))
                while pending:
                    yield pending.popleft().result()
                    frame = next(frames, None)
                    if frame is not None:
                        pending.append(executor.submit(register_frame, frame))

        progress = _FrameProgress(f"{log_label}: registration", frame_capacity)
        try:
            for frame_index, raw_field in registered_frames():
                raw_field[analysis_mask == 0] = 0.0
                previous_registration_field = raw_field.copy()

                for median_field in median_buffer.push(raw_field):
                    append_temporally_filtered(median_field)

                progress(frame_index, frame_capacity)

            for median_field in median_buffer.finish():
                append_temporally_filtered(median_field)
        except BaseException:
            caches.close()
            raise
        finally:
            caches.flush()

    try:
        output_video = magnitude_video_path or (args.output_dir / "displacement_magnitude.avi")
        output_video.parent.mkdir(parents=True, exist_ok=True)
        display_minimum, display_maximum = (None, None)
        if render_video:
            display_minimum, display_maximum = write_magnitude_avi(
                output_video, caches.magnitude, frame_count=caches.count,
                valid_mask=analysis_mask > 0, fps=sequence.fps,
                normalization=args.normalization,
                low_percentile=args.low_percentile, high_percentile=args.high_percentile,
                fixed_maximum=args.fixed_max_px, gamma=args.gamma,
                visualization_sigma=args.visualization_sigma,
                on_progress=_FrameProgress(f"{log_label}: magnitude video", caches.count),
            )

        metadata = {
            "input": str(args.input.resolve()),
            "input_kind": "hdf5" if sequence.is_h5 else "video",
            "h5_dataset": sequence.h5_dataset if sequence.is_h5 else None,
            "h5_frame_axis": sequence.h5_frame_axis if sequence.is_h5 else None,
            "h5_display_range": (
                [sequence.h5_value_min, sequence.h5_value_max] if sequence.is_h5 else None
            ),
            "recommended_video_input": "global_stabilized.mp4",
            "frames": int(caches.count),
            "width": int(sequence.width),
            "height": int(sequence.height),
            "fps": float(sequence.fps),
            "reference": "temporal mean of the input sequence",
            "reference_frame_count": int(actual_reference_count),
            "analysis_mask": "analysis_mask.png",
            "registration_method": args.registration_method,
            "algorithm": registration_method_class_name(args.registration_method),
            "registration_initialization": args.registration_initialization,
            "independent_zero_initial_field_per_frame": (args.registration_initialization == "zero"),
            "previous_raw_field_initialization": (args.registration_initialization == "previous"),
            "iterations": int(args.iterations),
            "calculation_scale": float(args.scale),
            "field_sigma": float(args.field_sigma),
            "update_sigma": float(args.update_sigma),
            "registration_metric": (
                "native PDE metric"
                if args.registration_method in PDE_REGISTRATION_METHODS
                else "ANTSNeighborhoodCorrelation"
            ),
            "registration_metric_radius": int(args.registration_metric_radius),
            "registration_learning_rate": float(args.registration_learning_rate),
            "temporal_median_window": int(args.temporal_median_window),
            "temporal_median": ("centered component-wise median on dx and dy with replicated edges"),
            "temporal_ema": "optional causal EMA applied after the centered median",
            "temporal_alpha": float(args.temporal_alpha),
            "temporal_smoothing_enabled": bool(
                args.temporal_median_window > 1 or args.temporal_alpha < 1.0
            ),
            "photometric_normalization": args.photometric_normalization,
            "illumination_sigma": float(args.illumination_sigma),
            "local_contrast_sigma": float(args.local_contrast_sigma),
            "photometric_clip": float(args.photometric_clip),
            "photometric_confidence_floor": float(args.photometric_confidence_floor),
            "photometric_dark_percentile": float(args.photometric_dark_percentile),
            "structure_edge_weight": float(args.structure_edge_weight),
            "registration_confidence": "registration_confidence.png",
            "registration_signal": "reference_registration.png",
            "magnitude_formula": "sqrt(dx**2 + dy**2)",
            "magnitude_file": "displacement_magnitude.npy",
            "magnitude_shape": [int(caches.count), int(sequence.height), int(sequence.width)],
            "magnitude_units": "pixels",
            "magnitude_is_normalized": False,
            "magnitude_reference": (
                "displacement from each input frame to the temporal mean reference"
            ),
            "dense_field_file": "displacement_field.npy" if args.save_field else None,
            "dense_field_shape": (
                [int(caches.count), int(sequence.height), int(sequence.width), 2]
                if args.save_field
                else None
            ),
            "dense_field_components": ["dx", "dy"] if args.save_field else None,
            "normalization": args.normalization,
            "display_range_source": "temporal_median_magnitude_image",
            "display_minimum_px": display_minimum,
            "display_maximum_px": display_maximum,
            "gamma": float(args.gamma),
            "visualization_sigma": float(args.visualization_sigma),
            "saved_dense_field": bool(args.save_field),
            "stabilization_modified": False,
        }
        write_json_file(args.output_dir / "metadata.json", metadata)

        outputs = {
            "output_dir": args.output_dir,
            "magnitude_video": output_video,
            "magnitude": args.output_dir / "displacement_magnitude.npy",
            "metadata": args.output_dir / "metadata.json",
        }
        if args.save_field:
            outputs["displacement_field"] = args.output_dir / "displacement_field.npy"
        return outputs
    finally:
        caches.close()


def render_retinal_motion_map_video(
    args: MotionMapConfig, *, magnitude_video_path: Path, log_label: str = "Displacement",
) -> Path:
    """Encode closed numerical backing; no input HDF5 readers or writers are used."""
    metadata_path = args.output_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    magnitude = np.load(args.output_dir / "displacement_magnitude.npy", mmap_mode="r")
    try:
        mask = cv2.imread(str(args.output_dir / "analysis_mask.png"), cv2.IMREAD_GRAYSCALE)
        if mask is None:
            raise ValueError("Displacement artifact mask is unavailable.")
        minimum, maximum = write_magnitude_avi(
            magnitude_video_path, magnitude, frame_count=metadata["frames"],
            valid_mask=mask > 0, fps=metadata["fps"], normalization=args.normalization,
            low_percentile=args.low_percentile, high_percentile=args.high_percentile,
            fixed_maximum=args.fixed_max_px, gamma=args.gamma,
            visualization_sigma=args.visualization_sigma,
            on_progress=_FrameProgress(f"{log_label}: magnitude video", metadata["frames"]),
        )
    finally:
        magnitude._mmap.close()
    metadata.update(display_minimum_px=float(minimum), display_maximum_px=float(maximum))
    write_json_file(metadata_path, metadata)
    return magnitude_video_path


class _FrameProgress:
    """Emit a small number of app-log updates for frame-based work."""

    def __init__(self, label: str, total: int, *, interval_percent: int = 20):
        self.label = str(label)
        self.total = max(0, int(total))
        self.interval_percent = int(interval_percent)
        self._next_percent = self.interval_percent
        Logger.log_info(f"{self.label} started ({self.total} frames).")

    def __call__(self, completed: int, total: int) -> None:
        total = max(0, int(total))
        completed = max(0, min(int(completed), total))
        if total == 0:
            return
        percent = min(100, int(completed * 100 / total))
        if percent < self._next_percent and completed < total:
            return
        Logger.log_info(
            f"{self.label}: {percent}% ({completed}/{total} frames)."
        )
        self._next_percent = min(
            101,
            (percent // self.interval_percent + 1) * self.interval_percent,
        )
