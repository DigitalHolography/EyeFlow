"""High-level orchestration for retinal displacement-map calculation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

from .constants import (
    CURVATURE_CONSTRAINT_WEIGHT,
    CURVATURE_TIME_STEP,
    PDE_REGISTRATION_METHODS,
)
from .filtering import CenteredMedianBuffer
from .outputs import OutputCaches, select_display_range, write_magnitude_video
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
from .sources import FrameSequence, compute_mean_reference, load_binary_mask


def create_retinal_motion_map(
    args: MotionMapConfig,
    *,
    analysis_mask_array: np.ndarray | None = None,
    magnitude_video_path: Path | None = None,
) -> dict[str, Path]:
    if cv2 is None:
        raise RuntimeError("OpenCV is required to create a displacement map.")
    if tqdm is None:
        raise RuntimeError("tqdm is required to create a displacement map.")
    require_registration_backend(args.registration_method)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequence = FrameSequence(
        args.input,
        args.h5_dataset,
        args.h5_frame_axis,
        args.h5_fps,
        args.h5_low_percentile,
        args.h5_high_percentile,
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
    cv2.imwrite(str(args.output_dir / "analysis_mask.png"), analysis_mask)

    print("Calcul de la référence moyenne...")
    reference_raw, actual_reference_count = compute_mean_reference(sequence, args.max_frames)
    frame_capacity = min(frame_capacity, actual_reference_count)
    confidence = photometric_confidence(reference_raw, valid_analysis, photometric_config)
    registration_weight = analysis_soft * confidence
    reference = structural_registration_signal(
        reference_raw, photometric_config, valid_analysis, args.structure_edge_weight
    )

    cv2.imwrite(
        str(args.output_dir / "reference_mean.png"),
        np.round(np.clip(reference_raw, 0.0, 1.0) * 255.0).astype(np.uint8),
    )
    cv2.imwrite(
        str(args.output_dir / "reference_registration.png"),
        np.round(np.clip(reference, 0.0, 1.0) * 255.0).astype(np.uint8),
    )
    cv2.imwrite(
        str(args.output_dir / "registration_confidence.png"),
        np.round(np.clip(confidence, 0.0, 1.0) * 255.0).astype(np.uint8),
    )

    caches = OutputCaches(
        args.output_dir,
        frame_capacity,
        sequence.height,
        sequence.width,
        args.save_field,
        analysis_mask > 0,
    )

    previous_registration_field = np.zeros((sequence.height, sequence.width, 2), dtype=np.float32)
    previous_smoothed_field: np.ndarray | None = None
    median_buffer = CenteredMedianBuffer(args.temporal_median_window)
    raw_frame_count = 0

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

    progress = tqdm(total=frame_capacity, desc="Registration magnitude", unit="frame")
    try:
        for frame in sequence.iter_frames(frame_capacity):
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
                bspline_mesh_size=args.bspline_mesh_size,
            )
            raw_field[analysis_mask == 0] = 0.0
            previous_registration_field = raw_field.copy()
            raw_frame_count += 1

            for median_field in median_buffer.push(raw_field):
                append_temporally_filtered(median_field)

            progress.update(1)

        for median_field in median_buffer.finish():
            append_temporally_filtered(median_field)
    finally:
        progress.close()
        caches.flush()

    display_minimum, display_maximum = select_display_range(
        caches,
        args.normalization,
        args.low_percentile,
        args.high_percentile,
        args.fixed_max_px,
    )
    output_video = magnitude_video_path or (args.output_dir / "displacement_magnitude.mp4")
    output_video.parent.mkdir(parents=True, exist_ok=True)
    write_magnitude_video(
        caches,
        output_video,
        args.codec,
        sequence.fps,
        (sequence.width, sequence.height),
        display_minimum,
        display_maximum,
        args.gamma,
        args.visualization_sigma,
        analysis_mask > 0,
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
        "bspline_mesh_size": int(args.bspline_mesh_size),
        "curvature_time_step": float(CURVATURE_TIME_STEP),
        "curvature_constraint_weight": float(CURVATURE_CONSTRAINT_WEIGHT),
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
        "display_minimum_px": float(display_minimum),
        "display_maximum_px": float(display_maximum),
        "gamma": float(args.gamma),
        "visualization_sigma": float(args.visualization_sigma),
        "saved_dense_field": bool(args.save_field),
        "stabilization_modified": False,
    }
    (args.output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    caches.close()

    print(f"Vidéo créée : {output_video}")
    print(f"Magnitude numérique créée : {args.output_dir / 'displacement_magnitude.npy'}")
    if args.save_field:
        print(f"Champ dense créé : {args.output_dir / 'displacement_field.npy'}")
    print("La stabilisation d'entrée n'a pas été recalculée ni modifiée.")
    outputs = {
        "output_dir": args.output_dir,
        "magnitude_video": output_video,
        "magnitude": args.output_dir / "displacement_magnitude.npy",
        "metadata": args.output_dir / "metadata.json",
    }
    if args.save_field:
        outputs["displacement_field"] = args.output_dir / "displacement_field.npy"
    return outputs
