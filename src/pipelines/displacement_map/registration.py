"""SimpleITK and ITK displacement-registration backends."""

from __future__ import annotations

from typing import Any

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None

try:
    import SimpleITK as sitk
except ImportError:
    sitk = None

try:
    import itk
except ImportError:
    itk = None

from .constants import (
    CURVATURE_CONSTRAINT_WEIGHT,
    CURVATURE_TIME_STEP,
    ITK_REGISTRATION_METHODS,
    PDE_REGISTRATION_METHODS,
    RegistrationMethod,
)


def require_registration_backend(method: RegistrationMethod) -> None:
    """Raise a useful error when the selected backend is unavailable."""

    if cv2 is None:
        raise RuntimeError("OpenCV is required: install the EyeFlow 'pipelines' dependency extra.")
    if sitk is None:
        raise RuntimeError(
            "SimpleITK is required: install the EyeFlow 'pipelines' dependency extra."
        )
    if method in ITK_REGISTRATION_METHODS:
        _require_itk_method(method)


def create_pde_registration_filter(method: RegistrationMethod) -> Any:
    if sitk is None:
        raise RuntimeError("SimpleITK is not installed.")
    filters = {
        "classic_demons": sitk.DemonsRegistrationFilter,
        "symmetric_forces_demons": sitk.SymmetricForcesDemonsRegistrationFilter,
        "fast_symmetric_demons": sitk.FastSymmetricForcesDemonsRegistrationFilter,
        "diffeomorphic_demons": sitk.DiffeomorphicDemonsRegistrationFilter,
        "level_set_motion": sitk.LevelSetMotionRegistrationFilter,
    }
    try:
        return filters[method]()
    except KeyError as exc:
        raise ValueError(f"Unsupported PDE registration method: {method}") from exc


def registration_method_class_name(method: RegistrationMethod) -> str:
    return {
        "classic_demons": "SimpleITK DemonsRegistrationFilter",
        "symmetric_forces_demons": ("SimpleITK SymmetricForcesDemonsRegistrationFilter"),
        "fast_symmetric_demons": ("SimpleITK FastSymmetricForcesDemonsRegistrationFilter"),
        "diffeomorphic_demons": ("SimpleITK DiffeomorphicDemonsRegistrationFilter"),
        "level_set_motion": "SimpleITK LevelSetMotionRegistrationFilter",
        "displacement_field": ("SimpleITK ImageRegistrationMethod + DisplacementFieldTransform"),
        "bspline": "SimpleITK ImageRegistrationMethod + BSplineTransform",
        "affine": "SimpleITK ImageRegistrationMethod + AffineTransform",
        "curvature_registration": "ITK CurvatureRegistrationFilter",
        "gpu_demons": "ITK GPUDemonsRegistrationFilter (GPU/OpenCL)",
    }[method]


def make_sitk_scalar(array: np.ndarray) -> Any:
    _require_simpleitk()
    image = sitk.GetImageFromArray(np.asarray(array, dtype=np.float32), isVector=False)
    image.SetSpacing((1.0, 1.0))
    return image


def make_sitk_vector(array: np.ndarray) -> Any:
    _require_simpleitk()
    field = np.asarray(array, dtype=np.float64)
    if field.ndim != 3 or field.shape[-1] != 2:
        raise ValueError(f"Displacement field must have shape (y, x, 2), got {field.shape}.")
    image = sitk.GetImageFromArray(field, isVector=True)
    image.SetSpacing((1.0, 1.0))
    return image


def resize_scalar(
    array: np.ndarray,
    shape: tuple[int, int],
    *,
    interpolation: int | None = None,
) -> np.ndarray:
    _require_opencv()
    source = np.asarray(array, dtype=np.float32)
    target_height, target_width = (int(size) for size in shape)
    if source.shape == (target_height, target_width):
        return source.astype(np.float32, copy=True)
    if interpolation is None:
        interpolation = (
            cv2.INTER_AREA
            if target_height < source.shape[0] or target_width < source.shape[1]
            else cv2.INTER_LINEAR
        )
    return cv2.resize(
        source,
        (target_width, target_height),
        interpolation=interpolation,
    ).astype(np.float32)


def resize_displacement_to_full(
    field: np.ndarray,
    shape: tuple[int, int],
) -> np.ndarray:
    """Resize a field and preserve its displacement units in pixels."""

    _require_opencv()
    source = np.asarray(field, dtype=np.float32)
    if source.ndim != 3 or source.shape[-1] != 2:
        raise ValueError(f"Displacement field must have shape (y, x, 2), got {source.shape}.")
    target_height, target_width = (int(size) for size in shape)
    source_height, source_width = source.shape[:2]
    if (source_height, source_width) == (target_height, target_width):
        return source.astype(np.float32, copy=True)

    resized = cv2.resize(
        source,
        (target_width, target_height),
        interpolation=cv2.INTER_LINEAR,
    ).astype(np.float32)
    resized[..., 0] *= np.float32(target_width / source_width)
    resized[..., 1] *= np.float32(target_height / source_height)
    return resized


def transform_to_displacement_array(transform: Any, reference: Any) -> np.ndarray:
    _require_simpleitk()
    field = sitk.TransformToDisplacementField(
        transform,
        sitk.sitkVectorFloat64,
        reference.GetSize(),
        reference.GetOrigin(),
        reference.GetSpacing(),
        reference.GetDirection(),
    )
    return np.asarray(sitk.GetArrayFromImage(field), dtype=np.float32)


def make_displacement_transform(field: np.ndarray) -> Any:
    _require_simpleitk()
    return sitk.DisplacementFieldTransform(make_sitk_vector(field))


def run_pde_registration(
    fixed: Any,
    moving: Any,
    initial_field: np.ndarray,
    method: RegistrationMethod,
    iterations: int,
    field_sigma: float,
) -> tuple[np.ndarray, float]:
    registration = create_pde_registration_filter(method)
    registration.SetNumberOfIterations(int(iterations))
    if hasattr(registration, "SetStandardDeviations"):
        registration.SetStandardDeviations(float(field_sigma))
    result = registration.Execute(fixed, moving, make_sitk_vector(initial_field))
    metric = _numeric_filter_metric(registration)
    return np.asarray(sitk.GetArrayFromImage(result), dtype=np.float32), metric


def run_image_registration_method(
    fixed: Any,
    moving: Any,
    initial_field: np.ndarray,
    method: RegistrationMethod,
    iterations: int,
    field_sigma: float,
    update_sigma: float,
    metric_radius: int,
    learning_rate: float,
    bspline_mesh_size: int,
) -> tuple[np.ndarray, float]:
    _require_simpleitk()
    registration = sitk.ImageRegistrationMethod()
    registration.SetMetricAsANTSNeighborhoodCorrelation(int(metric_radius))
    registration.SetInterpolator(sitk.sitkLinear)
    registration.SetOptimizerAsGradientDescent(
        learningRate=float(learning_rate),
        numberOfIterations=int(iterations),
        convergenceMinimumValue=1.0e-6,
        convergenceWindowSize=max(1, min(10, int(iterations))),
    )
    registration.SetOptimizerScalesFromPhysicalShift()

    if method == "displacement_field":
        transform = make_displacement_transform(initial_field)
        transform.SetSmoothingGaussianOnUpdate(
            varianceForUpdateField=float(update_sigma) ** 2,
            varianceForTotalField=float(field_sigma) ** 2,
        )
    elif method == "bspline":
        transform = sitk.BSplineTransformInitializer(
            fixed,
            [int(bspline_mesh_size)] * 2,
            order=3,
        )
    elif method == "affine":
        transform = sitk.CenteredTransformInitializer(
            fixed,
            moving,
            sitk.AffineTransform(2),
            sitk.CenteredTransformInitializerFilter.GEOMETRY,
        )
    else:
        raise ValueError(f"Unsupported ImageRegistrationMethod transform: {method}")

    registration.SetInitialTransform(transform, inPlace=False)
    final_transform = registration.Execute(fixed, moving)
    field = transform_to_displacement_array(final_transform, fixed)
    return field, float(registration.GetMetricValue())


def _require_itk_method(kind: RegistrationMethod) -> None:
    if itk is None:
        raise RuntimeError(f"Registration method '{kind}' requires the optional itk package.")
    required_class = {
        "curvature_registration": "CurvatureRegistrationFilter",
        "gpu_demons": "GPUDemonsRegistrationFilter",
    }.get(kind)
    if required_class is None:
        raise ValueError(f"Unsupported ITK registration method: {kind}")
    if not hasattr(itk, required_class):
        raise RuntimeError(f"The installed itk package does not provide {required_class}.")
    if kind == "gpu_demons" and not hasattr(itk, "GPUImage"):
        raise RuntimeError("The installed itk package has no GPU/OpenCL image support.")


def _itk_cpu_types() -> tuple[Any, Any]:
    if itk is None:
        raise RuntimeError("itk is not installed.")
    scalar_type = itk.Image[itk.F, 2]
    vector_type = itk.Image[itk.Vector[itk.F, 2], 2]
    return scalar_type, vector_type


def _itk_typed_image_from_array(array: np.ndarray, image_type: Any) -> Any:
    if itk is None:
        raise RuntimeError("itk is not installed.")
    value = np.asarray(array, dtype=np.float32)
    image = itk.image_from_array(value, is_vector=value.ndim == 3)
    image.SetSpacing((1.0, 1.0))
    return _itk_cast_image(image, image_type)


def _itk_cast_image(image: Any, output_type: Any) -> Any:
    if itk is None:
        raise RuntimeError("itk is not installed.")
    try:
        return itk.cast_image_filter(image, ttype=(type(image), output_type))
    except (KeyError, TypeError):
        caster = itk.CastImageFilter[type(image), output_type].New(Input=image)
        caster.Update()
        return caster.GetOutput()


def _run_itk_curvature_registration(
    fixed: np.ndarray,
    moving: np.ndarray,
    initial_field: np.ndarray,
    iterations: int,
) -> tuple[np.ndarray, float]:
    scalar_type, vector_type = _itk_cpu_types()
    fixed_image = _itk_typed_image_from_array(fixed, scalar_type)
    moving_image = _itk_typed_image_from_array(moving, scalar_type)
    initial_image = _itk_typed_image_from_array(initial_field, vector_type)
    filter_type = itk.CurvatureRegistrationFilter[
        scalar_type,
        scalar_type,
        vector_type,
    ]
    registration = filter_type.New(
        FixedImage=fixed_image,
        MovingImage=moving_image,
        InitialDisplacementField=initial_image,
    )
    registration.SetNumberOfIterations(int(iterations))
    registration.SetTimeStep(float(CURVATURE_TIME_STEP))
    registration.SetConstraintWeight(float(CURVATURE_CONSTRAINT_WEIGHT))
    registration.Update()
    field = np.asarray(itk.array_from_image(registration.GetOutput()), dtype=np.float32)
    return field, _numeric_filter_metric(registration)


def _run_itk_gpu_demons_registration(
    fixed: np.ndarray,
    moving: np.ndarray,
    initial_field: np.ndarray,
    iterations: int,
    field_sigma: float,
    update_sigma: float,
) -> tuple[np.ndarray, float]:
    _require_itk_method("gpu_demons")
    cpu_scalar_type, cpu_vector_type = _itk_cpu_types()
    gpu_scalar_type = itk.GPUImage[itk.F, 2]
    gpu_vector_type = itk.GPUImage[itk.Vector[itk.F, 2], 2]

    fixed_cpu = _itk_typed_image_from_array(fixed, cpu_scalar_type)
    moving_cpu = _itk_typed_image_from_array(moving, cpu_scalar_type)
    initial_cpu = _itk_typed_image_from_array(initial_field, cpu_vector_type)
    fixed_gpu = _itk_cast_image(fixed_cpu, gpu_scalar_type)
    moving_gpu = _itk_cast_image(moving_cpu, gpu_scalar_type)
    initial_gpu = _itk_cast_image(initial_cpu, gpu_vector_type)

    filter_type = itk.GPUDemonsRegistrationFilter[
        gpu_scalar_type,
        gpu_scalar_type,
        gpu_vector_type,
    ]
    registration = filter_type.New(
        FixedImage=fixed_gpu,
        MovingImage=moving_gpu,
        InitialDisplacementField=initial_gpu,
    )
    registration.SetNumberOfIterations(int(iterations))
    if hasattr(registration, "SetStandardDeviations"):
        registration.SetStandardDeviations(float(field_sigma))
    if hasattr(registration, "SetUpdateFieldStandardDeviations"):
        registration.SetUpdateFieldStandardDeviations(float(update_sigma))
    registration.Update()

    output_cpu = _itk_cast_image(registration.GetOutput(), cpu_vector_type)
    field = np.asarray(itk.array_from_image(output_cpu), dtype=np.float32)
    return field, _numeric_filter_metric(registration)


def run_itk_registration(
    fixed: np.ndarray,
    moving: np.ndarray,
    initial_field: np.ndarray,
    method: RegistrationMethod,
    iterations: int,
    field_sigma: float,
    update_sigma: float,
) -> tuple[np.ndarray, float]:
    if method == "curvature_registration":
        return _run_itk_curvature_registration(
            fixed,
            moving,
            initial_field,
            iterations,
        )
    if method == "gpu_demons":
        return _run_itk_gpu_demons_registration(
            fixed,
            moving,
            initial_field,
            iterations,
            field_sigma,
            update_sigma,
        )
    raise ValueError(f"Unsupported ITK registration method: {method}")


def estimate_registration_field(
    fixed_full: np.ndarray,
    moving_full: np.ndarray,
    mask_soft_full: np.ndarray,
    initial_full: np.ndarray | None,
    method: RegistrationMethod,
    scale: float,
    iterations: int,
    field_sigma: float,
    update_sigma: float,
    metric_radius: int,
    learning_rate: float,
    bspline_mesh_size: int,
) -> tuple[np.ndarray, float]:
    fixed_full = np.asarray(fixed_full, dtype=np.float32)
    moving_full = np.asarray(moving_full, dtype=np.float32)
    mask_soft_full = np.asarray(mask_soft_full, dtype=np.float32)
    if fixed_full.ndim != 2:
        raise ValueError(f"Fixed image must be 2-D, got {fixed_full.shape}.")
    if moving_full.shape != fixed_full.shape or mask_soft_full.shape != fixed_full.shape:
        raise ValueError(
            "Fixed, moving, and soft-mask shapes must match: "
            f"{fixed_full.shape}, {moving_full.shape}, {mask_soft_full.shape}."
        )

    height, width = fixed_full.shape
    small_shape = (
        max(1, int(round(height * float(scale)))),
        max(1, int(round(width * float(scale)))),
    )
    fixed_small = resize_scalar(fixed_full, small_shape)
    moving_small = resize_scalar(moving_full, small_shape)
    mask_small = np.clip(
        resize_scalar(
            mask_soft_full,
            small_shape,
            interpolation=cv2.INTER_LINEAR,
        ),
        0.0,
        1.0,
    ).astype(np.float32)
    fixed_small = fixed_small * mask_small + 0.5 * (1.0 - mask_small)
    moving_small = moving_small * mask_small + 0.5 * (1.0 - mask_small)

    if initial_full is None:
        initial_small = np.zeros((*small_shape, 2), dtype=np.float32)
    else:
        initial = np.asarray(initial_full, dtype=np.float32)
        if initial.shape != (*fixed_full.shape, 2):
            raise ValueError(f"Initial displacement shape must be (y, x, 2), got {initial.shape}.")
        initial_small = resize_displacement_to_full(initial, small_shape)
        initial_small *= mask_small[..., None]

    fixed_image = make_sitk_scalar(fixed_small)
    moving_image = make_sitk_scalar(moving_small)
    if method in ITK_REGISTRATION_METHODS:
        field_small, metric = run_itk_registration(
            fixed_small,
            moving_small,
            initial_small,
            method,
            iterations,
            field_sigma,
            update_sigma,
        )
    elif method in PDE_REGISTRATION_METHODS:
        field_small, metric = run_pde_registration(
            fixed_image,
            moving_image,
            initial_small,
            method,
            iterations,
            field_sigma,
        )
    else:
        field_small, metric = run_image_registration_method(
            fixed_image,
            moving_image,
            initial_small,
            method,
            iterations,
            field_sigma,
            update_sigma,
            metric_radius,
            learning_rate,
            bspline_mesh_size,
        )

    field = resize_displacement_to_full(field_small, fixed_full.shape)
    field *= mask_soft_full[..., None]
    return field.astype(np.float32), metric


def _numeric_filter_metric(registration: Any) -> float:
    for getter_name in ("GetMetric", "GetMetricValue", "GetRMSChange"):
        getter = getattr(registration, getter_name, None)
        if getter is None:
            continue
        try:
            return float(getter())
        except (TypeError, ValueError, RuntimeError):
            continue
    return float("nan")


def _require_simpleitk() -> None:
    if sitk is None:
        raise RuntimeError("SimpleITK is not installed.")


def _require_opencv() -> None:
    if cv2 is None:
        raise RuntimeError("OpenCV is not installed.")
