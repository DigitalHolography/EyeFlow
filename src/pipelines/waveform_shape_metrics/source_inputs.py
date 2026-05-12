"""Read HD/DV inputs for the waveform-shape metrics pipeline."""

from pipeline_engine.imports import np

from .constants import DOPPLERVIEW_DEFAULT_LOCAL_BACKGROUND_DIST


def dopplerview_cache_from_h5(ctx) -> dict[str, object]:
    return {
        "moment0": _read_hd_moment(ctx, "moment0"),
        "moment2": _read_hd_moment(ctx, "moment2"),
        "retinal_artery_mask": ctx.dv.array("retinal_artery_mask", dtype=bool),
        "retinal_vein_mask": ctx.dv.array("retinal_vein_mask", dtype=bool),
        "retinal_labeled_vessels": ctx.dv.array(
            "retinal_labeled_vessels",
            dtype=np.int32,
            default=None,
        ),
    }


def local_background_dist(ctx) -> int:
    return int(
        ctx.dv.config(
            "local_background_dist",
            default=DOPPLERVIEW_DEFAULT_LOCAL_BACKGROUND_DIST,
        )
    )


def _read_hd_moment(ctx, key: str) -> np.ndarray:
    return _coerce_dopplerview_moment(ctx.hd.array(key, dtype=np.float32))


def _coerce_dopplerview_moment(value) -> np.ndarray:
    squeezed = np.squeeze(np.asarray(value, dtype=np.float32))
    if squeezed.ndim != 3:
        raise ValueError(
            "Holodoppler moment datasets must become 3-D after squeeze, "
            f"got shape {squeezed.shape}."
        )
    return np.transpose(squeezed, (0, 2, 1))
