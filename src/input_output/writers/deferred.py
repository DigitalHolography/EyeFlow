"""Scoped figure capture; the pipeline engine owns scheduling and completion."""

from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

from input_output.artifact_ownership import register_artifact_path

_CAPTURE = ContextVar("eyeflow_artifact_capture", default=None)


@contextmanager
def capture_artifact_writes(capture):
    token = _CAPTURE.set(capture)
    try:
        yield
    finally:
        _CAPTURE.reset(token)


def defer_write(kind, payload, paths, **kwargs) -> bool:
    for path in paths:
        register_artifact_path(path, writing=True)
    capture = _CAPTURE.get()
    return capture is not None and capture(kind, payload, [Path(p) for p in paths], kwargs)


def save_figure_file(path, fig, **kwargs) -> Path:
    """Save through the managed exporter, or immediately outside a runtime."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not defer_write("figure", fig, [path], **kwargs):
        fig.savefig(path, **kwargs)
    return path


def snapshot_figure_write(path, kind, payload, relative_paths, kwargs):
    """Freeze artists/data and formatting now, before the caller mutates/closes them."""
    from joblib.externals import cloudpickle
    import matplotlib

    rc = dict(matplotlib.rcParams) if kind != "array" else None
    if rc is not None:
        rc.pop("backend", None)
        rc.pop("interactive", None)
    with Path(path).open("wb") as stream:
        cloudpickle.dump((kind, payload, relative_paths, kwargs, rc), stream)


def export_figure_writes(snapshot, manager, *, errors=()):
    """Render one pipeline's frozen figures in an isolated Agg worker process."""
    if errors:
        raise RuntimeError("\n".join(errors))
    from joblib.externals import cloudpickle
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib import pyplot as plt

    from .article import render_article_figure_pair
    from .png import write_png_file

    paths = []
    for source in sorted(Path(snapshot).glob("*.pickle")):
        with source.open("rb") as stream:
            kind, payload, relatives, kwargs, rc = cloudpickle.load(stream)
        targets = [manager.layout.ef_dir / relative for relative in relatives]
        for target in targets:
            target.parent.mkdir(parents=True, exist_ok=True)
        if kind == "array":
            write_png_file(targets[0], payload)
        else:
            try:
                FigureCanvasAgg(payload)
                with matplotlib.rc_context(rc):
                    if kind == "article_pair":
                        render_article_figure_pair(payload, targets, **kwargs)
                    else:
                        payload.savefig(targets[0], **kwargs)
            finally:
                plt.close(payload)
        paths.extend(targets)
    return paths
