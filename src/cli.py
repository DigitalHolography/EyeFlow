"""
Command-line interface to run EyeFlow pipelines over HOLO selections.

Usage example:
    python cli.py --data data/ --pipelines pipelines.txt --output ./results --zip --zip-name my_run.zip

Inputs:
    --data / -d        Directory, single .holo, .txt stem list, or .zip archive of HOLO data.
    --pipelines / -p   Text file listing pipeline target names (one per line, '#' and blank lines ignored).
    --output / -o      Base directory where results will be written (input subfolder layout is preserved).
    --zip / -z         When set, compress the outputs into a .zip archive after completion.
    --zip-name         Optional filename for the archive (default: outputs.zip).
"""

import argparse
import os
import shutil
import sys
import tempfile
import time
import zipfile
from collections.abc import Callable, Sequence
from pathlib import Path

from runtime_limits import configure_numeric_threads

configure_numeric_threads()

from input_output import (
    INPUT_LIST_SUFFIX,
    create_zip_from_tree,
    resolve_selected_run_layouts,
)
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipelines import (
    PipelineDescriptor,
    load_pipeline_catalog,
)
from pipeline_engine import (
    PipelineDAG,
    PipelineExecutionPlan,
    run_pipelines_to_output,
)

HOLO_SUFFIX = ".holo"


def _build_pipeline_registry() -> dict[str, PipelineDescriptor]:
    available, missing = load_pipeline_catalog()
    return {pipeline.name: pipeline for pipeline in (*available, *missing)}


def _load_pipeline_plan(
    path: Path, registry: dict[str, PipelineDescriptor]
) -> PipelineExecutionPlan:
    raw_lines = path.read_text(encoding="utf-8").splitlines()
    selected_names: list[str] = []
    missing: list[str] = []
    for line in raw_lines:
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        pipeline = registry.get(name)
        if pipeline is None:
            missing.append(name)
        else:
            selected_names.append(pipeline.name)
    if missing:
        available = ", ".join(registry.keys())
        raise ValueError(
            f"Unknown pipeline(s): {', '.join(missing)}. Available: {available}"
        )
    if not selected_names:
        raise ValueError(
            "No pipelines selected (file is empty or only contains comments)."
        )

    plan = PipelineDAG(registry.values()).resolve_targets(selected_names)
    unavailable = [pipeline for pipeline in plan.descriptors if not pipeline.available]
    if unavailable:
        details = []
        for pipeline in unavailable:
            reason = ", ".join(pipeline.missing_deps or pipeline.requires)
            details.append(f"{pipeline.name}" + (f" ({reason})" if reason else ""))
        raise ValueError(
            "The DAG requires unavailable pipeline(s): " + ", ".join(details)
        )
    return plan


def _find_holo_inputs(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() in {HOLO_SUFFIX, INPUT_LIST_SUFFIX}:
            return [path]
        raise ValueError(
            f"File is not a {HOLO_SUFFIX} or {INPUT_LIST_SUFFIX} file: {path}"
        )
    if path.is_dir():
        return sorted(path.rglob(f"*{HOLO_SUFFIX}"))
    raise FileNotFoundError(f"Input path does not exist: {path}")


def _prepare_data_root(
    data_path: Path,
) -> tuple[Path, tempfile.TemporaryDirectory | None]:
    """Return a directory containing HOLO files; extract zip archives when needed."""
    if data_path.is_file() and data_path.suffix.lower() == ".zip":
        tempdir = tempfile.TemporaryDirectory()
        with zipfile.ZipFile(data_path, "r") as zf:
            zf.extractall(tempdir.name)
        return Path(tempdir.name), tempdir
    return data_path, None


def _unique_output_root(path: Path) -> Path:
    if not path.exists():
        return path
    suffix = 1
    while True:
        candidate = path.with_name(f"{path.name}_{suffix}")
        if not candidate.exists():
            return candidate
        suffix += 1


def _batch_root(holo_paths: Sequence[Path]) -> Path:
    if not holo_paths:
        return Path.cwd()
    if len(holo_paths) == 1:
        return holo_paths[0].parent
    try:
        return Path(os.path.commonpath([str(path.parent) for path in holo_paths]))
    except ValueError:
        return Path.cwd()


def _relative_to_batch(holo_path: Path, batch_root: Path) -> Path:
    try:
        return holo_path.relative_to(batch_root)
    except ValueError:
        anchor = Path(holo_path.anchor)
        drive_token = holo_path.drive.rstrip(":\\/") or "root"
        tail = holo_path.relative_to(anchor) if anchor != holo_path else Path()
        return Path(drive_token) / tail


def _run_pipelines_on_input(
    run_layout: HoloRunLayout,
    relative_holo_path: Path,
    plan: PipelineExecutionPlan,
    output_root: Path,
) -> Path:
    target_dir = output_root / relative_holo_path.parent
    target_dir.mkdir(parents=True, exist_ok=True)
    output_layout = HoloRunLayout.from_holo(
        run_layout.holo_path,
        output_root=target_dir,
    )
    output_layout = output_layout.with_root_dir(
        _unique_output_root(output_layout.root_dir)
    )
    output_manager = OutputManager(output_layout)
    output_h5_path = output_manager.path_for(OutputType.H5)
    print(f"[INPUT] HOLO -> {run_layout.holo_path}")
    print(f"[RESOLVED] HD -> {run_layout.hd_h5}")
    print(f"[RESOLVED] DV -> {run_layout.dv_h5}")
    run_pipelines_to_output(
        output_manager=output_manager,
        pipelines=plan.descriptors,
        target_names=plan.targets,
        holodoppler_h5=run_layout.hd_h5,
        doppler_vision_h5=run_layout.dv_h5,
        on_pipeline_success=lambda name: print(
            f"[OK] {run_layout.holo_path.name} -> {name}"
        ),
    )
    print(f"[OK] {run_layout.holo_path.name}: output -> {output_h5_path}")
    return output_h5_path


def _zip_output_dir(
    folder: Path,
    target_path: Path | None = None,
    progress_callback: Callable[[int, int, Path], None] | None = None,
) -> Path:
    folder = folder.expanduser().resolve()
    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Output folder does not exist: {folder}")
    if target_path is None:
        zip_name = f"{folder.name}_outputs.zip" if folder.name else "outputs.zip"
        zip_path = folder.parent / zip_name
    else:
        zip_path = target_path.expanduser().resolve()
    if zip_path.exists():
        zip_path.unlink()
    create_zip_from_tree(folder, zip_path, progress_callback=progress_callback)
    return zip_path


def run_cli(
    data_path: Path,
    pipelines_file: Path,
    output_dir: Path,
    zip_outputs: bool = False,
    zip_name: str | None = None,
) -> int:
    registry = _build_pipeline_registry()
    plan = _load_pipeline_plan(pipelines_file, registry)
    data_root, tempdir = _prepare_data_root(data_path)
    work_tempdir_path: Path | None = None
    clean_work_output = False
    try:
        inputs = _find_holo_inputs(data_root)
        if not inputs:
            raise ValueError(f"No {HOLO_SUFFIX} files found under {data_path}")
        run_layouts = resolve_selected_run_layouts(inputs)
        batch_root = _batch_root(inputs)

        output_root = output_dir.expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)

        work_root = output_root
        if zip_outputs:
            work_tempdir_path = Path(tempfile.mkdtemp(dir=output_root))
            work_root = work_tempdir_path

        failures: list[str] = []
        processed_outputs: list[Path] = []
        for run_layout in run_layouts:
            relative_holo_path = _relative_to_batch(run_layout.holo_path, batch_root)
            try:
                combined_output = _run_pipelines_on_input(
                    run_layout,
                    relative_holo_path,
                    plan,
                    work_root,
                )
                processed_outputs.append(combined_output)
            except Exception as exc:  # noqa: BLE001
                failures.append(f"{run_layout.holo_path}: {exc}")
                print(
                    f"[FAIL] {run_layout.holo_path.name}: {exc}",
                    file=sys.stderr,
                )

        if zip_outputs:
            try:
                final_name = (zip_name or "outputs.zip").strip() or "outputs.zip"
                if not final_name.lower().endswith(".zip"):
                    final_name += ".zip"
                print("[ZIP] Preparing archive...")
                last_progress_log = 0.0

                def _zip_progress(done: int, total: int, _rel_path: Path) -> None:
                    nonlocal last_progress_log
                    now = time.monotonic()
                    if done == total or (now - last_progress_log) >= 0.5:
                        pct = 100 if total == 0 else int((done * 100) / total)
                        print(f"[ZIP] {done}/{total} files ({pct}%)")
                        last_progress_log = now

                zip_path = _zip_output_dir(
                    work_root,
                    target_path=output_root / final_name,
                    progress_callback=_zip_progress,
                )
                print(f"[ZIP] Archive created: {zip_path}")
                summary_msg = f"ZIP archive: {zip_path}"
                clean_work_output = True
            except Exception as exc:  # noqa: BLE001
                print(
                    f"[ZIP FAIL] Could not create ZIP archive: {exc}", file=sys.stderr
                )
                summary_msg = f"Outputs stored under: {work_root}"
        else:
            if len(processed_outputs) == 1:
                summary_msg = f"Output file: {processed_outputs[0]}"
            else:
                summary_msg = f"Outputs stored under: {work_root}"

        print(f"Completed. {summary_msg}")

        if failures:
            print(f"{len(failures)} failure(s):", file=sys.stderr)
            for msg in failures:
                print(f" - {msg}", file=sys.stderr)
            return 1
        return 0
    finally:
        if tempdir is not None:
            tempdir.cleanup()
        if clean_work_output and work_tempdir_path is not None:
            shutil.rmtree(work_tempdir_path, ignore_errors=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run EyeFlow pipelines over one or more HOLO selections."
    )
    parser.add_argument(
        "-d",
        "--data",
        required=True,
        type=Path,
        help=(
            "Directory containing .holo files, a single .holo file, "
            "a .txt stem list, or a .zip archive."
        ),
    )
    parser.add_argument(
        "-p",
        "--pipelines",
        required=True,
        type=Path,
        help="Text file with pipeline targets to resolve through the DAG.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Base output directory. Input subfolder layout is preserved for output files.",
    )
    parser.add_argument(
        "-z",
        "--zip",
        action="store_true",
        help="Zip the outputs after processing (only the archive is kept).",
    )
    parser.add_argument(
        "--zip-name",
        type=str,
        default="outputs.zip",
        help="Archive filename to place inside the output directory (default: outputs.zip).",
    )
    args = parser.parse_args(argv)

    try:
        return run_cli(
            args.data,
            args.pipelines,
            args.output,
            zip_outputs=args.zip,
            zip_name=args.zip_name,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
