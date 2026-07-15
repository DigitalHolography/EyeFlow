"""
Command-line interface to run EyeFlow pipelines over HOLO selections.

Usage example:
    python cli.py --data data/

Inputs:
    --data / -d        Directory, single .holo, .txt path list, or .zip archive of HOLO data.
    --pipelines / -p   Optional target list; enabled user/default settings are used otherwise.
    --output / -o      Optional output root; source-adjacent GUI defaults are used otherwise.
    --zip / -z         When set, compress the outputs into a .zip archive after completion.
    --zip-name         Optional filename for the archive (default: outputs.zip).
"""

import argparse
import shutil
import sys
import tempfile
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from uuid import uuid4

from app_settings import AppSettingsStore
from runtime_limits import configure_numeric_threads

configure_numeric_threads()

from input_output import (
    create_zip_from_tree,
)
from pipelines import (
    PipelineDescriptor,
    load_pipeline_catalog,
)
from pipeline_engine import (
    execute_run,
    expand_run_inputs,
    resolve_run_spec,
    selectable_pipeline_registry,
)

def _build_pipeline_registry() -> dict[str, PipelineDescriptor]:
    available, missing = load_pipeline_catalog()
    return {pipeline.name: pipeline for pipeline in (*available, *missing)}


def _load_pipeline_targets(
    path: Path, registry: dict[str, PipelineDescriptor]
) -> tuple[str, ...]:
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

    return tuple(dict.fromkeys(selected_names))


def _load_configured_pipeline_targets(
    registry: dict[str, PipelineDescriptor],
    settings_store: AppSettingsStore | None = None,
) -> tuple[str, ...]:
    store = settings_store or AppSettingsStore()
    visibility = store.load_pipeline_visibility()
    selected_names = tuple(
        name for name in registry if visibility.get(name, False)
    )
    if selected_names:
        return selected_names
    raise ValueError(
        "No visible pipelines are enabled in EyeFlow settings. "
        "Enable at least one pipeline in the GUI or provide --pipelines."
    )


def _default_cli_output_root(data_path: Path) -> Path:
    source = data_path.expanduser().resolve()
    return source if source.is_dir() else source.parent


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
    staging_zip = zip_path.with_name(
        f".{zip_path.name}.eyeflow-staging-{uuid4().hex}"
    )
    try:
        create_zip_from_tree(
            folder,
            staging_zip,
            progress_callback=progress_callback,
        )
        staging_zip.replace(zip_path)
    finally:
        if staging_zip.exists():
            staging_zip.unlink()
    return zip_path


def run_cli(
    data_path: Path,
    pipelines_file: Path | None = None,
    output_dir: Path | None = None,
    zip_outputs: bool = False,
    zip_name: str | None = None,
) -> int:
    registry = _build_pipeline_registry()
    target_registry = selectable_pipeline_registry(registry.values())
    target_names = (
        _load_pipeline_targets(pipelines_file, target_registry)
        if pipelines_file is not None
        else _load_configured_pipeline_targets(target_registry)
    )
    work_tempdir_path: Path | None = None
    clean_work_output = False
    with expand_run_inputs(data_path) as expanded_inputs:
        output_root = (
            output_dir.expanduser().resolve() if output_dir is not None else None
        )
        input_is_zip = data_path.suffix.lower() == ".zip"
        if output_root is None and (input_is_zip or zip_outputs):
            output_root = _default_cli_output_root(data_path)
        if output_root is not None:
            output_root.mkdir(parents=True, exist_ok=True)

        work_root = output_root
        if zip_outputs:
            assert output_root is not None
            work_tempdir_path = Path(tempfile.mkdtemp(dir=output_root))
            work_root = work_tempdir_path

        spec = resolve_run_spec(
            input_paths=expanded_inputs.paths,
            target_names=target_names,
            pipelines=registry.values(),
            output_root=work_root,
            batch_root=expanded_inputs.batch_root,
        )
        result = execute_run(spec, on_log=print)
        processed_outputs = list(result.outputs)
        archive_failed = False

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

                assert work_root is not None
                assert output_root is not None
                zip_path = _zip_output_dir(
                    work_root,
                    target_path=output_root / final_name,
                    progress_callback=_zip_progress,
                )
                print(f"[ZIP] Archive created: {zip_path}")
                summary_msg = f"ZIP archive: {zip_path}"
                clean_work_output = True
            except Exception as exc:  # noqa: BLE001
                archive_failed = True
                print(
                    f"[ZIP FAIL] Could not create ZIP archive: {exc}", file=sys.stderr
                )
                summary_msg = f"Outputs stored under: {work_root}"
        else:
            if len(processed_outputs) == 1:
                summary_msg = f"Output file: {processed_outputs[0]}"
            elif work_root is None:
                summary_msg = (
                    f"{len(processed_outputs)} output files written beside their inputs"
                )
            else:
                summary_msg = f"Outputs stored under: {work_root}"

        print(f"Completed. {summary_msg}")

        if result.failures:
            print(f"{len(result.failures)} failure(s):", file=sys.stderr)
            for failure in result.failures:
                print(f" - {failure}", file=sys.stderr)
        exit_code = 1 if result.failures or archive_failed else 0
        if clean_work_output and work_tempdir_path is not None:
            shutil.rmtree(work_tempdir_path, ignore_errors=True)
        return exit_code


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
            "a .txt path list, or a .zip archive."
        ),
    )
    parser.add_argument(
        "-p",
        "--pipelines",
        type=Path,
        help=(
            "Optional text file with pipeline targets. When omitted, enabled "
            "targets are loaded from user settings or default_settings.json."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "Optional output root. By default, normal inputs use GUI-style "
            "source-adjacent outputs; ZIP inputs use the ZIP's parent folder."
        ),
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
