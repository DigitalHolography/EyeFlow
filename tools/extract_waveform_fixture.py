#!/usr/bin/env python3
"""Extract compact, metadata-stripped waveform fixtures from EyeFlow HDF5 outputs.

The resulting files intentionally contain only allow-listed numeric waveform,
beat-timing, and optional per-beat datasets. Source filenames, HDF5 attributes,
images, masks, maps, and acquisition metadata are not copied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np


SCHEMA_NAME = "eyeflow-waveform-fixture"
SCHEMA_VERSION = 1
DEFAULT_MAX_OUTPUT_MB = 64.0
H5_SUFFIXES = frozenset({".h5", ".hdf5"})


@dataclass(frozen=True)
class DatasetSpec:
    output_path: str
    aliases: tuple[str, ...]
    required: bool = False


@dataclass(frozen=True)
class PreparedFixture:
    arrays: dict[str, np.ndarray]
    source_paths: dict[str, str]
    dt_seconds: float
    content_sha256: str
    estimated_bytes: int


@dataclass(frozen=True)
class ExtractionResult:
    input_path: Path
    output_path: Path
    created: bool
    dataset_count: int
    estimated_bytes: int


DATASET_SPECS = (
    DatasetSpec(
        "waveforms/artery/raw",
        (
            "analysis/retinal_artery_velocity_signal",
            "artery/velocity/signal/value",
        ),
        required=True,
    ),
    DatasetSpec(
        "waveforms/vein/raw",
        (
            "analysis/retinal_vein_velocity_signal",
            "vein/velocity/signal/value",
        ),
        required=True,
    ),
    DatasetSpec(
        "beats/boundary_indices",
        (
            "analysis/beat_indices",
            "perbeat/beat_indices/value",
        ),
        required=True,
    ),
    DatasetSpec(
        "beats/period_seconds",
        (
            "analysis/time_per_beat",
            "perbeat/time_per_beat/value",
            "Artery/VelocityPerBeat/beatPeriodSeconds/value",
            "perbeat/beat_period_seconds/value",
        ),
        required=True,
    ),
    DatasetSpec(
        "beats/period_samples",
        (
            "Artery/VelocityPerBeat/beatPeriodIdx/value",
            "perbeat/beat_period_idx/value",
        ),
    ),
    DatasetSpec(
        "waveforms/artery/filtered",
        (
            "analysis/velocitysignal_filtered",
            "artery/velocity/filtered_signal/value",
        ),
    ),
    DatasetSpec(
        "per_beat/artery/filtered",
        (
            "analysis/velocitysignal_per_beat",
            "artery/velocity/perbeat/filtered_signal/value",
        ),
    ),
    DatasetSpec(
        "per_beat/artery/raw",
        (
            "Artery/VelocityPerBeat/VelocitySignalPerBeat/value",
            "artery/velocity/perbeat/signal/value",
        ),
    ),
    DatasetSpec(
        "per_beat/artery/band_limited",
        (
            "Artery/VelocityPerBeat/VelocitySignalPerBeatBandLimited/value",
            "artery/velocity/perbeat/band_limited/value",
        ),
    ),
    DatasetSpec(
        "per_beat/vein/raw",
        (
            "Vein/VelocityPerBeat/VelocitySignalPerBeat/value",
            "vein/velocity/perbeat/signal/value",
        ),
    ),
    DatasetSpec(
        "per_beat/vein/band_limited",
        (
            "Vein/VelocityPerBeat/VelocitySignalPerBeatBandLimited/value",
            "vein/velocity/perbeat/band_limited/value",
        ),
    ),
    DatasetSpec(
        "per_beat/artery/segments/raw",
        (
            "Artery/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegment/value",
            "artery/velocity/perbeat/segments/signal/value",
        ),
    ),
    DatasetSpec(
        "per_beat/artery/segments/band_limited",
        (
            "Artery/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegmentBandLimited/value",
            "artery/velocity/perbeat/segments/band_limited/value",
        ),
    ),
    DatasetSpec(
        "per_beat/vein/segments/raw",
        (
            "Vein/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegment/value",
            "vein/velocity/perbeat/segments/signal/value",
        ),
    ),
    DatasetSpec(
        "per_beat/vein/segments/band_limited",
        (
            "Vein/VelocityPerBeat/Segments/"
            "VelocitySignalPerBeatPerSegmentBandLimited/value",
            "vein/velocity/perbeat/segments/band_limited/value",
        ),
    ),
)


class FixtureExtractionError(ValueError):
    """Raised when an HDF5 file cannot produce a valid waveform fixture."""


def prepare_fixture(input_path: Path, max_output_mb: float) -> PreparedFixture:
    """Read and validate allow-listed datasets without writing an output file."""
    arrays: dict[str, np.ndarray] = {}
    source_paths: dict[str, str] = {}
    missing: list[str] = []

    try:
        with h5py.File(input_path, "r") as source:
            for spec in DATASET_SPECS:
                source_path = _first_existing_path(source, spec.aliases)
                if source_path is None:
                    if spec.required:
                        missing.append(f"{spec.output_path} ({' or '.join(spec.aliases)})")
                    continue
                arrays[spec.output_path] = _read_numeric_dataset(
                    source[source_path],
                    source_path,
                )
                source_paths[spec.output_path] = source_path
    except OSError as exc:
        raise FixtureExtractionError(f"Could not open HDF5 file: {exc}") from exc

    if missing:
        raise FixtureExtractionError(
            "Missing required EyeFlow waveform datasets: " + "; ".join(missing)
        )

    dt_seconds = _validate_and_normalize(arrays)
    estimated_bytes = sum(int(array.nbytes) for array in arrays.values())
    max_bytes = int(max_output_mb * 1024 * 1024)
    if estimated_bytes > max_bytes:
        raise FixtureExtractionError(
            f"Selected datasets total {estimated_bytes / 1024 / 1024:.2f} MB, "
            f"above the {max_output_mb:.2f} MB safety limit."
        )

    content_sha256 = _content_hash(arrays, dt_seconds)
    return PreparedFixture(
        arrays=arrays,
        source_paths=source_paths,
        dt_seconds=dt_seconds,
        content_sha256=content_sha256,
        estimated_bytes=estimated_bytes,
    )


def extract_fixture(
    input_path: Path,
    output_dir: Path,
    *,
    max_output_mb: float = DEFAULT_MAX_OUTPUT_MB,
    overwrite: bool = False,
) -> ExtractionResult:
    """Create one compact fixture and return its output metadata."""
    prepared = prepare_fixture(input_path, max_output_mb=max_output_mb)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"waveform_fixture_{prepared.content_sha256[:12]}.h5"

    if output_path.exists() and not overwrite:
        if _existing_fixture_matches(output_path, prepared.content_sha256):
            return ExtractionResult(
                input_path=input_path,
                output_path=output_path,
                created=False,
                dataset_count=len(prepared.arrays),
                estimated_bytes=prepared.estimated_bytes,
            )
        raise FixtureExtractionError(
            f"Output path already exists with different content: {output_path}"
        )

    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=".waveform_fixture_",
        suffix=".tmp.h5",
        dir=output_dir,
    )
    os.close(file_descriptor)
    temporary_path = Path(temporary_name)
    try:
        _write_fixture(temporary_path, prepared)
        temporary_path.replace(output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    return ExtractionResult(
        input_path=input_path,
        output_path=output_path,
        created=True,
        dataset_count=len(prepared.arrays),
        estimated_bytes=prepared.estimated_bytes,
    )


def discover_h5_files(
    input_paths: Sequence[Path],
    *,
    recursive: bool,
    output_dir: Path,
) -> list[Path]:
    """Resolve files and directories into a stable, duplicate-free input list."""
    discovered: dict[Path, None] = {}
    resolved_output = output_dir.resolve()

    for raw_path in input_paths:
        path = raw_path.expanduser()
        if not path.exists():
            raise FixtureExtractionError(f"Input path does not exist: {path}")
        if path.is_file():
            if path.suffix.lower() not in H5_SUFFIXES:
                raise FixtureExtractionError(f"Input is not an HDF5 file: {path}")
            resolved = path.resolve()
            if not _is_generated_fixture(resolved, resolved_output):
                discovered[resolved] = None
            continue

        iterator: Iterable[Path]
        iterator = path.rglob("*") if recursive else path.glob("*")
        for candidate in iterator:
            if not candidate.is_file() or candidate.suffix.lower() not in H5_SUFFIXES:
                continue
            resolved = candidate.resolve()
            if not _is_generated_fixture(resolved, resolved_output):
                discovered[resolved] = None

    return sorted(discovered, key=lambda item: str(item).casefold())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create compact, metadata-stripped waveform fixtures from EyeFlow "
            "HDF5 outputs."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="EyeFlow H5 file(s) or directories containing EyeFlow H5 files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for deidentified waveform fixtures.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search input directories recursively.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate compatible inputs without writing fixtures.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing fixture with the same content-derived name.",
    )
    parser.add_argument(
        "--max-output-mb",
        type=float,
        default=DEFAULT_MAX_OUTPUT_MB,
        help=(
            "Maximum uncompressed allow-listed data per fixture "
            f"(default: {DEFAULT_MAX_OUTPUT_MB:g} MB)."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not np.isfinite(args.max_output_mb) or args.max_output_mb <= 0:
        parser.error("--max-output-mb must be a positive finite number.")

    try:
        files = discover_h5_files(
            args.inputs,
            recursive=args.recursive,
            output_dir=args.output_dir,
        )
    except FixtureExtractionError as exc:
        parser.error(str(exc))

    if not files:
        print("No HDF5 files found.", file=sys.stderr)
        return 1

    created = 0
    reused = 0
    failed = 0
    for input_path in files:
        try:
            if args.dry_run:
                prepared = prepare_fixture(input_path, args.max_output_mb)
                print(
                    f"compatible: {input_path} "
                    f"({len(prepared.arrays)} datasets, "
                    f"{prepared.estimated_bytes / 1024:.1f} KiB)"
                )
                continue
            result = extract_fixture(
                input_path,
                args.output_dir,
                max_output_mb=args.max_output_mb,
                overwrite=args.overwrite,
            )
            status = "created" if result.created else "reused"
            print(
                f"{status}: {result.output_path.name} "
                f"({result.dataset_count} datasets, "
                f"{result.estimated_bytes / 1024:.1f} KiB)"
            )
            created += int(result.created)
            reused += int(not result.created)
        except FixtureExtractionError as exc:
            print(f"failed: {input_path}: {exc}", file=sys.stderr)
            failed += 1

    if args.dry_run:
        print(f"Validated {len(files) - failed} compatible file(s); {failed} failed.")
    else:
        print(f"Created {created}, reused {reused}, failed {failed}.")
    return 1 if failed else 0


def _first_existing_path(source: h5py.File, aliases: Sequence[str]) -> str | None:
    return next((path for path in aliases if path in source), None)


def _read_numeric_dataset(dataset: h5py.Dataset, source_path: str) -> np.ndarray:
    if dataset.dtype.kind not in "biufc":
        raise FixtureExtractionError(
            f"Dataset {source_path} must be numeric, got dtype {dataset.dtype}."
        )
    return np.asarray(dataset[()])


def _validate_and_normalize(arrays: dict[str, np.ndarray]) -> float:
    artery = np.asarray(arrays["waveforms/artery/raw"]).reshape(-1)
    vein = np.asarray(arrays["waveforms/vein/raw"]).reshape(-1)
    if artery.size < 3 or vein.size < 3:
        raise FixtureExtractionError("Artery and vein waveforms must contain at least 3 samples.")
    if artery.size != vein.size:
        raise FixtureExtractionError(
            f"Artery and vein waveform lengths differ ({artery.size} vs {vein.size})."
        )
    if not np.any(np.isfinite(artery)) or not np.any(np.isfinite(vein)):
        raise FixtureExtractionError("Artery and vein waveforms must contain finite samples.")

    raw_indices = np.asarray(arrays["beats/boundary_indices"]).reshape(-1)
    if raw_indices.size < 2 or not np.all(np.isfinite(raw_indices)):
        raise FixtureExtractionError("Beat boundaries must contain at least two finite indexes.")
    rounded_indices = np.rint(raw_indices).astype(np.int64)
    if not np.allclose(raw_indices, rounded_indices):
        raise FixtureExtractionError("Beat boundary indexes must be integer-valued.")
    if rounded_indices[0] < 0 or rounded_indices[-1] > artery.size:
        raise FixtureExtractionError(
            "Beat boundary indexes must stay within the waveform sample range."
        )
    periods_samples = np.diff(rounded_indices)
    if np.any(periods_samples <= 0):
        raise FixtureExtractionError("Beat boundary indexes must be strictly increasing.")

    periods_seconds = np.asarray(arrays["beats/period_seconds"], dtype=np.float64).reshape(-1)
    if periods_seconds.size != periods_samples.size:
        raise FixtureExtractionError(
            "Beat period count must equal the number of boundary intervals."
        )
    if not np.all(np.isfinite(periods_seconds)) or np.any(periods_seconds <= 0):
        raise FixtureExtractionError("Beat periods must be positive and finite.")

    dt_values = periods_seconds / periods_samples
    dt_seconds = float(np.median(dt_values))
    if not np.isfinite(dt_seconds) or dt_seconds <= 0:
        raise FixtureExtractionError("Could not derive a positive sample interval.")

    arrays["waveforms/artery/raw"] = artery
    arrays["waveforms/vein/raw"] = vein
    arrays["beats/boundary_indices"] = rounded_indices.astype(np.int32)
    arrays["beats/period_seconds"] = periods_seconds.astype(np.float32)
    arrays["beats/period_samples"] = periods_samples.astype(np.int32)
    return dt_seconds


def _content_hash(arrays: dict[str, np.ndarray], dt_seconds: float) -> str:
    digest = hashlib.sha256()
    digest.update(f"{SCHEMA_NAME}:{SCHEMA_VERSION}".encode("ascii"))
    digest.update(np.float64(dt_seconds).tobytes())
    for output_path in sorted(arrays):
        array = np.ascontiguousarray(arrays[output_path])
        digest.update(output_path.encode("utf-8"))
        digest.update(array.dtype.str.encode("ascii"))
        digest.update(json.dumps(array.shape).encode("ascii"))
        digest.update(array.tobytes())
    return digest.hexdigest()


def _existing_fixture_matches(path: Path, content_sha256: str) -> bool:
    try:
        with h5py.File(path, "r") as fixture:
            return (
                fixture.attrs.get("fixture_schema") == SCHEMA_NAME
                and int(fixture.attrs.get("fixture_schema_version", -1)) == SCHEMA_VERSION
                and fixture.attrs.get("content_sha256") == content_sha256
            )
    except (OSError, TypeError, ValueError):
        return False


def _write_fixture(path: Path, prepared: PreparedFixture) -> None:
    with h5py.File(path, "w") as fixture:
        fixture.attrs["fixture_schema"] = SCHEMA_NAME
        fixture.attrs["fixture_schema_version"] = np.int32(SCHEMA_VERSION)
        fixture.attrs["content_sha256"] = prepared.content_sha256
        fixture.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        fixture.attrs["metadata_stripped"] = True
        fixture.attrs["contains_patient_derived_data"] = True
        fixture.attrs["privacy_warning"] = (
            "Clinical waveforms remain potentially sensitive even after metadata removal."
        )

        dt_dataset = fixture.create_dataset(
            "beats/dt_seconds",
            data=np.float64(prepared.dt_seconds),
        )
        dt_dataset.attrs["unit"] = "s"
        source_map = fixture.create_dataset(
            "metadata/source_dataset_paths_json",
            data=json.dumps(prepared.source_paths, sort_keys=True),
        )
        source_map.attrs["description"] = (
            "Schema paths used in the source HDF5; no filename or source attributes."
        )

        for output_path, array in prepared.arrays.items():
            options = _compression_options(array)
            dataset = fixture.create_dataset(output_path, data=array, **options)
            if output_path == "beats/boundary_indices":
                dataset.attrs["index_base"] = np.int32(0)
                dataset.attrs["unit"] = "sample"
            elif output_path == "beats/period_samples":
                dataset.attrs["unit"] = "sample"
            elif output_path == "beats/period_seconds":
                dataset.attrs["unit"] = "s"


def _compression_options(array: np.ndarray) -> dict[str, object]:
    if array.ndim == 0 or array.size == 0:
        return {}
    return {"compression": "gzip", "compression_opts": 4, "shuffle": True}


def _is_inside(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _is_generated_fixture(path: Path, output_dir: Path) -> bool:
    if not _is_inside(path, output_dir):
        return False
    return path.name.startswith(("waveform_fixture_", ".waveform_fixture_"))


if __name__ == "__main__":
    raise SystemExit(main())
