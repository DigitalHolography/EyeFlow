"""Benchmark robust beat handling on local, untracked waveform fixtures.

The fixture files and optional expectation manifest are inputs only.  This tool
does not copy them into the repository or include case-specific identifiers in
source control.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py
import numpy as np

from calculations.blood_flow_velocity import (
    BeatQualitySettings,
    PerBeatAnalysisInput,
    run_per_beat_analysis,
)
from calculations.blood_flow_velocity.analysis_preparation.beat_detection import (
    find_systole_index,
)
from pipelines.waveform_shape_metrics.metrics.calculator import (
    WaveformShapeMetricsCalculator,
)
from pipelines.waveform_shape_metrics.metrics.runner import (
    run_waveform_shape_metric_calculations,
)
from pipelines.waveform_shape_metrics.velocity.constants import (
    LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
)
from pipelines.waveform_shape_metrics.velocity.outputs import (
    pack_velocity_per_beat_outputs,
)

REQUIRED_DATASETS = (
    "beats/boundary_indices",
    "beats/dt_seconds",
    "waveforms/artery/raw",
    "waveforms/vein/raw",
)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    expectations = _load_expectations(args.expectations)
    fixture_paths = sorted(args.fixtures.glob("*.h5"))
    if not fixture_paths:
        raise ValueError(f"No .h5 waveform fixtures found in {args.fixtures}.")

    cases = [_benchmark_fixture(path) for path in fixture_paths]
    failures = _check_expectations(cases, expectations)
    payload = {
        "fixture_count": len(cases),
        "summary": _aggregate(cases),
        "cases": cases,
        "expectation_failures": failures,
    }
    if args.json_output is not None:
        args.json_output.write_text(
            json.dumps(payload, indent=2, allow_nan=False),
            encoding="utf-8",
        )
    _print_summary(payload)
    return 1 if failures else 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixtures", type=Path, help="Directory of local fixture H5 files.")
    parser.add_argument(
        "--expectations",
        type=Path,
        help="Optional local JSON manifest keyed by fixture filename.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        help="Optional path for the full benchmark result.",
    )
    return parser.parse_args(argv)


def _load_expectations(path: Path | None) -> dict[str, dict[str, object]]:
    if path is None:
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    cases = payload.get("cases", payload)
    if not isinstance(cases, dict):
        raise ValueError("Expectation manifest must contain an object keyed by filename.")
    return cases


def _benchmark_fixture(path: Path) -> dict[str, object]:
    artery, vein, stored_boundaries, dt_seconds = _read_fixture(path)
    detection = find_systole_index(artery, dt=np.float32(dt_seconds))
    legacy = _run_analysis(
        artery,
        vein,
        stored_boundaries,
        dt_seconds,
    )
    robust = _run_analysis(
        artery,
        vein,
        detection.systole_indexes,
        dt_seconds,
        period_ratio=detection.interval_period_ratio,
        recovered_boundaries=detection.recovered_systole_indexes,
        nominal_period_samples=float(detection.nominal_period_samples),
        quality_settings=BeatQualitySettings(),
    )

    legacy_metrics = _metric_medians(legacy)
    robust_metrics = _metric_medians(robust)
    metric_deltas = {
        key: _finite_or_none(robust_metrics[key] - legacy_metrics[key])
        for key in sorted(legacy_metrics.keys() & robust_metrics.keys())
    }
    quality = robust.quality
    return {
        "fixture": path.name,
        "stored_boundary_count": int(stored_boundaries.size),
        "candidate_boundary_count": int(detection.systole_indexes.size),
        "boundaries_unchanged": bool(
            np.array_equal(stored_boundaries, detection.systole_indexes)
        ),
        "recovered_boundary_count": int(detection.recovered_systole_indexes.size),
        "candidate_beat_count": int(detection.systole_indexes.size - 1),
        "accepted_beat_count": int(np.sum(quality.accepted_mask)),
        "rejected_beat_count": int(np.sum(~quality.accepted_mask)),
        "rejection_flags": [int(value) for value in quality.rejection_flags],
        "legacy_heart_rate_bpm": _heart_rate(stored_boundaries, dt_seconds),
        "robust_heart_rate_bpm": _heart_rate(
            detection.systole_indexes,
            dt_seconds,
            interval_mask=quality.accepted_mask,
        ),
        "legacy_period_cv": _period_cv(stored_boundaries),
        "robust_period_cv": _period_cv(
            detection.systole_indexes,
            interval_mask=quality.accepted_mask,
        ),
        "maximum_period_ratio": _finite_or_none(
            float(np.max(detection.interval_period_ratio))
        ),
        "metric_median_deltas": metric_deltas,
    }


def _read_fixture(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    with h5py.File(path, "r") as fixture:
        if fixture.attrs.get("fixture_schema") != "eyeflow-waveform-fixture":
            raise ValueError(f"{path.name}: unsupported or missing fixture schema.")
        if int(fixture.attrs.get("fixture_schema_version", -1)) != 1:
            raise ValueError(f"{path.name}: unsupported fixture schema version.")
        missing = [dataset for dataset in REQUIRED_DATASETS if dataset not in fixture]
        if missing:
            raise ValueError(f"{path.name}: missing datasets {missing}.")
        artery = np.asarray(fixture["waveforms/artery/raw"], dtype=np.float32).reshape(-1)
        vein = np.asarray(fixture["waveforms/vein/raw"], dtype=np.float32).reshape(-1)
        boundaries = np.asarray(
            fixture["beats/boundary_indices"],
            dtype=np.int32,
        ).reshape(-1)
        dt_seconds = float(np.asarray(fixture["beats/dt_seconds"]).reshape(-1)[0])
    return artery, vein, boundaries, dt_seconds


def _run_analysis(
    artery: np.ndarray,
    vein: np.ndarray,
    boundaries: np.ndarray,
    dt_seconds: float,
    *,
    period_ratio: np.ndarray | None = None,
    recovered_boundaries: np.ndarray | None = None,
    nominal_period_samples: float | None = None,
    quality_settings: BeatQualitySettings | None = None,
):
    return run_per_beat_analysis(
        PerBeatAnalysisInput(
            arterial_velocity_signal=artery,
            venous_velocity_signal=vein,
            systolic_acceleration_peak_indexes=boundaries,
            band_limited_signal_harmonic_count=(
                LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT
            ),
            dt_seconds=dt_seconds,
            index_base=0,
            interval_period_ratio=period_ratio,
            recovered_boundary_indexes=recovered_boundaries,
            nominal_period_samples=nominal_period_samples,
            quality_settings=quality_settings,
        )
    )


def _metric_medians(result) -> dict[str, float]:
    metrics = run_waveform_shape_metric_calculations(
        pack_velocity_per_beat_outputs(result)
    )
    metric_names = set(WaveformShapeMetricsCalculator._metric_names())
    return {
        key: float(np.nanmedian(np.asarray(value, dtype=np.float64)))
        for key, value in metrics.items()
        if "/global/" in key and key.rsplit("/", maxsplit=1)[-1] in metric_names
    }


def _heart_rate(
    boundaries: np.ndarray,
    dt_seconds: float,
    *,
    interval_mask: np.ndarray | None = None,
) -> float | None:
    periods = np.diff(boundaries).astype(np.float64) * float(dt_seconds)
    if interval_mask is not None:
        periods = periods[np.asarray(interval_mask, dtype=np.bool_)]
    if periods.size == 0:
        return None
    return _finite_or_none(60.0 / float(np.mean(periods)))


def _period_cv(
    boundaries: np.ndarray,
    *,
    interval_mask: np.ndarray | None = None,
) -> float | None:
    periods = np.diff(boundaries).astype(np.float64)
    if interval_mask is not None:
        periods = periods[np.asarray(interval_mask, dtype=np.bool_)]
    if periods.size == 0 or np.mean(periods) == 0:
        return None
    return _finite_or_none(float(np.std(periods) / np.mean(periods)))


def _aggregate(cases: list[dict[str, object]]) -> dict[str, int]:
    return {
        "recovered_boundaries": sum(
            int(case["recovered_boundary_count"]) for case in cases
        ),
        "rejected_beats": sum(int(case["rejected_beat_count"]) for case in cases),
        "unchanged_boundary_cases": sum(
            bool(case["boundaries_unchanged"]) for case in cases
        ),
    }


def _check_expectations(
    cases: list[dict[str, object]],
    expectations: dict[str, dict[str, object]],
) -> list[str]:
    failures: list[str] = []
    by_name = {str(case["fixture"]): case for case in cases}
    for filename, expected in expectations.items():
        if filename not in by_name:
            failures.append(f"{filename}: fixture was not found")
            continue
        case = by_name[filename]
        checks = (
            ("expected_recovered_min", "recovered_boundary_count", lambda a, b: a >= b),
            ("expected_rejected_min", "rejected_beat_count", lambda a, b: a >= b),
            ("expected_recovered_count", "recovered_boundary_count", lambda a, b: a == b),
            ("expected_rejected_count", "rejected_beat_count", lambda a, b: a == b),
            ("require_boundary_invariance", "boundaries_unchanged", lambda a, b: a == b),
        )
        for expectation_key, result_key, predicate in checks:
            if expectation_key not in expected:
                continue
            if not predicate(case[result_key], expected[expectation_key]):
                failures.append(
                    f"{filename}: {result_key}={case[result_key]!r}, "
                    f"expected {expected[expectation_key]!r}"
                )
    return failures


def _finite_or_none(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


def _print_summary(payload: dict[str, object]) -> None:
    summary = payload["summary"]
    print(
        "fixtures={fixture_count} recovered={recovered_boundaries} "
        "rejected={rejected_beats} unchanged={unchanged_boundary_cases}".format(
            fixture_count=payload["fixture_count"],
            **summary,
        )
    )
    for case in payload["cases"]:
        print(
            "{fixture}: boundaries {stored_boundary_count}->{candidate_boundary_count}, "
            "recovered={recovered_boundary_count}, accepted={accepted_beat_count}/"
            "{candidate_beat_count}, rejected={rejected_beat_count}, "
            "HR={legacy_heart_rate_bpm}->{robust_heart_rate_bpm}".format(**case)
        )
    for failure in payload["expectation_failures"]:
        print(f"EXPECTATION FAILED: {failure}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
