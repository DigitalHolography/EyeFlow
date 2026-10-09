"""End-to-end frequency-band pipeline coverage with physical velocity outputs."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import h5py
import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.run_service import execute_run, resolve_run_spec
from pipelines import load_pipeline_catalog


def test_frequency_band_run_produces_physical_downstream_outputs() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        holo = _write_frequency_band_run(Path(temp_dir))
        available, missing = load_pipeline_catalog()
        assert not {
            "velocity",
            "topology_core",
            "velocity_analysis",
            "waveform_shape_metrics",
            "absolute_waveform_metrics",
            "blood_volume_rate",
        } & {descriptor.name for descriptor in missing}

        spec = resolve_run_spec(
            input_paths=[holo],
            target_names=[
                "velocity_analysis",
                "waveform_shape_metrics",
                "absolute_waveform_metrics",
                "blood_volume_rate",
            ],
            pipelines=available,
            pipeline_options={
                "velocity_analysis": ("segments",),
                "waveform_shape_metrics": (),
                "absolute_waveform_metrics": (),
                "blood_volume_rate": ("masked_edges",),
            },
            band_ratio_frequency_scale_hz=2.0,
        )

        result = execute_run(spec)

        assert result.succeeded, result.failures
        from input_output.h5_access import PipelineH5Output

        schema = EyeFlowOutputPaths.active()
        with h5py.File(result.outputs[0], "r") as file:
            assert "Processing" not in file
            assert "doppler_moments" in json.loads(
                file.attrs["execution_variant_failures"]
            )
            output = PipelineH5Output(file, processing_root="ProcessingAlt")
            assert output.attrs["velocity_estimation_method"] == "frequency_bands"
            assert output.attrs["velocity_quantity"] == "physical_velocity"
            assert output.attrs["velocity_unit"] == "mm/s"
            assert output.attrs["band_ratio_frequency_scale_hz"] == 2.0
            assert output.attrs["band_lf_zero_sample_count"] == 1
            assert output.attrs["band_lf_near_zero_sample_count"] == 1

            artery_velocity = output.get(schema.analysis.retinal_artery_velocity_signal)
            assert artery_velocity.attrs["unit"] == "mm/s"
            assert artery_velocity.attrs["velocity_estimation_method"] == "frequency_bands"
            assert artery_velocity.attrs["band_ratio_frequency_scale_hz"] == 2.0
            assert np.any(np.isfinite(artery_velocity[:]))

            artery_per_beat = output.get(schema.artery_per_beat.velocity_signal)
            assert artery_per_beat.attrs["unit"] == "mm/s"
            assert artery_per_beat.attrs["velocity_estimation_method"] == "frequency_bands"
            assert artery_per_beat.attrs["band_ratio_frequency_scale_hz"] == 2.0

            cycle_durations = output.get(schema.cardiac_cycle.systolic_cycle_duration_seconds)
            assert cycle_durations.attrs["unit"] == "s"
            assert cycle_durations.ndim == 1
            assert cycle_durations.size > 0

            frequency_map = output.get(schema.analysis.fRMS_avg)
            assert frequency_map.attrs["unit"] == "Hz"
            assert frequency_map.attrs["band_ratio_frequency_scale_hz"] == 2.0

            velocity_average_masked = output.get(schema.analysis.velocity_map_avg_masked)
            assert velocity_average_masked.attrs["unit"] == "mm/s"

            artery_flow = output.get(schema.blood_volume_rate.artery.masked_edges)
            assert artery_flow.attrs["unit"] == "mm^3/s"

            absolute_root = output.get(schema.absolute_waveform_metrics_root)
            assert _dataset_count(absolute_root) > 0


def _write_frequency_band_run(root: Path, stem: str = "band_scan") -> Path:
    frames, height, width = 160, 64, 64
    root.mkdir(parents=True, exist_ok=True)
    holo = root / f"{stem}.holo"
    holo.write_text("holo", encoding="utf-8")
    hd_path = root / stem / f"{stem}_HD" / "h5" / f"{stem}_HD_output.h5"
    dv_path = root / stem / f"{stem}_DV" / "h5" / f"{stem}_DV.h5"
    hd_path.parent.mkdir(parents=True)
    dv_path.parent.mkdir(parents=True)

    artery = np.zeros((height, width), dtype=bool)
    vein = np.zeros_like(artery)
    artery[29:32, 38:60] = True
    vein[38:60, 34:37] = True
    disc_y, disc_x = np.ogrid[:height, :width]
    optic_disc = (disc_y - 32) ** 2 + (disc_x - 32) ** 2 <= 6**2

    time = np.arange(frames, dtype=np.float32) / np.float32(20.0)
    artery_ratio = 4.0 + 0.8 * np.sin(2.0 * np.pi * 1.0 * time)
    vein_ratio = 3.0 + 0.4 * np.sin(2.0 * np.pi * 1.0 * time + 0.5)
    low = np.ones((frames, height, width), dtype=np.float32)
    high = np.ones_like(low)
    high[:, artery] = artery_ratio[:, None]
    high[:, vein] = vein_ratio[:, None]
    low[0, 0, 0] = 0.0
    low[0, 0, 1] = 1e-8

    with h5py.File(hd_path, "w") as hd:
        hd.create_dataset("band_0_3000_9000", data=low)
        hd.create_dataset("band_1_9000_18000", data=high)
        hd.create_dataset(
            "HD_parameters",
            data=json.dumps(
                {
                    "pixel_pitch": [20e-6, 20e-6],
                    "sampling_freq": 20.0,
                    "batch_stride": 1.0,
                }
            ),
        )
        hd.attrs["number_of_radii_in_FOV"] = 8

    with h5py.File(dv_path, "w") as dv:
        retina = dv.create_group("segmentation/Retina")
        retina.create_dataset("artery_mask", data=artery)
        retina.create_dataset("vein_mask", data=vein)
        disc = dv.create_group("segmentation/OpticDisc")
        disc.create_dataset("mask", data=optic_disc)
        disc.create_dataset("center", data=np.asarray([32.0, 32.0]))
        disc.create_dataset("width", data=np.float32(12.0))
        disc.create_dataset("height", data=np.float32(12.0))

    return holo


def _dataset_count(group: h5py.Group) -> int:
    count = 0

    def visitor(_name: str, value: h5py.Group | h5py.Dataset) -> None:
        nonlocal count
        if isinstance(value, h5py.Dataset):
            count += 1

    group.visititems(visitor)
    return count
