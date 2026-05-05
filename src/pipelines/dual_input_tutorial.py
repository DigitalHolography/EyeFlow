from __future__ import annotations

import h5py
import numpy as np

from input_output import (
    DOPPLER_VIEW_SCHEMA,
    HOLODOPPLER_SCHEMA,
    H5DatasetSpec,
    H5SourceSchema,
    JsonConfigValueSpec,
)
from pipelines.utils.input_access import (
    read_nested_int_setting,
    resolve_holodoppler_timing,
    resolve_required_source_array,
)

from pipeline_engine import ProcessPipeline, ProcessResult, registerPipeline, with_attrs


@registerPipeline(
    name="dual_input_tutorial",
    dag_produces=["dual_input_summary"],
)
class DualInputTutorial(ProcessPipeline):
    """
    Tutorial pipeline showing how to consume HD and DV inputs at the same time.

    This example deliberately uses the schema objects instead of hard-coded paths.
    Those schemas are Pydantic models, so dataset paths and config keys are
    declared once, normalized at construction time, and reused by pipelines.

    When launched from the EyeFlow UI runtime, the incoming `h5file` object exposes:
    - `h5file.hd`: the holodoppler input handle
    - `h5file.dv`: the doppler vision input handle
    - `h5file.work`: the current EyeFlow output/work file
    - `h5file.ef`: current root-level EyeFlow outputs, searchable via `h5file.ef["path"]`
    - `h5file.hd_config`: Holodoppler sidecar JSON values, when present
    - `h5file.dv_config`: DopplerView sidecar JSON values, when present
    """

    description = "Tutorial: read HD and DV inputs simultaneously in one pipeline."
    input_slot = "both"

    HD_EXAMPLE_DATASET_KEY = "moment0"
    DV_EXAMPLE_DATASET_KEY = "retinal_artery_mask"
    DV_EXAMPLE_CONFIG_KEY = "local_background_dist"

    @staticmethod
    def _read_schema_dataset(
        source: h5py.File,
        *,
        schema: H5SourceSchema,
        dataset_key: str,
    ) -> tuple[H5DatasetSpec, str, np.ndarray]:
        spec = schema.dataset(dataset_key)
        resolved = resolve_required_source_array(
            source,
            source_name=schema.label,
            logical_name=spec.key,
            path=spec.path,
        )
        return spec, resolved.path, resolved.value

    @staticmethod
    def _numeric_dataset_summary(
        values: np.ndarray,
    ) -> tuple[np.float32, np.int32]:
        if values.size == 0:
            return np.float32(np.nan), np.int32(0)

        try:
            numeric = np.asarray(values, dtype=float)
        except (TypeError, ValueError):
            return np.float32(np.nan), np.int32(values.size)

        return np.float32(np.nanmean(numeric)), np.int32(numeric.size)

    @staticmethod
    def _read_dv_int_config(
        h5file,
        *,
        key: str,
    ) -> tuple[JsonConfigValueSpec, np.int32]:
        spec = DOPPLER_VIEW_SCHEMA.config_value(key)
        if spec.section is None:
            raw_value = spec.read_json_config(h5file.dv_config)
            value = spec.default if raw_value is None else raw_value
            return spec, np.int32(value)

        value = read_nested_int_setting(
            h5file.dv_config,
            spec.section,
            spec.json_key,
            default=int(spec.default or 0),
        )
        return spec, np.int32(value)

    def run(self, h5file) -> ProcessResult:
        hd_h5 = h5file.hd
        dv_h5 = h5file.dv
        if hd_h5 is None or dv_h5 is None:
            raise ValueError(
                "dual_input_tutorial requires both HD and DV inputs from the UI runtime."
            )

        timing = resolve_holodoppler_timing(h5file)
        hd_spec, hd_path, hd_values = self._read_schema_dataset(
            hd_h5,
            schema=HOLODOPPLER_SCHEMA,
            dataset_key=self.HD_EXAMPLE_DATASET_KEY,
        )
        dv_spec, dv_path, dv_values = self._read_schema_dataset(
            dv_h5,
            schema=DOPPLER_VIEW_SCHEMA,
            dataset_key=self.DV_EXAMPLE_DATASET_KEY,
        )
        dv_config_spec, local_background_dist = self._read_dv_int_config(
            h5file,
            key=self.DV_EXAMPLE_CONFIG_KEY,
        )

        hd_mean, hd_size = self._numeric_dataset_summary(hd_values)
        dv_mean, dv_size = self._numeric_dataset_summary(dv_values)
        mean_delta = (
            np.float32(hd_mean - dv_mean)
            if np.isfinite(hd_mean) and np.isfinite(dv_mean)
            else np.float32(np.nan)
        )

        hd_root_keys = sorted(str(key) for key in hd_h5.keys())
        dv_root_keys = sorted(str(key) for key in dv_h5.keys())
        shared_root_keys = sorted(set(hd_root_keys) & set(dv_root_keys))

        metrics = {
            "summary/hd_root_group_count": with_attrs(
                np.int32(len(hd_root_keys)),
                {"unit": ["count"]},
            ),
            "summary/dv_root_group_count": with_attrs(
                np.int32(len(dv_root_keys)),
                {"unit": ["count"]},
            ),
            "summary/shared_root_group_count": with_attrs(
                np.int32(len(shared_root_keys)),
                {"unit": ["count"]},
            ),
            "summary/hd_example_size": with_attrs(hd_size, {"unit": ["samples"]}),
            "summary/dv_example_size": with_attrs(dv_size, {"unit": ["samples"]}),
            "summary/hd_example_mean": with_attrs(hd_mean, {"unit": ["a.u."]}),
            "summary/dv_example_mean": with_attrs(dv_mean, {"unit": ["a.u."]}),
            "summary/hd_minus_dv_example_mean": with_attrs(
                mean_delta,
                {"unit": ["a.u."]},
            ),
            "summary/both_inputs_available": with_attrs(
                np.uint8(1),
                {"unit": ["bool"]},
            ),
            "schema/hd_sampling_freq_hz": with_attrs(
                np.float32(timing.sampling_freq),
                {"unit": ["Hz"]},
            ),
            "schema/hd_batch_stride": with_attrs(
                np.float32(timing.batch_stride),
                {"unit": ["frames"]},
            ),
            "schema/hd_dt_seconds": with_attrs(
                np.float32(timing.dt_seconds),
                {"unit": ["s"]},
            ),
            "schema/dv_local_background_dist": with_attrs(
                local_background_dist,
                {"unit": ["pixels"]},
            ),
        }
        attrs = {
            "hd_source_file": str(hd_h5.filename or ""),
            "dv_source_file": str(dv_h5.filename or ""),
            "hd_schema_label": HOLODOPPLER_SCHEMA.label,
            "dv_schema_label": DOPPLER_VIEW_SCHEMA.label,
            "hd_example_key": hd_spec.key,
            "dv_example_key": dv_spec.key,
            "dv_config_key": dv_config_spec.key,
            "hd_example_path": hd_path,
            "dv_example_path": dv_path,
            "hd_example_dtype": hd_spec.dtype or "",
            "dv_example_dtype": dv_spec.dtype or "",
            "hd_example_dims": list(hd_spec.dims),
            "dv_example_dims": list(dv_spec.dims),
            "shared_root_groups": shared_root_keys,
        }
        return ProcessResult(metrics=metrics, attrs=attrs)
