"""Tests for runtime pipeline discovery paths."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from app_settings import PIPELINES_DIR_ENV, app_version_dir_name  # noqa: E402
import pipelines as pipeline_package  # noqa: E402


class PipelineDiscoveryTests(unittest.TestCase):
    def tearDown(self) -> None:
        sys.modules.pop("pipelines.runtime_sample", None)
        pipeline_package.load_pipeline_catalog()

    def test_version_dir_name_prefixes_plain_versions(self) -> None:
        self.assertEqual("v0.11.0", app_version_dir_name("0.11.0"))
        self.assertEqual("v0.11.0", app_version_dir_name("v0.11.0"))

    def test_load_catalog_includes_runtime_pipeline_folder(self) -> None:
        original_env = os.environ.get(PIPELINES_DIR_ENV)
        original_package_path = list(pipeline_package.__path__)

        with tempfile.TemporaryDirectory() as temp_dir:
            runtime_pipelines = Path(temp_dir) / "pipelines"
            sample_pipeline = runtime_pipelines / "runtime_sample"
            sample_pipeline.mkdir(parents=True)
            (sample_pipeline / "__init__.py").write_text(
                "\n".join(
                    [
                        "from pipeline_engine.imports import pipeline",
                        "",
                        '@pipeline(name="runtime_sample", description="Runtime sample")',
                        "def run(ctx):",
                        "    return None",
                    ]
                ),
                encoding="utf-8",
            )

            os.environ[PIPELINES_DIR_ENV] = str(runtime_pipelines)
            try:
                available, missing = pipeline_package.load_pipeline_catalog()
            finally:
                if original_env is None:
                    os.environ.pop(PIPELINES_DIR_ENV, None)
                else:
                    os.environ[PIPELINES_DIR_ENV] = original_env
                pipeline_package.__path__[:] = original_package_path

        names = {descriptor.name for descriptor in [*available, *missing]}
        self.assertIn("runtime_sample", names)


if __name__ == "__main__":
    unittest.main()
