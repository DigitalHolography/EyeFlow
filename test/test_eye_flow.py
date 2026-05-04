from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from eye_flow import ProcessApp  # noqa: E402
from pipelines import PipelineDescriptor  # noqa: E402


class _Var:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value

    def set(self, value):
        self._value = value


class SingleRunTests(unittest.TestCase):
    @mock.patch("ui.run.messagebox.showwarning")
    @mock.patch("ui.run.messagebox.showerror")
    def test_run_process_uses_dag_plan_order(
        self,
        showerror,
        showwarning,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            holo_path = tmp_path / "subject.holo"
            output_dir = tmp_path / "subject" / "subject_EF"
            logs: list[str] = []
            captured: dict[str, object] = {}

            load_pipeline = SimpleNamespace(name="load_data", available=True)
            target_pipeline = SimpleNamespace(name="analyze", available=True)
            plan = SimpleNamespace(
                targets=("analyze",),
                names=("load_data", "analyze"),
                descriptors=(load_pipeline, target_pipeline),
            )

            def _run_pipelines_to_output(**kwargs):
                captured.update(kwargs)
                kwargs["output_h5_path"].parent.mkdir(parents=True, exist_ok=True)
                kwargs["output_h5_path"].write_text("result", encoding="utf-8")
                return kwargs["output_h5_path"]

            app = SimpleNamespace(
                _progress_primary_style="MinimalPrimary.Horizontal.TProgressbar",
                _progress_total_units=1.0,
                _reset_progress=lambda: None,
                _selected_holo_path=lambda: holo_path,
                _selected_target_pipeline_names=lambda: ["analyze"],
                _resolve_pipeline_plan=lambda _targets: plan,
                _validate_selected_input=lambda _holo_path: SimpleNamespace(
                    holo_path=holo_path,
                    data_dir=holo_path.parent / holo_path.stem,
                    hd_h5=tmp_path / "subject_holodoppler.h5",
                    dv_h5=tmp_path / "subject_doppler_vision.h5",
                ),
                _reset_run_log=lambda message="": logs.append(message.strip()),
                _log_run=logs.append,
                _prepare_default_output_dir=lambda _input_path: output_dir,
                _default_work_h5_name_for_input=lambda _input_path: "subject_eyeflow.h5",
                _start_progress=lambda total_units, **_kwargs: setattr(
                    app,
                    "_progress_total_units",
                    total_units,
                ),
                _run_pipelines_to_output=_run_pipelines_to_output,
                _set_progress_units=lambda *_args, **_kwargs: None,
                _set_minimal_status=lambda *_args, **_kwargs: None,
            )

            ProcessApp.run_process(app)

            self.assertEqual(
                [load_pipeline, target_pipeline],
                list(captured["pipelines"]),
            )
            self.assertEqual(("analyze",), captured["target_names"])
            self.assertTrue(
                any(
                    "[DAG] Execution order -> load_data, analyze" in line
                    for line in logs
                )
            )
            showwarning.assert_not_called()
            showerror.assert_not_called()

    @mock.patch("ui.run.messagebox.showerror")
    def test_run_process_stops_when_dag_requires_unavailable_pipeline(
        self,
        showerror,
    ) -> None:
        unavailable = SimpleNamespace(
            name="load_data",
            available=False,
            missing_deps=["scipy"],
            requires=[],
        )
        plan = SimpleNamespace(
            targets=("analyze",),
            names=("load_data", "analyze"),
            descriptors=(unavailable,),
        )
        app = SimpleNamespace(
            _reset_progress=lambda: None,
            _selected_holo_path=lambda: Path("subject.holo"),
            _selected_target_pipeline_names=lambda: ["analyze"],
            _resolve_pipeline_plan=lambda _targets: plan,
        )

        ProcessApp.run_process(app)

        self.assertEqual("Pipeline unavailable", showerror.call_args.args[0])


class InputHandlingTests(unittest.TestCase):
    def _make_bare_app(self):
        app = ProcessApp.__new__(ProcessApp)
        app.holo_input_var = _Var("")
        app._selected_holo_input_paths = []
        app._synchronizing_holo_input_var = False
        app._apply_input_defaults = lambda _path: None
        logs: list[str] = []
        app._log_run = logs.append
        return app, logs

    def test_handle_dropped_paths_assigns_single_holo_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            holo_path = tmp_path / "subject.holo"
            holo_path.write_text("dummy", encoding="utf-8")

            app, logs = self._make_bare_app()

            accepted = ProcessApp._handle_dropped_paths(app, [holo_path])

            self.assertTrue(accepted)
            self.assertEqual(str(holo_path), app.holo_input_var.get())
            self.assertEqual([holo_path], app._selected_holo_input_paths)
            self.assertEqual(1, len(logs))
            self.assertIn("Drag and drop HOLO", logs[0])

    @mock.patch("ui.run.messagebox.showwarning")
    def test_validate_selected_input_requires_holo_selection(
        self,
        showwarning,
    ) -> None:
        app = ProcessApp.__new__(ProcessApp)

        resolved = ProcessApp._validate_selected_input(app, None)

        self.assertIsNone(resolved)
        self.assertEqual("Missing input", showwarning.call_args.args[0])

    def test_default_work_h5_name_uses_current_input(self) -> None:
        app = ProcessApp.__new__(ProcessApp)

        output_name = ProcessApp._default_work_h5_name_for_input(
            app,
            Path("subject.holo"),
        )

        self.assertEqual("subject_eyeflow.h5", output_name)


class PipelineTargetTests(unittest.TestCase):
    def test_selected_target_pipeline_names_ignores_unavailable_pipelines(self) -> None:
        app = ProcessApp.__new__(ProcessApp)
        app.pipeline_rows = [
            SimpleNamespace(name="A", available=True),
            SimpleNamespace(name="B", available=False),
            SimpleNamespace(name="C", available=True),
        ]
        app.pipeline_visibility = {"A": True, "B": True, "C": False}

        self.assertEqual(["A"], ProcessApp._selected_target_pipeline_names(app))

    def test_resolve_pipeline_plan_adds_required_upstream_pipeline(self) -> None:
        app = ProcessApp.__new__(ProcessApp)
        app.pipeline_dag = None
        app.pipeline_rows = [
            PipelineDescriptor(
                name="prepare_velocity",
                description="",
                available=True,
                dag_produces=("velocity_per_beat",),
            ),
            PipelineDescriptor(
                name="waveform_metrics",
                description="",
                available=True,
                dag_requires=("velocity_per_beat",),
            ),
        ]

        plan = ProcessApp._resolve_pipeline_plan(app, ["waveform_metrics"])

        self.assertEqual(
            ("prepare_velocity", "waveform_metrics"),
            plan.names,
        )

    def test_set_pipeline_visibility_persists_and_updates_summary(self) -> None:
        app = ProcessApp.__new__(ProcessApp)
        pipeline = SimpleNamespace(name="B", available=True)
        app.pipeline_catalog = {"B": pipeline}
        app.pipeline_visibility = {"A": True, "B": False}
        app._persist_pipeline_visibility = mock.Mock()
        app._update_pipeline_library_summary = mock.Mock()

        ProcessApp._set_pipeline_visibility(app, "B", True)

        self.assertTrue(app.pipeline_visibility["B"])
        app._persist_pipeline_visibility.assert_called_once_with()
        app._update_pipeline_library_summary.assert_called_once_with()


class MouseWheelBindingTests(unittest.TestCase):
    def test_mousewheel_scroll_units_handles_delta_and_button_events(self) -> None:
        self.assertEqual(
            -1, ProcessApp._mousewheel_scroll_units(SimpleNamespace(delta=120))
        )
        self.assertEqual(
            1, ProcessApp._mousewheel_scroll_units(SimpleNamespace(delta=-120))
        )
        self.assertEqual(
            -2, ProcessApp._mousewheel_scroll_units(SimpleNamespace(delta=240))
        )
        self.assertEqual(
            -1, ProcessApp._mousewheel_scroll_units(SimpleNamespace(delta=1))
        )
        self.assertEqual(
            -1, ProcessApp._mousewheel_scroll_units(SimpleNamespace(delta=0, num=4))
        )
        self.assertEqual(
            1, ProcessApp._mousewheel_scroll_units(SimpleNamespace(delta=0, num=5))
        )

    def test_bind_vertical_mousewheel_registers_handlers_that_scroll_canvas(self):
        app = ProcessApp.__new__(ProcessApp)
        widget = mock.Mock()
        canvas = mock.Mock()

        ProcessApp._bind_vertical_mousewheel(app, widget, canvas)

        self.assertEqual(
            ["<MouseWheel>", "<Button-4>", "<Button-5>"],
            [call.args[0] for call in widget.bind.call_args_list],
        )
        self.assertTrue(
            all(call.kwargs.get("add") == "+" for call in widget.bind.call_args_list)
        )

        mousewheel_handler = widget.bind.call_args_list[0].args[1]
        result = mousewheel_handler(SimpleNamespace(delta=-120))

        canvas.yview_scroll.assert_called_once_with(1, "units")
        self.assertEqual("break", result)


if __name__ == "__main__":
    unittest.main()
