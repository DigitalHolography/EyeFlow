"""Tests for selectable velocity-analysis product groups."""

from __future__ import annotations

import unittest
from contextlib import nullcontext
from types import SimpleNamespace
from unittest.mock import Mock, patch

import h5py
import numpy as np
import pytest

from input_output.readers.h5 import RawH5SourceReader
from input_output.schema import HolodopplerSource
from pipeline_engine.planning.base import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.velocity_analysis import runner as analysis_builder
from pipelines.velocity_analysis.artifacts.figures import export, flow_maps, frequency_maps


class _State:
    def __init__(self, values=None):
        self.values = dict(values or {})

    def get(self, key, default=None):
        return self.values.get(key, default)

    def set(self, key, value):
        self.values[key] = value


def _unpublished_output():
    return SimpleNamespace(available=False,
        h5=SimpleNamespace(build_arrays=lambda: nullcontext(Mock())))


def _context(options, state_values=None, scheduled=None, targeted=None):
    scheduled = set(
        scheduled
        or {
            "velocity",
            "velocity_analysis",
            "waveform_shape_metrics",
        }
    )
    return SimpleNamespace(
        state=_State(state_values),
        options_for=lambda pipeline: frozenset(options.get(pipeline, ())),
        pipeline_scheduled=lambda pipeline: pipeline in scheduled,
        pipeline_targeted=lambda pipeline: pipeline in (
            scheduled if targeted is None else set(targeted)
        ),
        option_enabled=lambda name, pipeline=None: name
        in options.get(pipeline or "velocity_analysis", ()),
    )


@pytest.mark.parametrize("source_name", ["moment0ff", "M0FF", None])
def test_pulse_velocity_overlays_use_the_flat_field_mean(tmp_path, monkeypatch, source_name):
    shape = (2, 4)
    expected_mean = np.arange(8, dtype=np.float32).reshape(shape) * 3
    volume = np.stack([expected_mean + offset for offset in range(-4, 5)])
    volume[4, 0, 1] = np.nan
    vessel_mask = np.zeros(shape, dtype=bool)
    vessel_mask[1, 1] = True
    vein_mask = np.zeros(shape, dtype=bool)
    background_average = expected_mean[:, ::-1].copy()
    frequency_average = expected_mean[::-1].copy() * 1000
    velocity_values = np.ones(shape, dtype=np.float32)
    timing = SimpleNamespace(dt_seconds=0.1)
    source = SimpleNamespace(
        holodoppler=SimpleNamespace(timing=timing),
        segmentation=SimpleNamespace(vessels=SimpleNamespace(artery=vessel_mask, vein=vein_mask)),
    )
    velocity = SimpleNamespace(
        continuous=lambda *args, **kwargs: np.ones(volume.shape[0]),
        maps=SimpleNamespace(
            background_average=background_average,
            frms_average=frequency_average,
            frms_background_average=frequency_average,
            delta_frms_average=frequency_average,
            velocity=None,
            velocity_average_masked=velocity_values,
            section_mask=np.ones(shape, dtype=bool),
        ),
    )
    analysis = SimpleNamespace(source_data=SimpleNamespace(source=source), velocity=velocity, per_beat_result=None)
    images = {}

    def record_overlay(writer, fig, suffix, **kwargs):
        background, overlay = fig.axes[0].images
        images[suffix] = np.array(background.get_array())
        if suffix.startswith("artery"):
            np.testing.assert_array_equal(overlay.get_array()[..., 3] > 0, vessel_mask)
        export._plt().close(fig)
        return [tmp_path / suffix]

    monkeypatch.setattr(export, "EXPORTERS", (frequency_maps._export_frequency_maps, flow_maps._export_flow_maps))
    monkeypatch.setattr(export, "export_cross_section_profile_artifacts", lambda *_: [])
    monkeypatch.setattr(export.FigureArtifactWriter, "path", lambda self, suffix: tmp_path / suffix)
    def record_image(writer, image, suffix):
        images[suffix] = image
        return tmp_path / suffix

    monkeypatch.setattr(export.FigureArtifactWriter, "save_image", record_image)
    monkeypatch.setattr(export.FigureArtifactWriter, "save_figure_pair", record_overlay)
    monkeypatch.setattr(frequency_maps, "_heatmap_with_colorbar", lambda *a, **kw: [])
    monkeypatch.setattr(frequency_maps, "_image_map", lambda *a, **kw: tmp_path / "frequency.png")
    monkeypatch.setattr(flow_maps, "_image_map", lambda *a, **kw: tmp_path / "flow.png")
    monkeypatch.setattr(flow_maps, "_colorbar", lambda *a, **kw: tmp_path / "colorbar.png")
    monkeypatch.setattr(flow_maps, "_histogram_plot", lambda *a: [])
    monkeypatch.setattr(flow_maps, "_combined_plot", lambda *a: tmp_path / "combined.png")

    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        if source_name is not None:
            hd.create_dataset(source_name, data=volume)
        adapter = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
        ctx = SimpleNamespace(
            inputs=SimpleNamespace(hd=SimpleNamespace(as_holodoppler=lambda: adapter)),
            output=SimpleNamespace(available=True),
            dag=SimpleNamespace(require=lambda name: analysis),
            state=_State(),
        )
        analysis_builder._pulse_figures(ctx)

    flow_background = expected_mean if source_name else background_average
    separate_background = expected_mean if source_name else frequency_average / 1000
    outside = images["v_map.png"][~vessel_mask]
    np.testing.assert_allclose(outside, np.repeat((flow_background[~vessel_mask] / 21)[:, None], 3, axis=1))
    assert np.ptp(images["v_map.png"][vessel_mask]) > 0
    for vessel in ("artery", "vein"):
        np.testing.assert_allclose(images[f"{vessel}_seg_map_bkg.png"], separate_background)


@pytest.mark.parametrize("scheduled", [False, True])
@pytest.mark.parametrize(
    ("options", "expected"),
    [
        ((), False),
        (("segments", "quadrants"), False),
        (("segment_velocity_maps",), False),
        (("velocity_profiles",), True),
        (("velocity_profile_analysis",), True),
        (("velocity_profile_fft",), True),
    ],
)
def test_longitudinal_velocity_profiles_require_a_profile_consumer(
    options, expected, scheduled,
):
    source = SimpleNamespace(
        source=SimpleNamespace(
            segmentation=SimpleNamespace(
                vessels=SimpleNamespace(artery="artery_mask", vein="vein_mask"),
                optic_disc="optic_disc",
            ),
        ),
        profile_settings="settings",
        provenance={"beat_index_base": 0},
    )
    ctx = _context(
        {"velocity_analysis": options},
        scheduled={"velocity_analysis"} if scheduled else {"waveform_shape_metrics"},
    )
    ctx.output = _unpublished_output()
    with patch.object(
        analysis_builder,
        "analyze_velocity_segment_profiles",
        return_value={"artery": "artery", "vein": "vein"},
    ) as analyze:
        analysis_builder._segment_velocity_inputs(
            "velocity_map", source, "rings", ctx, cycle_boundary_indexes=(0, 5, 10),
        )
    assert analyze.call_args.kwargs["compute_longitudinal_profiles"] is (
        expected and scheduled
    )
    assert analyze.call_args.kwargs["compute_transverse_profiles"] is (
        expected and scheduled
    )


@pytest.mark.parametrize(
    ("blood_volume_rate_options", "expected"),
    [
        ((), False),
        (("masked_edges",), False),
        (("gradient_edges",), True),
        (("gradient_edges", "masked_edges"), True),
    ],
)
def test_only_gradient_edge_flow_requires_transverse_velocity_profiles(
    blood_volume_rate_options, expected,
):
    source = SimpleNamespace(
        source=SimpleNamespace(
            segmentation=SimpleNamespace(
                vessels=SimpleNamespace(artery="artery_mask", vein="vein_mask"),
                optic_disc="optic_disc",
            ),
        ),
        profile_settings="settings",
        provenance={"beat_index_base": 0},
    )
    ctx = _context(
        {
            "velocity_analysis": ("segments",),
            "blood_volume_rate": blood_volume_rate_options,
        },
        scheduled={"velocity_analysis", "blood_volume_rate"},
    )
    ctx.output = _unpublished_output()
    with patch.object(
        analysis_builder,
        "analyze_velocity_segment_profiles",
        return_value={"artery": "artery", "vein": "vein"},
    ) as analyze:
        analysis_builder._segment_velocity_inputs(
            "velocity_map", source, "rings", ctx, cycle_boundary_indexes=(0, 5, 10),
        )
    assert analyze.call_args.kwargs["compute_transverse_profiles"] is expected
    assert analyze.call_args.kwargs["compute_longitudinal_profiles"] is False


class VelocityAnalysisPipelineTests(unittest.TestCase):
    def test_profile_options_are_disabled_by_default(self) -> None:
        load_pipeline_catalog()
        options = {
            option.name: option
            for option in PIPELINE_REGISTRY["velocity_analysis"].options
        }
        profiles = options["velocity_profiles"]
        self.assertEqual(("segments",), options["segment_velocity_maps"].requires)
        self.assertEqual(("segments",), options["velocity_profiles"].requires)
        self.assertEqual(("segments",), options["quadrants"].requires)
        analysis = options["velocity_profile_analysis"]
        fft = options["velocity_profile_fft"]
        self.assertFalse(profiles.default_enabled)
        self.assertFalse(analysis.default_enabled)
        self.assertFalse(fft.default_enabled)
        self.assertEqual(("velocity_profiles",), analysis.requires)
        self.assertEqual(("velocity_profiles",), fft.requires)
        descriptor = PIPELINE_REGISTRY["velocity_analysis"]
        for selected, present, absent in (
            ((), (), ("segments", "maps", "profiles", "profile_fits", "profile_fft", "quadrants")),
            (("segments",), ("segments",), ("maps", "profiles", "profile_fft")),
            (("velocity_profiles",), ("segments", "profiles"), ("maps", "profile_fits", "profile_fft")),
            (("velocity_profile_analysis",), ("segments", "profiles", "profile_fits"), ("maps", "profile_fft")),
            (("velocity_profile_fft",), ("segments", "profiles", "profile_fft"), ("maps", "profile_fits")),
            (("segment_velocity_maps",), ("segments", "maps", "map_videos"), ("profiles", "profile_fft")),
        ):
            stages = descriptor.build_graph(selected).plan().order
            self.assertTrue(set(present).issubset(stages))
            self.assertTrue(set(absent).isdisjoint(stages))


    def test_core_segment_requirement_uses_synchronized_segment_selection(self) -> None:
        ctx = _context(
            {
                "velocity_analysis": ("segments",),
                "waveform_shape_metrics": ("segments",),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": ("segment_velocity_maps",),
                "waveform_shape_metrics": (),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": ("segments",),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": ("quadrants",),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": (),
            }
        )

        self.assertFalse(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": ("segments",),
            },
            scheduled={"velocity", "velocity_analysis"},
        )

        self.assertFalse(analysis_builder._segments_required(ctx))


    def test_core_skips_segment_extraction_when_not_required(self) -> None:
        analysis = SimpleNamespace(
            cycle_boundary_indexes=np.asarray([0, 1], dtype=np.int32),
            cardiac_cycle=SimpleNamespace(spectral="cardiac_cycle"),
            continuous=lambda vessel, raw=False: np.asarray(
                [1.0, 2.0], dtype=np.float32
            ),
        )
        source = SimpleNamespace(
            timing=SimpleNamespace(dt_seconds=0.1),
            provenance={"beat_index_base": 0},
        )
        ctx = SimpleNamespace()

        with patch.object(analysis_builder, "_segment_velocity_inputs") as extract:
            artery, vein = analysis_builder._build_segments(
                ctx,
                analysis,
                source,
                segments_required=False,
            )

        extract.assert_not_called()
        self.assertIsNone(artery)
        self.assertIsNone(vein)

    def test_core_plain_profiles_do_not_stream_fft_or_retain_velocity_maps(self) -> None:
        source = SimpleNamespace(
            source=SimpleNamespace(
                segmentation=SimpleNamespace(
                    vessels=SimpleNamespace(
                        artery="artery_mask",
                        vein="vein_mask",
                    ),
                    optic_disc="optic_disc",
                ),
            ),
            profile_settings="settings",
            provenance={"beat_index_base": 1},
        )
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name == "velocity_analysis",
            option_enabled=lambda name, pipeline=None: name == "velocity_profiles",
            inputs=SimpleNamespace(
                hd=SimpleNamespace(filename="hd.h5"),
                dv=SimpleNamespace(filename="dv.h5"),
            ),
            state=SimpleNamespace(raw={}),
            output=_unpublished_output(),
        )

        with patch.object(
            analysis_builder,
            "analyze_velocity_segment_profiles",
            return_value={"artery": "artery", "vein": "vein"},
        ) as analyze:
            artery, vein = analysis_builder._segment_velocity_inputs(
                "velocity_map",
                source,
                "rings",
                ctx,
                cycle_boundary_indexes=(1, 6, 11),
            )

        self.assertEqual(("artery", "vein"), (artery, vein))
        self.assertFalse(analyze.call_args.kwargs["retain_velocity_maps"])
        self.assertFalse(analyze.call_args.kwargs["velocity_profile_fft"])
        self.assertEqual(
            (1, 6, 11),
            analyze.call_args.kwargs["cycle_boundary_indexes"],
        )
        self.assertEqual(1, analyze.call_args.kwargs["index_base"])

    def test_core_explicit_fft_option_streams_fft_without_retaining_maps(self) -> None:
        source = SimpleNamespace(
            source=SimpleNamespace(
                segmentation=SimpleNamespace(
                    vessels=SimpleNamespace(
                        artery="artery_mask",
                        vein="vein_mask",
                    ),
                    optic_disc="optic_disc",
                ),
            ),
            profile_settings="settings",
            provenance={"beat_index_base": 0},
        )
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name == "velocity_analysis",
            option_enabled=lambda name, pipeline=None: name
            in {"velocity_profiles", "velocity_profile_fft"},
            inputs=SimpleNamespace(
                hd=SimpleNamespace(filename="hd.h5"),
                dv=SimpleNamespace(filename="dv.h5"),
            ),
            state=SimpleNamespace(raw={}),
            output=_unpublished_output(),
        )
        with patch.object(
            analysis_builder,
            "analyze_velocity_segment_profiles",
            return_value={"artery": "artery", "vein": "vein"},
        ) as analyze:
            analysis_builder._segment_velocity_inputs(
                "velocity_map",
                source,
                "rings",
                ctx,
                cycle_boundary_indexes=(0, 5, 10),
            )
        self.assertTrue(analyze.call_args.kwargs["velocity_profile_fft"])
        self.assertFalse(analyze.call_args.kwargs["retain_velocity_maps"])


    def test_pdf_report_requires_pulse_pngs(self) -> None:
        ctx = _context(
            {"velocity_analysis": (), "waveform_shape_metrics": ()},
            scheduled={
                "velocity",
                "velocity_analysis",
                "waveform_shape_metrics",
                "pdf_report",
            },
        )

        self.assertTrue(analysis_builder._pulse_pngs_required(ctx))

    def test_downstream_dependency_does_not_require_pulse_pngs(self) -> None:
        ctx = _context(
            {"velocity_analysis": ()},
            scheduled={"velocity", "velocity_analysis", "blood_volume_rate"},
            targeted={"blood_volume_rate"},
        )

        self.assertFalse(analysis_builder._pulse_pngs_required(ctx))


if __name__ == "__main__":
    unittest.main()
