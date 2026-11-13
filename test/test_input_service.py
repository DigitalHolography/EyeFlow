"""Real HDF5 and spawn transport, readiness and exclusive buffer ownership."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from multiprocessing.shared_memory import SharedMemory
import multiprocessing as mp
from threading import Event
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from pipeline_engine.io.input.service import InputPolicy, InputService
from pipeline_engine.data.storage import ProductStorage
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.stage_graph import InputAccess, InputUse, StageGraph, StageSpec


@contextmanager
def inputs(tmp_path, values=None, *, policy=None, budget=None, location="network"):
    values = np.arange(96 * 4 * 8, dtype=np.float32).reshape(96, 4, 8) if values is None else values
    with h5py.File(tmp_path / "source.h5", "w") as source:
        source["moment0"] = values
        source["M0"] = source["moment0"]
        source["moment2"] = values + 1
        source["band_8_2000_4000"] = values + 2
        source["band_0_9000_18000"] = values + 3
        source["HD_parameters"] = '{"sampling_freq":40000,"batch_stride":1}'
        source.create_group("spectrograms")
    with h5py.File(tmp_path / "source.h5", "r") as source:
        with ProductStorage(directory=tmp_path / "scratch", budget=budget) as storage:
            with patch("pipeline_engine.io.input.service.storage_location", return_value=location):
                service = InputService({"hd": source}, storage, mp.get_context("spawn"),
                                       policy=policy or InputPolicy(lookahead=False, input_backend="replay"))
            clients = []
            def client():
                item = service.client()
                clients.append(item)
                return service.files(item)["hd"]
            try:
                yield service, client, values
            finally:
                for item in clients:
                    item.close()
                service.close()


def _never_run(ctx):
    raise AssertionError("Manifest construction executed a scientific callback")


@pytest.mark.parametrize("method,expected", [
    ("doppler_moments", {"M0", "moment2"}),
    ("bandratio", {"band_8_2000_4000", "band_0_9000_18000"}),
    ("both", {"M0", "moment2", "band_8_2000_4000", "band_0_9000_18000"}),
])
def test_manifest_methods_aliases_dynamic_sorting_and_no_bulk_reads(tmp_path, method, expected):
    from pipelines.velocity.runner import velocity_graph
    with inputs(tmp_path) as (service, client, values):
        service.configure({"velocity": velocity_graph()}, {"velocity_workflow_selection": method})
        assert {key[1] for key in service.manifest if key[1] != "HD_parameters"} == expected
        assert service.source_reads == {}
        assert service._replays == {}
        assert len(service.pairs) == len(expected)
        assert service.aliases["hd", "M0"] == service.aliases["hd", "moment0"]


def test_metadata_optional_and_legacy_spectrogram_group_do_not_read(tmp_path):
    with inputs(tmp_path) as (service, client, _):
        graph = StageGraph((StageSpec("source", _never_run, inputs=(
            InputUse("hd.moment0", InputAccess.METADATA),
            InputUse("hd.moment0ff", InputAccess.SEQUENTIAL, optional=True),
            InputUse("hd.spectrograms", InputAccess.METADATA, optional=True),
        )),))
        service.configure({"spectra": graph})
        assert service.manifest[service.aliases["hd", "moment0"]]["access"] == {InputAccess.METADATA}
        assert len(service.unresolved) == 2
        assert service.source_reads == {}
        assert not service._replays


def test_first_region_is_visible_without_full_staging_and_reused_by_alias(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        first, second = client(), client()
        actual = first["moment0"].read_slice(np.s_[2:7])
        np.testing.assert_array_equal(actual, values[2:7])
        assert not actual.flags.writeable
        assert actual.base.shape[0] == 5  # Cannot inspect unfinished data through .base.
        replay = service._replays[service.aliases["hd", "moment0"]]
        assert replay.states == ["ready", "not_requested", "not_requested"]
        calls = service.storage.metrics.snapshot()["source_reads"]
        np.testing.assert_array_equal(second["M0"][2:7], values[2:7])
        assert service.storage.metrics.snapshot()["source_reads"] == calls


@pytest.mark.parametrize("selection", [
    np.s_[31:65, 1:3, ::2], np.s_[-1], np.s_[-40:-2], np.s_[100:101],
    np.s_[70:12:-3], np.s_[..., 1], (), np.s_[None, :3],
    ([4, 30, 2], slice(None), 1), np.s_[32:32],
])
def test_indexing_and_full_fallback_are_exact(tmp_path, selection):
    with inputs(tmp_path) as (service, client, values):
        actual = client()["moment0"][selection]
        np.testing.assert_array_equal(actual, values[selection])
        if not (isinstance(selection, tuple) and selection and isinstance(selection[0], list)):
            assert not actual.flags.writeable
        else:
            assert not np.shares_memory(actual, np.asarray(client()["moment0"]))


def test_multiple_clients_wait_for_same_region_and_never_get_partial_data(tmp_path):
    entered, proceed = Event(), Event()
    with inputs(tmp_path) as (service, client, values):
        original = service._read_direct
        def held(dataset, destination, selection):
            destination[:] = -999
            entered.set()
            assert proceed.wait(5)
            original(dataset, destination, selection)
        service._read_direct = held
        a, b = client()["moment0"], client()["M0"]
        with ThreadPoolExecutor(2) as pool:
            first = pool.submit(a.read_slice, np.s_[40:43])
            assert entered.wait(5)
            second = pool.submit(b.read_slice, np.s_[40:43])
            assert not first.done() and not second.done()
            proceed.set()
            np.testing.assert_array_equal(first.result(5), values[40:43])
            np.testing.assert_array_equal(second.result(5), values[40:43])
        assert service.storage.metrics.snapshot()["source_reads"] == 1


def test_urgent_region_overtakes_speculation_at_read_boundary(tmp_path):
    entered, proceed = Event(), Event()
    with inputs(tmp_path, policy=InputPolicy(lookahead=True, input_backend="replay")) as (service, client, values):
        original = service._read_direct
        order = []
        def held(dataset, destination, selection):
            order.append(selection.start)
            if len(order) == 1:
                entered.set()
                assert proceed.wait(5)
            original(dataset, destination, selection)
        service._read_direct = held
        service.prefetch(("hd", "moment0"))
        assert entered.wait(5)
        service.prefetch(("hd", "moment0"), first=32)
        file = client()
        # Send control explicitly so arrival before releasing the slab is
        # deterministic, without relying on a thread's scheduling latency.
        file.client.send("read", "hd", "moment0", (64, 68))
        proceed.set()
        reference, error = file.client.responses.get(timeout=5)
        assert error is None
        array = reference.open()
        try:
            np.testing.assert_array_equal(array, values[64:68])
            assert order[:2] == [0, 64]
        finally:
            array._mmap.close()


def test_read_error_preserves_ready_regions_and_removes_incomplete_file(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        file = client()
        first = file["moment0"][:2]
        def fail(*args):
            raise OSError("injected read failure")
        service._read_direct = fail
        with pytest.raises(OSError, match="injected read failure"):
            file["moment0"][64:68]
        np.testing.assert_array_equal(first, values[:2])
        assert service._replays[service.aliases["hd", "moment0"]].states[-1] == "failed"
    assert not list((tmp_path / "scratch").glob("*.npy"))


def _spawn_consumer(file, results, crash=False):
    try:
        if crash:
            with file["moment0"].iter_blocks(frames=8, writable=True) as blocks:
                lease = next(blocks)
                with lease as block:
                    results.put((blocks.identity, block.shape))
                    # Parent terminates this worker while it owns the lease.
                    Event().wait(30)
        else:
            result = []
            with file["moment0"].iter_blocks(frames=13, writable=True) as blocks:
                for lease in blocks:
                    with lease as block:
                        result.append(float(block.sum(dtype=np.float64)))
                        block[:] = -1
            results.put(result)
    finally:
        file.client.close()


def test_spawn_worker_streams_exact_values_without_materializing_source(tmp_path):
    with inputs(tmp_path, policy=InputPolicy(input_backend="direct")) as (service, client, values):
        context = mp.get_context("spawn")
        results = context.Queue()
        file = client()
        process = context.Process(target=_spawn_consumer, args=(file, results))
        process.start()
        try:
            actual = results.get(timeout=15)
            expected = [float(values[start:start + 13].sum(dtype=np.float64)) for start in range(0, 96, 13)]
            assert actual == expected
            process.join(10)
            assert process.exitcode == 0
            assert not service._replays
            assert service.storage.metrics.snapshot()["input_stream_buffers_created"] == 2
            np.testing.assert_array_equal(service.sources["hd"]["moment0"][()], values)
        finally:
            if process.is_alive():
                process.terminate()
            process.join()
            process.close()
            service.release_client(file.client.identity)
            results.close()


def test_two_buffers_overlap_without_overwriting_lease_and_early_exit(tmp_path):
    second_ready = Event()
    budget = ResourceBudget(memory_bytes=64 * 2**20)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(input_backend="direct")) as (service, client, values):
        original = service._read_stream
        def measured(identity, index, *completion):
            original(identity, index, *completion)
            if service._streams[identity]["buffers"][index]["start"] == 8:
                second_ready.set()
        service._read_stream = measured
        file = client()
        with file["moment0"].iter_blocks(frames=8, writable=True) as blocks:
            lease = next(blocks)
            with lease as block:
                assert second_ready.wait(5)
                np.testing.assert_array_equal(block, values[:8])
                block[:] = -50
                with pytest.raises(RuntimeError, match="Release"):
                    next(blocks)
                assert service.storage.metrics.snapshot()["source_reads"] == 2
                names = [buf["memory"].name for buf in service._streams[blocks.identity]["buffers"]]
                assert budget.used["input"] == 2 * values[:8].nbytes
            lease = next(blocks)
            with lease as block:
                np.testing.assert_array_equal(block, values[8:16])
        file.client.rpc("dtype", "hd", "moment0")  # Barrier after stream_close.
        assert not service._streams
        assert budget.used["input"] == 0
        for name in names:
            with pytest.raises(FileNotFoundError):
                SharedMemory(name=name)


def test_worker_termination_reclaims_leased_buffers_and_budget(tmp_path):
    budget = ResourceBudget(memory_bytes=64 * 2**20)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(input_backend="direct")) as (service, client, _):
        file, barrier = client(), client()
        context = mp.get_context("spawn")
        results = context.Queue()
        process = context.Process(target=_spawn_consumer, args=(file, results, True))
        process.start()
        try:
            stream_id, shape = results.get(timeout=15)
            assert stream_id in service._streams
            process.terminate()
            process.join(10)
            service.release_client(file.client.identity)
            barrier.client.rpc("dtype", "hd", "moment0")
            assert not service._streams
            assert budget.used["input"] == 0
        finally:
            if process.is_alive():
                process.kill()
            process.join()
            process.close()
            results.close()


def test_stream_exception_cleanup_and_capacity_limit(tmp_path):
    with inputs(tmp_path, policy=InputPolicy(max_streams=1, input_backend="direct")) as (service, client, _):
        file = client()
        with pytest.raises(ValueError, match="consumer failure"):
            with file["moment0"].iter_blocks(frames=7) as blocks:
                with pytest.raises(OSError, match="capacity"):
                    file["moment2"].iter_blocks()
                with next(blocks) as block:
                    assert not block.flags.writeable
                    raise ValueError("consumer failure")
        file.client.rpc("dtype", "hd", "moment0")
        assert not service._streams


@pytest.mark.parametrize("location,mode", [("local", "direct"), ("network", "direct"), ("network", "replay")])
def test_spectral_streaming_preserves_native_gap_and_nonfinite_results(tmp_path, location, mode):
    from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms
    values = np.random.default_rng(5).normal(size=(39, 5, 4, 7)).astype(np.float32)
    values[0, 0, 0, 1] = np.nan
    values[18, 1, 2, 2] = np.inf
    mask = np.eye(4, 7, dtype=bool)
    with inputs(tmp_path, values, location=location, policy=InputPolicy(input_backend=mode)) as (service, client, original):
        expected = preprocess_spectrograms(original, mask)
        actual = preprocess_spectrograms(client()["moment0"], mask)
        assert actual.keys() == expected.keys()
        for key in expected:
            np.testing.assert_array_equal(actual[key][0], expected[key][0])
            assert actual[key][1] == expected[key][1]
        assert bool(service._replays) == (mode == "replay")


def test_paired_velocity_current_window_precedes_next_window(tmp_path):
    from pipelines.velocity.runner import velocity_graph
    finished = Event()
    with inputs(tmp_path, policy=InputPolicy(lookahead=True, input_backend="replay")) as (service, client, values):
        service.configure({"velocity": velocity_graph()}, {"velocity_workflow_selection": "doppler_moments"})
        order = []
        original = service._read_direct
        def measured(dataset, destination, selection):
            order.append((dataset.name.strip("/"), selection.start))
            original(dataset, destination, selection)
            if len(order) == 4:
                finished.set()
        service._read_direct = measured
        np.testing.assert_array_equal(client()["moment0"][:3], values[:3])
        assert finished.wait(5)
        assert order[:4] == [("M0", 0), ("moment2", 0), ("M0", 32), ("moment2", 32)]
        assert not any("band_" in name for name, _ in order)


def test_stream_and_random_consumers_reuse_progressive_backing(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        service.configure({"mixed": StageGraph((StageSpec("stage", _never_run, inputs=(
            InputUse("hd.moment0", InputAccess.RANDOM),
            InputUse("hd.moment0", InputAccess.STREAM),
        )),))})
        file = client()
        np.testing.assert_array_equal(file["moment0"][:32], values[:32])
        with file["M0"].iter_blocks(frames=16, writable=True) as blocks:
            for lease in blocks:
                with lease as block:
                    np.testing.assert_array_equal(block, values[lease.start:lease.stop])
                    block[:] = -1
        np.testing.assert_array_equal(file["moment0"][:], values)
        assert service.storage.metrics.snapshot()["source_reads"] == 3


def test_stream_read_failure_releases_buffers_and_reservations(tmp_path):
    budget = ResourceBudget(memory_bytes=64 * 2**20)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(input_backend="direct")) as (service, client, _):
        def fail(*args):
            raise OSError("injected stream failure")
        service._read_direct = fail
        file = client()
        with pytest.raises(OSError, match="injected stream failure"):
            with file["moment0"].iter_blocks() as blocks:
                next(blocks)
        file.client.rpc("dtype", "hd", "moment0")
        assert not service._streams
        assert budget.used["input"] == 0


def test_read_direct_and_repeated_threaded_requests_on_one_client(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        dataset = client()["moment0"]
        with ThreadPoolExecutor(4) as pool:
            results = list(pool.map(lambda i: dataset[i:i + 2].copy(), [1, 35, 73, 10]))
        for result, index in zip(results, [1, 35, 73, 10]):
            np.testing.assert_array_equal(result, values[index:index + 2])
        dest = np.zeros((3, 4, 8), dtype=np.float64)
        dataset.read_direct(dest, source_sel=np.s_[40:42], dest_sel=np.s_[1:])
        np.testing.assert_array_equal(dest[1:], values[40:42])


def test_shutdown_cancels_pending_full_replay_without_deadlock(tmp_path):
    entered, proceed = Event(), Event()
    with inputs(tmp_path) as (service, client, _):
        original = service._read_direct
        def held(*args):
            entered.set()
            assert proceed.wait(5)
            original(*args)
        service._read_direct = held
        file = client()
        with ThreadPoolExecutor(2) as pool:
            read = pool.submit(np.asarray, file["moment0"])
            assert entered.wait(5)
            close = pool.submit(service.close)
            proceed.set()
            close.result(5)
            with pytest.raises(OSError, match="cancelled|stopped"):
                read.result(5)
        assert not service.thread.is_alive()


def test_tight_stream_budget_degrades_to_one_smaller_buffer(tmp_path):
    budget = ResourceBudget(memory_bytes=4096)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(input_backend="direct")) as (service, client, _):
        file = client()
        with file["moment0"].iter_blocks(frames=32) as blocks:
            with next(blocks) as block:
                assert block.nbytes <= 3072
                assert len(service._streams[blocks.identity]["buffers"]) == 1
        file.client.rpc("dtype", "hd", "moment0")
        assert budget.used["input"] == 0


def test_cancelling_speculation_keeps_an_urgent_request(tmp_path):
    entered, proceed = Event(), Event()
    with inputs(tmp_path, policy=InputPolicy(lookahead=True, input_backend="replay")) as (service, client, values):
        original = service._read_direct
        def held(dataset, destination, selection):
            if selection.start == 0:
                entered.set()
                assert proceed.wait(5)
            original(dataset, destination, selection)
        service._read_direct = held
        service.prefetch(("hd", "moment0"))
        assert entered.wait(5)
        file = client()
        service.prefetch(("hd", "moment0"), first=32)
        file.client.send("read", "hd", "moment0", (64, 68))
        service.cancel_prefetch()
        proceed.set()
        reference, error = file.client.responses.get(timeout=5)
        assert error is None
        mapping = reference.open()
        try:
            np.testing.assert_array_equal(mapping, values[64:68])
        finally:
            mapping._mmap.close()
        assert service._replays[service.aliases["hd", "moment0"]].states[1] == "not_requested"


def test_metadata_only_spectra_inventory_never_creates_bulk_backing(tmp_path):
    from input_output.readers.h5 import RawH5SourceReader
    from input_output.schema.holodoppler import HolodopplerSource
    with h5py.File(tmp_path / "spectral.h5", "w") as source, ProductStorage() as storage:
        source.create_dataset("spectrograms", shape=(128, 16, 32, 32), dtype="f4")
        service = InputService({"hd": source}, storage, mp.get_context("spawn"))
        client = service.client()
        try:
            graph = StageGraph((StageSpec("source", _never_run, inputs=(
                InputUse("hd.spectrograms", InputAccess.METADATA),)),))
            service.configure({"spectra": graph})
            adapter = HolodopplerSource(RawH5SourceReader(h5file=service.files(client)["hd"], label="HD"))
            assert adapter.spectrograms().shape == (128, 16, 32, 32)
            assert service.manifest["hd", "spectrograms"]["access"] == {InputAccess.METADATA}
            assert not service.source_reads and not service._datasets and not service._replays
        finally:
            client.close()
            service.close()


def test_native_handle_and_chunk_cache_survive_between_reads(tmp_path):
    with inputs(tmp_path) as (service, client, _):
        file = client()
        file["moment0"][:2]
        dataset = service._datasets[service.aliases["hd", "moment0"]]
        file["M0"][40:42]
        assert service._datasets[service.aliases["hd", "moment0"]] is dataset
        assert dataset.id.valid
    assert not dataset.id.valid


def test_speculative_job_count_is_bounded(tmp_path):
    entered, proceed = Event(), Event()
    values = np.arange(1024 * 4 * 8, dtype="f4").reshape(1024, 4, 8)
    with inputs(tmp_path, values, policy=InputPolicy(max_prefetch_jobs=3, input_backend="replay")) as (service, client, _):
        original_read, original_enqueue = service._read_direct, service._enqueue
        counts = []
        def held(dataset, destination, selection):
            if selection.start == 0:
                entered.set()
                assert proceed.wait(5)
            original_read(dataset, destination, selection)
        def counted(*args):
            result = original_enqueue(*args)
            counts.append(sum(job[0] > 0 and job[2] == "replay" for job in service._jobs.values()))
            return result
        service._read_direct, service._enqueue = held, counted
        service.prefetch(("hd", "moment0"))
        assert entered.wait(5)
        for region in range(1, 31):
            service.prefetch(("hd", "moment0"), first=region * 32)
        proceed.set()
        client().client.rpc("dtype", "hd", "moment0")
        assert counts and max(counts) <= 3


def test_small_declared_dv_array_uses_one_bounded_read(tmp_path):
    from input_output.schema.doppler_view import DV_ARTERY_MASK_PATH, DV_VEIN_MASK_PATH
    with h5py.File(tmp_path / "dv.h5", "w") as source, ProductStorage() as storage:
        values = np.ones((128, 128), dtype="u1")
        source[DV_ARTERY_MASK_PATH], source[DV_VEIN_MASK_PATH] = values, values
        service = InputService({"dv": source}, storage, mp.get_context("spawn"),
                               policy=InputPolicy(input_backend="replay"))
        client = service.client()
        try:
            service.configure({"topology": StageGraph((StageSpec("inputs", _never_run, inputs=(
                InputUse("dv.vessels", InputAccess.SMALL),)),))})
            np.testing.assert_array_equal(service.files(client)["dv"][DV_ARTERY_MASK_PATH][:], values)
            assert service.storage.metrics.snapshot()["source_reads"] == 1
            assert len(service._replays["dv", DV_ARTERY_MASK_PATH].states) == 1
        finally:
            client.close()
            service.close()


def _prefetch_graph(*uses):
    return StageGraph((StageSpec("read", _never_run, inputs=uses),))


def test_future_plan_progresses_fairly_and_reuses_optional_alias(tmp_path):
    completed = Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay", prefetch_regions=2)) as (service, client, values):
        order = []
        original = service._read_direct
        def measured(dataset, destination, selection):
            order.append((dataset.name.strip("/"), selection.start, service._priority))
            original(dataset, destination, selection)
            if len(order) == 4:
                completed.set()
        service._read_direct = measured
        graphs = {
            "topology": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL, optional=True)),
            "velocity": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL),
                                        InputUse("hd.moment2", InputAccess.SEQUENTIAL)),
        }
        service.configure(graphs, execution_order=tuple(graphs))
        service.activate("topology")
        assert completed.wait(5)
        file = client()
        file.client.rpc("dtype", "hd", "moment0")
        assert order == [("M0", 0, 3), ("moment2", 0, 3), ("M0", 32, 3), ("moment2", 32, 3)]
        assert "velocity" not in service._active
        np.testing.assert_array_equal(file["moment0"][:32], values[:32])
        np.testing.assert_array_equal(file["M0"][:32], values[:32])
        assert sum(name == "M0" and first == 0 for name, first, _ in order) == 1


def test_plan_ignores_optional_disabled_inputs_but_prefetches_forced_replay_stream(tmp_path):
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay")) as (service, client, _):
        graph = _prefetch_graph(
            InputUse("hd.moment0", InputAccess.SEQUENTIAL, optional=True),
            InputUse("hd.moment2", InputAccess.STREAM),
            InputUse("hd.band_low", InputAccess.SEQUENTIAL, when="frequency_bands"))
        service.configure({"future": graph}, {"velocity_workflow_selection": "doppler_moments"},
                          execution_order=("future",))
        client().client.rpc("dtype", "hd", "moment0")
        # Explicit replay now supports early streaming inputs too. Optional
        # moment0 and the disabled band must still generate no work.
        assert set(service._replays) <= {service.aliases["hd", "moment2"]}
        np.testing.assert_array_equal(client()["moment2"][:2], np.arange(96 * 4 * 8, dtype="f4").reshape(96, 4, 8)[:2] + 1)
        assert set(service.source_reads) == {service.aliases["hd", "moment2"]}


def test_failure_eligibility_cancels_future_progress_at_boundary(tmp_path):
    entered, proceed = Event(), Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay")) as (service, client, _):
        original = service._read_direct
        order = []
        def held(dataset, destination, selection):
            order.append(selection.start)
            entered.set()
            assert proceed.wait(5)
            original(dataset, destination, selection)
        service._read_direct = held
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))},
                          execution_order=("future",))
        assert entered.wait(5)
        service.update_prefetch(())
        proceed.set()
        client().client.rpc("dtype", "hd", "moment0")
        assert order == [0]
        replay = service._replays[service.aliases["hd", "moment0"]]
        assert replay.states == ["ready", "not_requested", "not_requested"]


def test_continuous_plan_demand_preempts_future_and_preserves_pair_priority(tmp_path):
    entered, proceed, finished = Event(), Event(), Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay")) as (service, client, values):
        original = service._read_direct
        order = []
        def held(dataset, destination, selection):
            order.append((dataset.name.strip("/"), selection.start, service._priority))
            if len(order) == 1:
                entered.set()
                assert proceed.wait(5)
            original(dataset, destination, selection)
            if len(order) == 4:
                finished.set()
        service._read_direct = held
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL),
                                                     InputUse("hd.moment2", InputAccess.SEQUENTIAL))},
                          execution_order=("future",))
        assert entered.wait(5)
        file = client()
        # An explicit queued speculative region must be promoted, not duplicated.
        service.prefetch(("hd", "moment0"), first=64)
        file.client.send("read", "hd", "moment0", (64, 68))
        proceed.set()
        reference, error = file.client.responses.get(timeout=5)
        assert error is None
        array = reference.open()
        try:
            np.testing.assert_array_equal(array, values[64:68])
        finally:
            array._mmap.close()
        assert finished.wait(5)
        assert order[1:3] == [("M0", 64, 0), ("moment2", 64, 1)]
        assert service.storage.metrics.snapshot()["input_priority_changes"] >= 1
        assert sum(name == "M0" and first == 64 for name, first, _ in order) == 1


def test_plan_byte_limit_and_memory_decline_do_not_spin_or_starve_demand(tmp_path):
    budget = ResourceBudget(memory_bytes=4096)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(input_backend="replay", replay_frames=32)) as (service, client, _):
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))},
                          execution_order=("future",))
        client().client.rpc("dtype", "hd", "moment0")
        assert not service.source_reads
        assert budget.used["input"] == 0
        assert service.storage.metrics.snapshot().get("input_speculation_declined_memory", 0) <= 2
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay", max_prefetch_bytes=4096)) as (service, client, values):
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))},
                          execution_order=("future",))
        file = client()
        file.client.rpc("dtype", "hd", "moment0")
        # This demand is a barrier behind the first speculative region.
        np.testing.assert_array_equal(file["moment0"][:32], values[:32])
        service.cancel_prefetch()
        file.client.rpc("dtype", "hd", "moment0")
        assert service.storage.metrics.snapshot()["input_replay_bytes"] <= 8192


def test_future_switch_retains_active_prefetch(tmp_path):
    ready = Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay", future_prefetch=False)) as (service, client, _):
        original = service._read_direct
        def measured(*args):
            original(*args)
            ready.set()
        service._read_direct = measured
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))},
                          execution_order=("future",))
        client().client.rpc("dtype", "hd", "moment0")
        assert not service.source_reads
        service.activate("future")
        assert ready.wait(5)


def test_scheduler_excludes_failed_dependency_closure_from_prefetch():
    from types import SimpleNamespace
    from pipeline_engine.execution.scheduler import CoordinatorScheduler
    # Exercise the notification boundary without running scientific callbacks.
    pipelines = [SimpleNamespace(name=name) for name in ("root", "child", "grandchild", "independent")]
    parents = {"root": (), "child": ("root",), "grandchild": ("child",), "independent": ()}
    notifications = []
    scheduler = CoordinatorScheduler(pipelines=pipelines,
        dag=SimpleNamespace(dependencies_of=lambda name, **kw: parents[name]), options={}, executor=None,
        statuses={"root": "failed"}, failures={}, skips={}, on_input_eligibility=notifications.append)
    scheduler._refresh_input_eligibility()
    assert notifications == [("independent",)]


def test_planner_byte_bound_is_strict_without_consumers(tmp_path):
    bounded = Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay", max_prefetch_bytes=4096)) as (service, client, _):
        original = service._plan_prefetch
        def planned():
            original()
            if service._idle_reason == "resource" and not service._jobs:
                bounded.set()
        service._plan_prefetch = planned
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL),
                                                     InputUse("hd.moment2", InputAccess.SEQUENTIAL))},
                          execution_order=("future",))
        assert bounded.wait(5)
        client().client.rpc("dtype", "hd", "moment0")
        assert service.storage.metrics.snapshot()["input_replay_bytes"] == 4096
        assert sum(state == "ready" for replay in service._replays.values() for state in replay.states) == 1


def test_speculative_replay_rejects_nonlocal_scratch_but_preserves_service(tmp_path):
    declined = Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay")) as (service, client, _):
        original = service._plan_prefetch
        def planned():
            original()
            if service._prefetch_storage_declined:
                declined.set()
        service._plan_prefetch = planned
        with patch("scratch_location.is_local_path", return_value=False):
            service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))},
                              execution_order=("future",))
            assert declined.wait(5)
            file = client()
            assert file.client.rpc("dtype", "hd", "moment0") == np.dtype("f4")
            with pytest.raises(OSError, match="scratch must be local"):
                file["moment0"][:2]
        assert not service._replays and not service.source_reads


def test_stage_topological_order_guides_first_use_without_callbacks(tmp_path):
    from pipeline_engine.stage_graph import ProductSpec
    read = Event()
    with inputs(tmp_path, policy=InputPolicy(input_backend="replay", prefetch_regions=1)) as (service, client, _):
        order = []
        original = service._read_direct
        def measured(dataset, *args):
            order.append(dataset.name.strip("/"))
            original(dataset, *args)
            if len(order) == 2:
                read.set()
        service._read_direct = measured
        graph = StageGraph((
            StageSpec("second", _never_run, requires=("prerequisite",), inputs=(InputUse("hd.moment2", InputAccess.RANDOM),)),
            StageSpec("first", _never_run, products=(ProductSpec("prerequisite"),), inputs=(InputUse("hd.moment0", InputAccess.RANDOM),)),
        ))
        service.configure({"future": graph}, execution_order=("future",))
        assert read.wait(5)
        assert order == ["M0", "moment2"]


@pytest.mark.parametrize("dtype", ["f4", "f8"])
def test_local_native_read_direct_avoids_slice_copy_and_keeps_metrics(tmp_path, dtype):
    with inputs(tmp_path, location="local", policy=InputPolicy()) as (service, client, values):
        file = client()
        dataset = file["moment0"]
        destination = np.full((5, 4, 8), -1, dtype=dtype)
        with patch.object(dataset, "read_slice", side_effect=AssertionError("intermediate array")):
            dataset.read_direct(destination, np.s_[40:42], np.s_[1:3])
        np.testing.assert_array_equal(destination[1:3], values[40:42])
        assert np.all(destination[[0, 3, 4]] == -1)
        assert destination.flags.writeable
        assert file.client._metrics["input_local_native_direct_reads"] == 1
        assert file.client._metrics["input_local_logical_bytes"] == values[40:42].nbytes
        assert len(file.client._local_io) == 1
        assert not service._replays and not service.source_reads


@pytest.mark.parametrize("selection", [np.s_[40:42], np.s_[2:2], np.s_[10:4:-1], [5, 1, 5]])
def test_local_read_direct_fallbacks_and_strided_destination(tmp_path, selection):
    with inputs(tmp_path, location="local", policy=InputPolicy()) as (service, client, values):
        expected = values[selection]
        storage = np.zeros((*expected.shape[:-1], expected.shape[-1] * 2), dtype="f8")
        destination = storage[..., ::2]
        client()["moment0"].read_direct(destination, selection)
        np.testing.assert_array_equal(destination, expected)


def test_local_astype_and_array_convert_in_hdf5_without_original_array(tmp_path):
    with inputs(tmp_path, location="local", policy=InputPolicy()) as (service, client, values):
        dataset = client()["moment0"]
        with patch.object(dataset, "read_slice", side_effect=AssertionError("original dtype allocation")):
            converted = dataset.astype("f8")
            array = np.asarray(dataset, dtype="f8")
        np.testing.assert_array_equal(converted, values.astype("f8"))
        np.testing.assert_array_equal(array, converted)
        assert converted.flags.writeable
        assert not array.flags.writeable
        assert not service._replays


def test_local_narrowing_casts_preserve_numpy_overflow_and_nonfinite_policy(tmp_path):
    values = np.array([np.nan, np.inf, -np.inf, -2.5, 300., -300.], dtype="f4").reshape(6, 1, 1)
    with inputs(tmp_path, values, location="local", policy=InputPolicy()) as (service, client, _):
        dataset = client()["moment0"]
        with np.errstate(invalid="ignore", over="ignore"):
            expected = values.astype("i1")
            converted = dataset.astype("i1")
            destination = np.empty(values.shape, dtype="i1")
            dataset.read_direct(destination)
        np.testing.assert_array_equal(converted, expected)
        np.testing.assert_array_equal(destination, expected)
        assert "input_local_native_direct_reads" not in dataset.file.client._metrics


def test_local_native_read_direct_checks_lifetime_and_reuses_handle(tmp_path):
    with inputs(tmp_path, location="local", policy=InputPolicy()) as (service, client, values):
        file = client()
        dataset = file["moment0"]
        dest = np.empty_like(values)
        dataset.read_direct(dest)
        handle = file.client._local.datasets[file.filename, "moment0"]
        dataset.read_direct(dest)
        assert file.client._local.datasets[file.filename, "moment0"] is handle
        file.client.close()
        assert not handle.id.valid
        with pytest.raises(OSError, match="closed"):
            dataset.read_direct(dest)


@pytest.mark.parametrize("mode,location,expected", [("auto", "network", "replay"), ("auto", "local", "local_direct"),
                                                    ("direct", "network", "nas_direct"), ("replay", "local", "replay")])
def test_access_policy_single_sequential_consumer_and_overrides(tmp_path, mode, location, expected):
    values = np.zeros((96, 64, 64), dtype="f4")
    with inputs(tmp_path, values, location=location, policy=InputPolicy(input_backend=mode)) as (service, client, _):
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))})
        item = service.manifest[service.aliases["hd", "moment0"]]
        assert item["expected_consumers"] == 1
        assert item["backend"] == expected
        assert not service.source_reads


def test_auto_optional_random_use_changes_representation_not_speculative_eligibility(tmp_path):
    ready = Event()
    with inputs(tmp_path, policy=InputPolicy(prefetch_regions=1)) as (service, client, values):
        original = service._read_direct
        def measured(*args):
            original(*args)
            ready.set()
        service._read_direct = measured
        graphs = {"topology": _prefetch_graph(InputUse("hd.moment0", InputAccess.RANDOM, optional=True),
                                             InputUse("hd.moment2", InputAccess.RANDOM, optional=True)),
                  "velocity": _prefetch_graph(InputUse("hd.moment0", InputAccess.SEQUENTIAL))}
        service.configure(graphs, execution_order=tuple(graphs))
        assert ready.wait(5)
        file = client()
        file.client.rpc("dtype", "hd", "moment0")
        item = service.manifest[service.aliases["hd", "moment0"]]
        assert item["expected_consumers"] == 1 and item["potential_consumers"] == 2
        assert item["backend"] == "replay" and item["expected_reuse"]
        assert set(service._replays) == {service.aliases["hd", "moment0"]}
        for _ in range(2):
            np.testing.assert_array_equal(np.asarray(file["moment0"]), values)
        assert sum(value["bytes"] for value in service.dataset_metrics.values()) == values.nbytes


@pytest.mark.parametrize("mode", ["auto", "direct"])
def test_future_stream_adopts_existing_exclusive_buffers_without_replay(tmp_path, mode):
    ready = Event()
    budget = ResourceBudget(memory_bytes=16 * 2**20)
    values = np.arange(96 * 64 * 64, dtype="f4").reshape(96, 64, 64)
    with inputs(tmp_path, values, budget=budget, policy=InputPolicy(input_backend=mode, stream_frames=8)) as (service, client, _):
        original = service._read_stream
        def measured(*args):
            original(*args)
            if sum(b["state"] == "ready" for stream in service._streams.values() for b in stream["buffers"]) == 2:
                ready.set()
        service._read_stream = measured
        service.configure({"spectra": _prefetch_graph(InputUse("hd.moment0", InputAccess.STREAM))}, execution_order=("spectra",))
        assert ready.wait(5)
        assert "spectra" not in service._active
        assert not service._replays
        initial = next(iter(service._streams))
        memories = [b["memory"].name for b in service._streams[initial]["buffers"]]
        assert service._streams[initial]["owner"] is None
        file = client()
        with file["moment0"].iter_blocks(writable=True) as blocks:
            assert blocks.identity == initial
            for lease in blocks:
                with lease as block:
                    np.testing.assert_array_equal(block, values[lease.start:lease.stop])
                    block[:] = -1
        file.client.rpc("dtype", "hd", "moment0")
        assert not service._streams and budget.used["input"] == 0
        assert service.storage.metrics.snapshot()["input_planned_streams_adopted"] == 1
        assert sum(item["bytes"] for item in service.dataset_metrics.values()) == values.nbytes
        for name in memories:
            with pytest.raises(FileNotFoundError):
                SharedMemory(name=name)


def test_future_stream_cancellation_releases_memory_and_keeps_demand_capacity(tmp_path):
    ready = Event()
    budget = ResourceBudget(memory_bytes=32768)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(stream_frames=16, max_streams=1)) as (service, client, values):
        original = service._read_stream
        def measured(*args):
            original(*args)
            ready.set()
        service._read_stream = measured
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.STREAM))}, execution_order=("future",))
        assert ready.wait(5)
        assert budget.used["input"] > 0
        # An unrelated demanded stream may reclaim an unopened planned pool.
        file = client()
        with file["moment2"].iter_blocks(frames=16) as blocks:
            with next(blocks) as block:
                np.testing.assert_array_equal(block, values[:16] + 1)
                assert sum(s["owner"] is not None for s in service._streams.values()) == 1
        service.update_prefetch(())
        file.client.rpc("dtype", "hd", "moment0")
        assert not service._streams and budget.used["input"] == 0


@pytest.mark.parametrize("chunks", [(8, 32, 32), (24, 8, 8)])
def test_spatial_requests_replay_each_source_region_once_with_chunked_compression(tmp_path, chunks):
    values = np.arange(96 * 32 * 32, dtype="f4").reshape(96, 32, 32)
    with h5py.File(tmp_path / "tiled.h5", "w") as source, ProductStorage(directory=tmp_path / "scratch") as storage:
        source.create_dataset("moment0ff", data=values, chunks=chunks, compression="gzip")
        with patch("pipeline_engine.io.input.service.storage_location", return_value="network"):
            service = InputService({"hd": source}, storage, mp.get_context("spawn"), policy=InputPolicy(replay_frames=32))
        client = service.client()
        try:
            service.configure({"topology": _prefetch_graph(InputUse("hd.moment0ff", InputAccess.RANDOM, optional=True))})
            file = service.files(client)["hd"]
            for y in range(0, 32, 8):
                for x in range(0, 32, 8):
                    np.testing.assert_array_equal(file["moment0ff"][:, y:y+8, x:x+8], values[:, y:y+8, x:x+8])
            assert service.manifest["hd", "moment0ff"]["backend"] == "replay"
            assert storage.metrics.snapshot()["input_logical_bytes"] == values.nbytes
            assert storage.metrics.snapshot()["source_reads"] <= 4
            assert service._replays["hd", "moment0ff"].frames % chunks[0] == 0
        finally:
            client.close()
            service.close()


def test_future_stream_queue_bound_cancellation_and_source_failure(tmp_path):
    read = Event()
    budget = ResourceBudget(memory_bytes=16384)
    with inputs(tmp_path, budget=budget, policy=InputPolicy(max_prefetch_jobs=1)) as (service, client, _):
        def fail(*args):
            read.set()
            raise OSError("planned source failure")
        service._read_direct = fail
        service.configure({"future": _prefetch_graph(InputUse("hd.moment0", InputAccess.STREAM))}, execution_order=("future",))
        assert read.wait(5)
        file = client()
        file.client.rpc("dtype", "hd", "moment0")
        assert len(next(iter(service._streams.values()))["buffers"]) == 1
        with pytest.raises(OSError, match="planned source failure"):
            with file["moment0"].iter_blocks() as blocks:
                next(blocks)
        service.update_prefetch(())
        file.client.rpc("dtype", "hd", "moment0")
        assert not service._streams and budget.used["input"] == 0


def test_small_array_planning_uses_actual_replay_region_size_for_byte_bound(tmp_path):
    from input_output.schema.doppler_view import DV_ARTERY_MASK_PATH, DV_VEIN_MASK_PATH
    with h5py.File(tmp_path / "dv.h5", "w") as source, ProductStorage(directory=tmp_path / "scratch") as storage:
        source[DV_ARTERY_MASK_PATH] = np.ones((128, 128), dtype="u1")
        source[DV_VEIN_MASK_PATH] = np.ones((128, 128), dtype="u1")
        service = InputService({"dv": source}, storage, mp.get_context("spawn"),
                               policy=InputPolicy(input_backend="replay", max_prefetch_bytes=4096))
        client = service.client()
        try:
            service.configure({"future": _prefetch_graph(InputUse("dv.vessels", InputAccess.SMALL))}, execution_order=("future",))
            client.rpc("dtype", "dv", DV_ARTERY_MASK_PATH)
            assert not service._replays and not service.source_reads
        finally:
            client.close()
            service.close()


def test_real_spectral_graph_allows_future_buffers_but_skips_absent_video(tmp_path):
    from pipelines.low_freq_spectra.runner import spectral_graph
    for present in (False, True):
        ready = Event()
        with h5py.File(tmp_path / f"spectra-{present}.h5", "w") as source, ProductStorage(directory=tmp_path / f"scratch-{present}") as storage:
            source["marker"] = 0
            if present:
                source.create_dataset("spectrograms", shape=(32, 4, 64, 64), dtype="f4")
            else:
                source.create_group("spectrograms")
            with patch("pipeline_engine.io.input.service.storage_location", return_value="network"):
                service = InputService({"hd": source}, storage, mp.get_context("spawn"))
            client = service.client()
            original = service._read_stream
            def measured(*args):
                original(*args)
                ready.set()
            service._read_stream = measured
            try:
                service.configure({"low_freq_spectra": spectral_graph()}, execution_order=("low_freq_spectra",))
                client.rpc("dtype", "hd", "marker")
                if present:
                    assert ready.wait(5)
                    assert not service._active and not service._replays
                    assert next(iter(service._streams.values()))["owner"] is None
                else:
                    assert not service._streams and not service.source_reads
            finally:
                client.close()
                service.close()
