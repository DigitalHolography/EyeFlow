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

from pipeline_engine.data.input_service import InputPolicy, InputService
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
            with patch("pipeline_engine.data.input_service.storage_location", return_value=location):
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
        original = service._read_direct
        def measured(dataset, destination, selection):
            original(dataset, destination, selection)
            if selection.start == 8:
                second_ready.set()
        service._read_direct = measured
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
