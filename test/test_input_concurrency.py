"""Replay authorization and coordinator/data-plane ownership under concurrent I/O."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from threading import Event
from time import perf_counter
from unittest.mock import patch
import multiprocessing as mp

import numpy as np
import pytest

from test_input_service import inputs, _spawn_consumer
from pipeline_engine.io.input.service import InputPolicy
from pipeline_engine.io.input.source_reader import SourceCompletion
from pipeline_engine.resources import ResourceBudget


@contextmanager
def held_source(service, dataset="moment2"):
    entered, proceed, finished = Event(), Event(), Event()
    original = service._read_direct
    def held(native, destination, selection):
        if native.name.strip('/') == dataset:
            entered.set()
            assert proceed.wait(5), "test source operation was not released"
        try:
            original(native, destination, selection)
        finally:
            if native.name.strip('/') == dataset:
                finished.set()
    service._read_direct = held
    try:
        yield entered, proceed, finished
    finally:
        proceed.set()
        service._read_direct = original


def test_authorized_superset_alias_and_full_mapping_need_no_more_rpcs(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        file = client()
        granted = file['moment0'][32:64, 1:3, :]
        with patch.object(file.client, 'rpc', side_effect=AssertionError('cache hit made an RPC')):
            for selection in (np.s_[35:40, :2, 1:4], np.s_[40], np.s_[-40:-32], np.s_[35:35]):
                np.testing.assert_array_equal(file['M0'][selection], values[selection])
        np.testing.assert_array_equal(granted, values[32:64,1:3,:])
        with patch.object(file.client, 'rpc', wraps=file.client.rpc) as rpc:
            np.testing.assert_array_equal(file['moment0'][65:68], values[65:68])
            assert rpc.call_count == 1
        full = np.asarray(file['M0'])
        with patch.object(file.client, 'rpc', side_effect=AssertionError('full coverage made an RPC')):
            for selection in (np.s_[..., 1], np.s_[70:10:-3], np.s_[::2], [7,2,7], values[:,0,0] % 2 == 0):
                np.testing.assert_array_equal(file['moment0'][selection], values[selection])
            destination = np.empty_like(values)
            file['M0'].read_direct(destination)
            np.testing.assert_array_equal(destination, values)
        assert file.client._replay_maps.retained_count == 1
        assert not full.flags.writeable
        with pytest.raises(ValueError):
            full.flags.writeable = True
        state = file.client.__getstate__()
        assert state['_replay_maps'] is None and state['cache'] == {}


def test_cache_eviction_does_not_close_a_live_view(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        file = client()
        file.client._replay_maps.limit = 1
        first = file['moment0'][:2]
        second = file['moment0'][40:42]
        np.testing.assert_array_equal(first, values[:2])
        np.testing.assert_array_equal(second, values[40:42])
        assert file.client._replay_maps.retained_count == 1
        with patch.object(file.client, 'rpc', side_effect=AssertionError('live authorized view was forgotten')):
            np.testing.assert_array_equal(file['M0'][1:2], values[1:2])
        file.client.close()
        assert file.client._replay_maps.live_count == 0
        with pytest.raises(OSError):
            file['moment0'][:2]


def test_ready_response_and_same_client_cache_hit_bypass_unrelated_source_read(tmp_path):
    with inputs(tmp_path, policy=InputPolicy(input_backend='replay')) as (service, client, values):
        file, newcomer = client(), client()
        np.testing.assert_array_equal(file['moment0'][:32], values[:32])
        with held_source(service) as (entered, proceed, _), ThreadPoolExecutor(3) as pool:
            service.prefetch(('hd','moment2'))
            assert entered.wait(5)
            blocked = pool.submit(file['moment2'].read_slice, np.s_[40:42])
            cached = pool.submit(file['M0'].read_slice, np.s_[4:8])
            ready = pool.submit(newcomer['moment0'].read_slice, np.s_[4:8])
            try:
                np.testing.assert_array_equal(cached.result(1), values[4:8])
                np.testing.assert_array_equal(ready.result(1), values[4:8])
                assert not blocked.done()
                assert newcomer.client.rpc('dtype','hd','moment0') == values.dtype
                assert service._inflight is not None
            finally:
                proceed.set()
            np.testing.assert_array_equal(blocked.result(5), values[40:42]+1)


def test_partial_tile_failure_and_adjacent_ready_reference(tmp_path):
    values = np.arange(4*16*16,dtype='f4').reshape(4,16,16)
    with inputs(tmp_path, values, policy=InputPolicy(input_backend='replay',lookahead=False,max_io_bytes=128,replay_frames=1)) as (service, client, _):
        file, adjacent = client(), client()
        np.testing.assert_array_equal(file['moment0'][:1],values[:1])
        entered, proceed = Event(), Event()
        original = service._read_direct
        def fail(native, destination, selection):
            if selection[0].start == 1 and selection[1].start == 2:
                destination[:] = -999
                entered.set()
                assert proceed.wait(5)
                raise OSError('partial tile failed')
            original(native,destination,selection)
        service._read_direct = fail
        with ThreadPoolExecutor(2) as pool:
            missing = pool.submit(file['moment0'].read_slice,np.s_[1:2])
            assert entered.wait(5)
            try:
                authorized = pool.submit(adjacent['M0'].read_slice,np.s_[:1])
                np.testing.assert_array_equal(authorized.result(1),values[:1])
                assert not missing.done()
            finally:
                proceed.set()
            with pytest.raises(OSError,match='partial tile failed'):
                missing.result(5)
        replay=service._replays[service.aliases['hd','moment0']]
        assert replay.states[:2] == ['ready','failed']
        np.testing.assert_array_equal(adjacent['moment0'][:1],values[:1])


def test_flush_failure_is_not_published_and_reader_survives_baseexception(tmp_path):
    with inputs(tmp_path) as (service, client, values):
        file, other = client(), client()
        np.testing.assert_array_equal(file['moment0'][:2],values[:2])
        replay=service._replays[service.aliases['hd','moment0']]
        with patch.object(replay.array,'flush',side_effect=OSError('flush failed')):
            with pytest.raises(OSError,match='flush failed'):
                file['moment0'][40:42]
        np.testing.assert_array_equal(other['M0'][:2],values[:2])
        original = service._read_direct
        def fail(native,*args):
            if native.name == '/moment2': raise SystemExit('reader task exited')
            original(native,*args)
        service._read_direct=fail
        with pytest.raises(OSError,match='reader task exited'):
            file['moment2'][:2]
        assert service.thread.is_alive() and service.reader.thread.is_alive()
        np.testing.assert_array_equal(other['moment0'][:2],values[:2])


def test_stale_completion_never_publishes_the_active_destination(tmp_path):
    with inputs(tmp_path,policy=InputPolicy(input_backend='replay')) as (service,client,values):
        file=client()
        with held_source(service) as (entered,proceed,_),ThreadPoolExecutor(1) as pool:
            request=pool.submit(file['moment2'].read_slice,np.s_[:2])
            assert entered.wait(5)
            token=service._inflight['token']
            service.reader.completed.put_nowait(SourceCompletion(token+1000,None,None,(),perf_counter(),perf_counter(),0.))
            service._wake_completion()
            barrier=client()
            barrier.client.rpc('dtype','hd','moment0')
            try:
                assert not request.done() and service._inflight['token']==token
                assert service.storage.metrics.snapshot()['input_stale_completions']==1
            finally:
                proceed.set()
            np.testing.assert_array_equal(request.result(5),values[:2]+1)


def test_shutdown_cancels_rpc_before_bounded_read_finishes(tmp_path):
    with inputs(tmp_path) as (service,client,_):
        file=client()
        with held_source(service) as (entered,proceed,_),ThreadPoolExecutor(2) as pool:
            request=pool.submit(file['moment2'].read_slice,np.s_[:2])
            assert entered.wait(5)
            closing=pool.submit(service.close)
            try:
                with pytest.raises(OSError,match='stopped|cancelled'):
                    request.result(1)
                file.client.close()
                assert not closing.done()
            finally:
                proceed.set()
            closing.result(5)
        assert not service.thread.is_alive() and not service.reader.thread.is_alive()


def test_reclaim_during_fill_defers_free_until_completion(tmp_path):
    budget=ResourceBudget(memory_bytes=64*2**20)
    with inputs(tmp_path,budget=budget,policy=InputPolicy(input_backend='direct')) as (service,client,_):
        file,barrier=client(),client()
        context=mp.get_context('spawn')
        results=context.Queue()
        process=context.Process(target=_spawn_consumer,args=(file,results,True))
        retired=Event()
        close=service._close_stream
        def observed(identity):
            close(identity)
            if identity not in service._streams: retired.set()
        service._close_stream=observed
        with held_source(service,'M0') as (entered,proceed,_):
            process.start()
            try:
                assert entered.wait(15)
                process.terminate(); process.join(5)
                service.release_client(file.client.identity)
                barrier.client.rpc('dtype','hd','moment0')
                assert next(iter(service._streams.values()))['closing']
                assert budget.used['input'] > 0
                proceed.set()
                assert retired.wait(5)
                assert not service._streams and budget.used['input']==0
            finally:
                proceed.set()
                if process.is_alive(): process.terminate()
                process.join(); process.close(); results.close()


def test_parallel_clients_deduplicate_regions_and_keep_values_exact(tmp_path):
    with inputs(tmp_path) as (service,client,values):
        files=[client() for _ in range(8)]
        def read(i):
            first=(i*7)%90
            actual=files[i%len(files)]['M0'][first:first+4]
            np.testing.assert_array_equal(actual,values[first:first+4])
        with ThreadPoolExecutor(12) as pool:
            list(pool.map(read,range(120)))
        assert service.storage.metrics.snapshot()['source_reads']==3
        assert all(f.client._replay_maps.retained_count<=32 for f in files)


def test_reaped_consumer_does_not_continue_an_unfinished_multitile_region(tmp_path):
    values=np.zeros((4,16,16),dtype='f4')
    with inputs(tmp_path,values,policy=InputPolicy(input_backend='replay',lookahead=False,max_io_bytes=128,replay_frames=1)) as (service,client,_):
        file,barrier=client(),client()
        completed=Event()
        finish=service._complete_replay
        def observed(*args):
            finish(*args)
            completed.set()
        service._complete_replay=observed
        with held_source(service,'M0') as (entered,proceed,_):
            file.client.send('read','hd','moment0',(0,1))
            assert entered.wait(5)
            service.release_client(file.client.identity)
            barrier.client.rpc('dtype','hd','moment0')
            proceed.set()
            assert completed.wait(5)
        barrier.client.rpc('dtype','hd','moment0')
        assert service.storage.metrics.snapshot()['source_reads']==1
        assert not service._jobs and service._inflight is None
        assert service._replays[service.aliases['hd','moment0']].states[0]=='not_requested'


@pytest.mark.parametrize('backend', ['auto', 'direct', 'replay'])
def test_repeated_startup_and_exception_cleanup_retire_all_input_ownership(tmp_path, backend):
    budget = ResourceBudget(memory_bytes=64*2**20)
    for cycle in range(4):
        directory = tmp_path / str(cycle)
        directory.mkdir()
        with inputs(directory, budget=budget, policy=InputPolicy(input_backend=backend)) as (service, client, values):
            file = client()
            np.testing.assert_array_equal(file['M0'][:2], values[:2])
            # Deliberately leave the stream after an exception with a live lease.
            with pytest.raises(RuntimeError, match='consumer failed'):
                with file['moment2'].iter_blocks(frames=4, writable=True) as blocks:
                    with next(blocks) as block:
                        np.testing.assert_array_equal(block, values[:4]+1)
                        raise RuntimeError('consumer failed')
            replay_paths = [item.reference.path for item in service._replays.values()]
        assert not service.thread.is_alive() and not service.reader.thread.is_alive()
        assert not service._streams and not service._transfers and not service._datasets
        assert file.client._replay_maps.live_count == 0
        assert budget.used['input'] == 0
        assert all(not Path(path).exists() for path in replay_paths)
