"""Spawned output ownership, accepted-only writes, backpressure and failures."""
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from functools import partial
from hashlib import sha256
import gc
import os
import multiprocessing as mp
from pathlib import Path
import pickle

import h5py
import numpy as np
import pytest

from input_output.writers.h5_publication import OutputPublisher
from pipeline_engine.data.storage import ProductStorage, ArrayReference, H5ArrayReference
from pipeline_engine.io.manager import IOManager
from pipeline_engine.io.output.protocol import OutputPolicy, OutputServiceError, WriteCommand, WriterConnection
from pipeline_engine.io.output.service import OutputService
from pipeline_engine.io.output.writer import writer_main
from pipeline_engine.output.output_model import OutputModel, Record
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.io.input.service import InputPolicy


@contextmanager
def output(tmp_path, *, policy=None, budget=None, target=None):
    with ProductStorage(directory=tmp_path/'scratch', budget=budget) as storage:
        model = OutputModel(storage, tmp_path/'final.h5')
        service = OutputService(model, policy=policy, budget=budget, _writer_target=target)
        try:
            yield model, service
        finally:
            service.close()


def gated_writer(gate, *args):
    if not gate.wait(15):
        raise RuntimeError('writer gate timed out')
    writer_main(*args)


def locked_hdf5_writer(entered, release, *args):
    with h5py._objects.phil:
        entered.set()
        if not release.wait(15): raise RuntimeError('writer lock gate timed out')
        writer_main(*args)


def test_explicit_network_scratch_is_rejected_before_ipc_allocation(tmp_path, monkeypatch):
    context = mp.get_context('spawn')
    with ProductStorage(directory=tmp_path/'scratch') as storage:
        model = OutputModel(storage, tmp_path/'final.h5')
        monkeypatch.setattr('pipeline_engine.io.output.service.is_local_path', lambda path: False)
        def unexpected_ipc(*args, **kwargs):
            raise AssertionError('Network output scratch must be rejected before allocating IPC')
        monkeypatch.setattr(context, 'Pipe', unexpected_ipc)
        with pytest.raises(OSError, match='output scratch must be local'):
            OutputService(model, context=context)
        assert not list(storage.directory.iterdir())
        assert model._commit_listener is None


def test_duplicate_service_does_not_spawn_or_detach_the_existing_writer(tmp_path, monkeypatch):
    with output(tmp_path) as (model, service):
        def unexpected_spawn(*args, **kwargs):
            raise AssertionError('A duplicate writer must not be spawned')
        monkeypatch.setattr(mp.process.BaseProcess, 'start', unexpected_spawn)
        with pytest.raises(RuntimeError, match='listener is already bound'):
            OutputService(model)
        assert model._commit_listener == service.submit
        with model.checkpoint():
            model.create_dataset('accepted', data=7)
        service.drain()
        assert service.written_paths == {'accepted'}
        with h5py.File(service.finish()) as physical:
            assert physical['accepted'][()] == 7


def test_manifest_interrupt_stops_the_batch_at_a_safe_health_boundary(tmp_path, monkeypatch):
    with output(tmp_path) as (model, service):
        interrupt = KeyboardInterrupt('stop output')
        def interrupted(*args, **kwargs):
            raise interrupt
        monkeypatch.setattr(service, '_enqueue', interrupted)
        with model.checkpoint():
            model.create_dataset('accepted', data=7)
        assert model.records['accepted'].committed
        with pytest.raises(KeyboardInterrupt) as caught:
            service.check_health()
        assert caught.value is interrupt
        assert not Path(model.filename).exists()


def test_committed_data_is_in_incomplete_file_before_analysis_finishes(tmp_path):
    with output(tmp_path) as (model, service):
        values = np.array([1., -0., np.nan], dtype='f4')
        with model.checkpoint(owner='early'):
            model.create_dataset('Early/value', data=values).attrs['unit']='mm/s'
        service.drain()
        assert service.completed_sequence == 1 and not Path(model.filename).exists()
        # Only the writer can open the incomplete HDF5 on Windows. Its ACK is
        # sent after creation and flush, before the later scientific work below.
        assert service.written_paths == {'Early/value'} and service.local_path.stat().st_size > 0
        with pytest.raises(RuntimeError), model.checkpoint(owner='later'):
            model.create_dataset('Rejected/value', data=9)
            raise RuntimeError('calculation failed')
        # An uncommitted bare edit must also stay out of the writer inventory.
        model.create_dataset('Uncommitted', data=17)
        local=service.finish()
        with h5py.File(local) as physical:
            assert 'Early/value' in physical and 'Rejected' not in physical
            assert 'Uncommitted' not in physical
            assert physical['Early/value'][:].tobytes()==values.tobytes()


def test_final_metadata_replaces_typed_values_and_removes_stale_keys(tmp_path):
    with output(tmp_path) as (model, service):
        with model.checkpoint():
            group=model.require_group('Shared')
            group.attrs['old']='remove'
            group.attrs.create('count',1,dtype=np.uint8)
            model.attrs['pipeline_status']='running'
            model.create_dataset('Shared/value',data=3)
        service.drain()
        del group.attrs['old']
        group.attrs['count']=300
        group.attrs.create('typed',7,dtype=np.uint16)
        model.attrs['pipeline_status']='failed'
        with h5py.File(service.finish()) as physical:
            assert 'old' not in physical['Shared'].attrs
            assert physical['Shared'].attrs['count']==300
            assert physical['Shared'].attrs.get_id('typed').dtype==np.dtype('u2')
            assert physical.attrs['pipeline_status']=='failed'


def test_full_backlog_does_not_hold_the_coordinator_or_single_cpu(tmp_path):
    budget=ResourceBudget(cpu_threads=1,memory_bytes=512*2**20)
    policy=OutputPolicy(max_pending_jobs=1,max_pending_bytes=256*1024)
    with budget.activate(), output(tmp_path,policy=policy,budget=budget) as (model,service):
        with budget.compute(1,0), model.checkpoint():
            for i in range(7):
                model.create_dataset(f'value{i}',data=np.full((256,256),i,dtype='f4'))
        assert service.pending_count <= 1
        assert model.storage.metrics.snapshot()['output_submission_deferred'] > 0
        with h5py.File(service.finish()) as physical:
            assert len(physical)==7
            for i in range(7): assert physical[f'value{i}'][0,0]==i
        assert budget.cpu_used==0 and sum(budget.used.values())==0


def test_ipc_contains_only_manifest_identity_and_large_arrays_keep_their_reference(tmp_path):
    gate=mp.get_context('spawn').Event()
    with output(tmp_path,target=partial(gated_writer,gate)) as (model,service):
        try:
            with model.checkpoint():
                model.create_dataset('large',data=np.arange(1024*1024,dtype='f4'))
            command=service._pending[0][0]
            assert isinstance(command,WriteCommand) and len(pickle.dumps(command))<1024
            with Path(command.manifest).open('rb') as stream:
                records=pickle.load(stream)
            assert isinstance(records[0][1].value,ArrayReference)
            assert Path(records[0][1].value.path).exists()
        finally:
            gate.set()
        service.finish()
        assert not model.storage.metrics.snapshot().get('array_copied_bytes',0) > 4*2**20


def test_failed_write_preserves_logical_commits_and_prevents_publication(tmp_path):
    with output(tmp_path) as (model,service):
        with model.checkpoint():
            model.create_dataset('accepted',data=[1,2])
        service.drain()
        with model.checkpoint():
            model.create_dataset('broken',data=np.ones((32,32)),compression='invalid-codec')
        assert model.records['accepted'].committed and model.records['broken'].committed
        with pytest.raises(OutputServiceError,match='compression|filter'):
            service.finish()
        assert not Path(model.filename).exists()


def test_writer_death_fails_pending_drain_and_close_releases_resources(tmp_path):
    gate=mp.get_context('spawn').Event()
    with output(tmp_path,target=partial(gated_writer,gate)) as (model,service):
        with model.checkpoint(): model.create_dataset('value',data=1)
        service.process.terminate()
        service.process.join(5)
        with pytest.raises(OutputServiceError,match='startup|exited'):
            service.finish()
    assert not service._thread.is_alive() and not service.local_path.exists()
    assert not list((tmp_path/'scratch').glob('*.pickle'))


def test_close_during_startup_is_bounded_and_idempotent(tmp_path):
    gate=mp.get_context('spawn').Event()
    with output(tmp_path,target=partial(gated_writer,gate)) as (_,service):
        service.close(); service.close()
        assert not service._thread.is_alive()
        assert service.commands is None and service.acknowledgements is None


@pytest.mark.parametrize('failure',['copy','replace'])
def test_publication_failure_keeps_existing_file_and_removes_incoming(tmp_path,monkeypatch,failure):
    local=tmp_path/'local.h5'
    with h5py.File(local,'w') as physical: physical['value']=7
    destination=tmp_path/'published.h5'
    destination.write_bytes(b'previous result')
    def fail(*args): raise OSError('destination failed')
    monkeypatch.setattr(f'input_output.writers.h5_publication.{"shutil.copyfile" if failure=="copy" else "os.replace"}',fail)
    with pytest.raises(OSError,match='destination failed'):
        OutputPublisher(destination).publish(local)
    assert destination.read_bytes()==b'previous result'
    assert not list(tmp_path.glob('.ef-*')) and local.exists()


def test_two_acquisitions_retain_independent_writers_after_input_closes(tmp_path):
    a=IOManager(OutputModel(ProductStorage(),tmp_path/'a.h5'))
    b=IOManager(OutputModel(ProductStorage(),tmp_path/'b.h5'))
    try:
        assert a.output.process.pid != b.output.process.pid
        with a.model.checkpoint(): a.model.create_dataset('value',data=1)
        a.close_input()
        assert a.output.process.is_alive() and not (tmp_path/'a.h5').exists()
        with b.model.checkpoint(): b.model.create_dataset('value',data=2)
        b.finish(); a.finish()
        assert a.finish()==tmp_path/'a.h5'
        with h5py.File(tmp_path/'a.h5') as physical: assert physical['value'][()]==1
        with h5py.File(tmp_path/'b.h5') as physical: assert physical['value'][()]==2
    finally:
        a.close();b.close()


@pytest.mark.parametrize('conflict',['sequence','dataset','inventory','gap'])
def test_duplicate_sequences_are_idempotent_and_conflicting_delivery_is_terminal(tmp_path,conflict):
    context=mp.get_context('spawn')
    receiver,commands=context.Pipe(duplex=False)
    acks,sender=context.Pipe(duplex=False)
    local=tmp_path/'incomplete.h5'
    process=context.Process(target=writer_main,args=(WriterConnection(receiver),WriterConnection(sender),str(local),2**20))
    process.start()
    receiver.close();sender.close()
    def receive():
        assert acks.poll(10)
        return acks.recv()
    try:
        assert receive().sequence==0
        manifest=tmp_path/'record.pickle'
        manifest.write_bytes(pickle.dumps([('value',Record('dataset',np.array([7]),committed=True))]))
        command=WriteCommand(1,str(manifest),sha256(manifest.read_bytes()).hexdigest())
        for _ in range(2):
            commands.send(command)
            ack=receive()
            assert ack.sequence==1 and not ack.error
        if conflict=='sequence':
            invalid=WriteCommand(1,str(manifest),'conflicting-digest')
        elif conflict=='gap':
            invalid=WriteCommand(3,str(manifest),command.digest)
        else:
            another=tmp_path/'invalid.pickle'
            contents=({'':Record(committed=True)},('value','missing')) if conflict=='inventory' else [('value',Record('dataset',np.array([8]),committed=True))]
            another.write_bytes(pickle.dumps(contents))
            invalid=WriteCommand(2,str(another),sha256(another.read_bytes()).hexdigest(),conflict=='inventory')
        commands.send(invalid)
        message=receive().error
        assert message and any(word in message for word in ('Conflicting','inventory','Out-of-order'))
        process.join(5)
        with h5py.File(local) as physical: assert physical['value'][0]==7
    finally:
        if process.is_alive(): process.terminate();process.join(5)
        process.close()
        for channel in (commands,acks): channel.close()


def test_missing_backing_and_final_metadata_errors_never_publish(tmp_path):
    gate=mp.get_context('spawn').Event()
    with output(tmp_path,target=partial(gated_writer,gate)) as (model,service):
        with model.checkpoint(): model.create_dataset('large',data=np.ones(100000,dtype='f4'))
        Path(model.records['large'].value.path).unlink()
        gate.set()
        with pytest.raises(OutputServiceError,match='FileNotFound|No such file'):
            service.finish()
        assert model.records['large'].committed and not Path(model.filename).exists()
    with output(tmp_path/'metadata') as (model,service):
        with model.checkpoint(): model.create_dataset('valid',data=7)
        service.drain()
        model.attrs['unsupported']={'not':'an HDF5 attribute'}
        with pytest.raises(OutputServiceError): service.finish()
        assert model.records['valid'].committed and not Path(model.filename).exists()


def test_input_progresses_while_output_process_holds_its_hdf5_lock(tmp_path,monkeypatch):
    context=mp.get_context('spawn')
    entered,release=context.Event(),context.Event()
    monkeypatch.setattr('pipeline_engine.io.output.service.writer_main',partial(locked_hdf5_writer,entered,release))
    source=tmp_path/'source.h5'
    values=np.arange(64*4*8,dtype='f4').reshape(64,4,8)
    with h5py.File(source,'w') as file: file['video']=values
    io=IOManager(OutputModel(ProductStorage(),tmp_path/'final.h5'))
    try:
        assert entered.wait(10)
        io.open_sources(source,None)
        io.start_input(InputPolicy(input_backend='replay',lookahead=False))
        with io.model.checkpoint(): io.model.create_dataset('result',data=1)
        np.testing.assert_array_equal(io.files['hd']['video'][:32],values[:32])
        other=io.input.client()
        try:
            np.testing.assert_array_equal(io.input.files(other)['hd']['video'][:32],values[:32])
            assert other.rpc('dtype','hd','video')==values.dtype
        finally: other.close()
        assert io.output.pending_count and not release.is_set()
        io.close_input()
        release.set()
        io.finish()
    finally:
        release.set();io.close()


def test_nonexposed_stage_serializes_before_later_failure_in_same_invocation(tmp_path):
    from input_output.output_manager import OutputManager
    from input_output.runs.layout import HoloRunLayout
    from pipeline_engine.context import PipelineContext, managed_invocation
    from pipeline_engine.execution.invocation import _pipeline_checkpoint
    from pipeline_engine.stage_graph import StageGraph,StageSpec,ProductSpec
    from pipeline_engine.stage_graph.executor import StageGraphExecutor
    with output(tmp_path) as (model,service):
        manager=OutputManager(HoloRunLayout.from_holo(tmp_path/'scan.holo'))
        context=PipelineContext(work_h5=model,holodoppler_h5=None,doppler_vision_h5=None,
            output_manager=manager,pipeline_name='science',variables={})
        def produce(ctx):
            ctx.dag.handoff('value',np.array([7.,-0.,np.nan]))
            service.drain()
            assert service.written_paths=={'Scientific/value'}
            assert not Path(model.filename).exists()
            raise ValueError('later failed')
        graph=StageGraph((StageSpec('produce',produce,products=(ProductSpec('value','Scientific/value',{'unit':'Hz'}),)),))
        with pytest.raises(RuntimeError,match='later failed'),managed_invocation(context):
            with _pipeline_checkpoint(model,context.state.raw,manager,owner='science') as checkpoint:
                StageGraphExecutor(graph,context,checkpoint,_pipeline_checkpoint).run()
        with h5py.File(service.finish()) as physical:
            np.testing.assert_array_equal(physical['Scientific/value'][:],[7.,-0.,np.nan])


def test_native_hdf5_backing_preserves_payload_attributes_without_dense_materialization(tmp_path):
    with output(tmp_path) as (model,service):
        backing=model.storage.path('.h5')
        shape=(1024,1024,256)
        with h5py.File(backing,'w') as file:
            data=file.create_dataset('data',shape=shape,dtype='f4',chunks=(2,2,2),compression='gzip',fillvalue=np.nan)
            data[:2,:2,:2]=3
            data.attrs['native_provenance']='retained'
        with model.checkpoint():
            model.create_dataset('value',data=H5ArrayReference(str(backing),shape,np.dtype('f4'))).attrs['unit']='Hz'
        with h5py.File(service.finish()) as file:
            assert file['value'].id.get_num_chunks()==1
            assert file['value'].attrs['native_provenance']=='retained'
            assert file['value'].attrs['unit']=='Hz' and file['value'][0,0,0]==3


def test_cancellation_during_serialization_wakes_finalizer_and_reaps_writer(tmp_path,monkeypatch):
    gate=mp.get_context('spawn').Event()
    monkeypatch.setattr('pipeline_engine.io.output.service.writer_main',partial(gated_writer,gate))
    io=IOManager(OutputModel(ProductStorage(),tmp_path/'final.h5'))
    with io.model.checkpoint(): io.model.create_dataset('value',data=7)
    try:
        with ThreadPoolExecutor(1) as pool:
            future=pool.submit(io.finish)
            io.close()
            with pytest.raises((OSError,RuntimeError),match='cancelled|closed|aborted'):
                future.result(5)
        assert io.closed and not io.output._thread.is_alive()
        assert not (tmp_path/'final.h5').exists()
    finally:
        gate.set();io.close()


def test_cancelled_publication_never_replaces_existing_destination(tmp_path,monkeypatch):
    import shutil
    local,final=tmp_path/'local.h5',tmp_path/'final.h5'
    with h5py.File(local,'w') as file: file['value']=1
    final.write_bytes(b'previous result')
    cancelled=[]
    copy=shutil.copyfile
    def copied(*args):
        copy(*args)
        cancelled.append(True)
    monkeypatch.setattr('input_output.writers.h5_publication.shutil.copyfile',copied)
    with pytest.raises(OSError,match='cancelled'):
        OutputPublisher(final,cancelled=lambda:bool(cancelled)).publish(local)
    assert final.read_bytes()==b'previous result' and not list(tmp_path.glob('.ef-*'))


def test_resource_admission_and_ipc_failure_are_terminal_without_losing_commits(tmp_path,monkeypatch):
    class RefusingBudget:
        memory_limit=128*2**20
        def try_acquire(self,*args): return False
        def release(self,*args): raise AssertionError('unacquired reservation was released')
    with output(tmp_path/'admission') as (model,service):
        assert service._ready.wait(10)
        service.budget=RefusingBudget()
        service.policy=replace(service.policy,operation_timeout=0.1)
        with model.checkpoint(): model.create_dataset('accepted',data=7)
        with pytest.raises(OutputServiceError,match='resource admission timed out'):
            service.finish()
        assert model.records['accepted'].committed
    with output(tmp_path/'ipc') as (model,service):
        assert service._ready.wait(10)
        def fail(*args,**kwargs): raise OSError('output IPC failed')
        monkeypatch.setattr(service.commands,'send',fail)
        with model.checkpoint(): model.create_dataset('accepted',data=8)
        with pytest.raises(OutputServiceError,match='IPC failed'): service.finish()
        assert model.records['accepted'].committed and not Path(model.filename).exists()


@pytest.mark.skipif(os.name!='nt',reason='Windows spawn handle ownership')
def test_aborted_startup_releases_transferred_handles_without_cyclic_gc(tmp_path):
    import psutil
    process=psutil.Process()
    retained=[]
    warm=IOManager(OutputModel(ProductStorage(),tmp_path/'warm.h5'))
    warm.finish()
    retained.append(warm)
    baseline=process.num_handles()
    enabled=gc.isenabled()
    gc.disable()
    try:
        for i in range(6):
            manager=IOManager(OutputModel(ProductStorage(),tmp_path/f'aborted{i}.h5'))
            with manager.model.checkpoint(): manager.model.create_dataset('value',data=i)
            manager.close()
            retained.append(manager)
            assert not manager.output._thread.is_alive()
        assert process.num_handles() <= baseline+2
    finally:
        if enabled: gc.enable()
