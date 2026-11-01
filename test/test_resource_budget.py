"""Resource scopes must remain safe across concurrent work and barriers."""

from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import pytest

from pipeline_engine.resources import ResourceBudget, suspend_compute
from runtime_limits import thread_allowance


def test_overlapping_calculation_and_publication_do_not_reconfigure_native_pools(monkeypatch):
    native_calls = []
    monkeypatch.setattr("pipeline_engine.resources.configure_native_threads", native_calls.append)
    budget = ResourceBudget(cpu_threads=4, memory_bytes=4096)
    rendezvous = Barrier(2)

    def calculate():
        with budget.compute(2, 128):
            rendezvous.wait(timeout=3)
            assert thread_allowance() == 2
            assert budget.cpu_used == 4
            rendezvous.wait(timeout=3)

    with budget.activate(), ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(calculate) for _ in range(2)]
        for future in futures:
            future.result(timeout=5)
    assert native_calls == [1]
    assert budget.cpu_used == 0
    assert not any(budget.used.values())


def test_compute_exception_and_nested_artifact_barrier_release_reservations():
    budget = ResourceBudget(cpu_threads=2, memory_bytes=1024)
    with pytest.raises(ValueError, match="calculation failed"):
        with budget.compute(2, 128):
            with suspend_compute():
                assert budget.cpu_used == 0
                assert budget.used["transient"] == 0
                assert budget.try_acquire(2, 128)
                budget.release(2, 128)
            assert thread_allowance() == 2
            assert budget.cpu_used == 2
            assert budget.used["transient"] == 128
            raise ValueError("calculation failed")
    assert budget.cpu_used == 0
    assert not any(budget.used.values())
