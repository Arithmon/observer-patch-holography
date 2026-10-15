"""Compatibility with the existing independently interpreted source history."""
from fractions import Fraction as F

import numpy as np
import pytest

from ..common_source import protocol
from ..common_source.open_environment import Environment, PhaseNoise
from ..common_source import open_protocol
from . import common, runtime


@pytest.mark.parametrize('dimension,current', [(2, -1), (2, 0), (2, 1), (4, 0)])
def test_complete_canonical_continuation_keeps_every_map_record_and_resource(dimension, current):
    program = protocol.continuation_program(dimension, 2, native_current=current)
    env = Environment((PhaseNoise(0, 0, F(1, 100), 4), PhaseNoise(3, 1, F(1, 200), 4)))
    report = common.compatibility_report(program, env)
    assert report['complete_cq_trace_distance'] < 2e-13
    assert report['visible_records'] == 16
    assert report['source_resources_agree'] is True


def test_common_source_retains_entangled_reference_and_suppresses_default_preparation_cost():
    p = protocol.Program(dimension=2, steps=2, native_current=1, operations=(
        protocol.Operation('gauge_read'), protocol.Operation('feedback', record=0),
        protocol.Operation('smooth_read')))
    env = Environment((PhaseNoise(0, 1, F(1, 20), 4), PhaseNoise(1, 0, F(1, 10), 4)))
    rng = np.random.default_rng(189)
    vector = rng.normal(size=p.work_dimension*3)+1j*rng.normal(size=p.work_dimension*3)
    vector /= np.linalg.norm(vector)
    initial = np.outer(vector, vector.conj())
    report = common.compatibility_report(p, env, initial=initial)
    assert report['reference_dimension'] == 3
    assert report['complete_cq_trace_distance'] < 2e-13
    output = common.execute(p, env, initial=initial)
    assert all(dict(b.resources)['work_preparations'] == 0 for b in output.branches)


def test_common_checkpoint_replays_then_continues_the_same_recorded_work_state():
    p = protocol.continuation_program(2, 2)
    env = Environment((PhaseNoise(2, 0, F(1, 100), 4),))
    compiled = common.compile_program(p, env)
    prefix = runtime.execute(compiled, stop_after=5)
    restored = runtime.checkpoint_from_json(compiled, runtime.checkpoint_to_json(prefix))
    resumed = runtime.execute(compiled, checkpoint=restored)
    full = runtime.execute(compiled)
    assert runtime.complete_distance(resumed, full) == 0
    assert any(branch.record for branch in prefix.branches)


def test_common_prediction_uses_shared_runtime_and_old_interpreter_only_for_comparison(monkeypatch):
    p = protocol.Program(operations=(protocol.Operation('gauge_read'),))
    def forbidden(*args, **kwargs):
        raise AssertionError('old mixed trajectory used as prediction')
    monkeypatch.setattr(open_protocol, 'execute', forbidden)
    result = common.execute(p)
    assert result.records == ((0,), (1,))
