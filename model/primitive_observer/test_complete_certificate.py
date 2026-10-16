"""Challenge the complete certificate through the actual shared executor."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest

from ..common_source import protocol, quantum_native, smooth_history
from ..common_source.open_environment import Environment
from ..common_source.output_certificate import compare_receipt
from . import common, complete_certificate as certificate, native, runtime


@pytest.fixture(scope='module')
def history():
    return native.execute(Environment(()), 'C', 2, 2)


def test_original_runtime_native_arrays_have_small_exact_reference_uniform_error(history):
    bound = certificate.native_error(history)
    assert type(bound['distance_upper']) is F
    assert 0 < bound['distance_upper'] < F(1, 10**8)
    assert bound['complete_raw_outcomes'] == 770
    assert 0 <= bound['abort_distance_upper'] <= bound['distance_upper']


def test_wrong_but_complete_logical_phase_is_not_an_accepted_native_certificate(monkeypatch):
    original = quantum_native._native_tree
    z = np.diag([1, -1])
    def wrong(law, mode):
        # A logical input Z preserves every branch probability/completeness,
        # but changes an entangled input. It must fail the endpoint guarantee.
        return tuple((c, raw, current, abort, first@z,
                      {bit: tuple((coin, record, stopped, matrix@z)
                                  for coin, record, stopped, matrix in ends)
                       for bit, ends in endings.items()})
                     for c, raw, current, abort, first, endings in original(law, mode))
    monkeypatch.setattr(quantum_native, '_native_tree', wrong)
    altered = native.execute(Environment(()), 'C', 2, 2)
    bound = certificate.native_error(altered)
    assert bound['distance_upper'] > F(1, 2)
    with pytest.raises(ArithmeticError, match='native'):
        certificate.report._gate(bound['distance_upper'], F(1, 10**8), 'native')


def test_all_raw_outcomes_and_history_integrity_are_required(history):
    with pytest.raises(ValueError, match='complete native'):
        certificate.native_error(replace(history, branches=history.branches[:-1]))
    changed = deepcopy(history)
    object.__setattr__(changed.branches[0], 'work_density', np.zeros_like(changed.branches[0].work_density))
    with pytest.raises(ValueError, match='detached'):
        certificate.native_error(changed)
    with pytest.raises(ValueError):
        certificate.native_error(replace(history, native_mode='full_carrier'))


def test_work_comparison_executes_shared_runtime_and_detects_wrong_cptp_feedback(monkeypatch):
    program = protocol.Program(native_current=1, operations=(
        protocol.Operation('gauge_read'), protocol.Operation('feedback', record=0)))
    finite = smooth_history.enclose(program, dimension=2)
    infinite = smooth_history.enclose(program)
    correct = certificate.work_error(common.execute(program), program, finite, infinite)
    assert correct['finite_upper'] < F(1, 10**10)
    original = protocol._operation_maps
    def wrong(p, operation, bit, backend):
        maps, costs = original(p, operation, bit, backend)
        return ((np.eye(p.work_dimension),), costs) if operation.kind == 'feedback' else (maps, costs)
    monkeypatch.setattr(protocol, '_operation_maps', wrong)
    changed = certificate.work_error(common.execute(program), program, finite, infinite)
    assert changed['mass_defect_upper'] < F(1, 10**10)
    assert changed['finite_upper'] > F(1, 2)


def test_work_certificate_rejects_partial_supplied_or_mismatched_histories():
    program = protocol.Program(operations=(protocol.Operation('gauge_read'),))
    finite = smooth_history.enclose(program, dimension=2)
    infinite = smooth_history.enclose(program)
    compiled = common.compile_program(program)
    for worker in (runtime.execute(compiled, stop_after=0),
                   runtime.execute(compiled, initial=compiled.initial_density)):
        with pytest.raises(ValueError, match='complete declared'):
            certificate.work_error(worker, program, finite, infinite)
    with pytest.raises(ValueError):
        certificate.work_error(runtime.execute(compiled), replace(program, native_current=1), finite, infinite)
    for first, second in ((finite, finite), (infinite, finite), (infinite, infinite)):
        with pytest.raises(ValueError, match='finite and infinite'):
            certificate.work_error(runtime.execute(compiled), program, first, second)


def test_weighted_composition_counts_no_record_multiplier_and_does_not_cap():
    assert certificate.compose(F(1, 10), (F(1), F(2), F(3))) == F(21, 10)
    assert certificate.compose(F(0), (F(0), F(1), F(0))) == F(1, 2)


@pytest.mark.parametrize('native_bound,work', [
    (True, (F(0),)*3), (0., (F(0),)*3), (F(-1), (F(0),)*3),
    (F(0), [F(0)]*3), (F(0), (F(0),)*2),
    (F(0), (F(0), True, F(0))), (F(0), (F(0), F(-1), F(0))),
])
def test_composition_refuses_malformed_or_negative_bounds(native_bound, work):
    with pytest.raises(ValueError):
        certificate.compose(native_bound, work)


@pytest.mark.parametrize('array', [np.array([[3., 4j], [0., 0.]]),
    np.array([[2**-1000, 0], [0, 0.]]), np.array([[3, 4], [0, 0]], dtype=np.uint64)])
def test_nuclear_bound_preserves_complex_tiny_and_integer_inputs(array):
    # Each has one nonzero row, so its exact nuclear norm is that row's norm.
    squared = sum((F(float(x.real))**2+F(float(x.imag))**2 for x in array[0]), F(0))
    upper = certificate.nuclear_upper(array)
    assert upper*upper >= squared
    assert upper*upper >= 2*squared
    assert upper > 0


@pytest.mark.parametrize('array', [np.empty((0, 0)), np.ones(3), np.ones((2, 3)),
    [[1]], np.array([[True]]), np.array([[float('nan')]]),
    np.ma.array(np.eye(2), mask=False)])
def test_nuclear_bound_refuses_invalid_original_arrays(array):
    with pytest.raises(ValueError):
        certificate.nuclear_upper(array)


def test_compiler_inputs_runtime_and_new_proof_sources_are_bound():
    from ..observer.source_scalar import compiler
    required = {compiler.RESPONSE, compiler.CONTROL_RECEIPT,
                'model/primitive_observer/runtime.py', 'model/primitive_observer/common.py',
                'model/primitive_observer/native.py', 'model/primitive_observer/exact_native.py',
                'model/primitive_observer/complete_certificate.py'}
    assert required <= set(certificate.SOURCE_PATHS)
    assert len(set(certificate.SOURCE_PATHS)) == len(certificate.SOURCE_PATHS)
    assert all('\\' not in path for path in certificate.SOURCE_PATHS)


def test_strict_portable_receipt_checks_every_scalar_leaf():
    # Run against the canonical producer artifact, not a self-fulfilling mock.
    from pathlib import Path
    from tools import strict_json
    path = Path(__file__).with_name('complete_certificate_results.json')
    if not path.exists():
        pytest.fail('canonical complete observer receipt is required')
    actual = strict_json.load(path)
    compare_receipt(actual, strict_json.loads(json.dumps(actual)))
    def leaves(value, prefix=()):
        if type(value) is dict:
            for key, item in value.items():
                yield from leaves(item, (*prefix, key))
        elif type(value) is list:
            for key, item in enumerate(value):
                yield from leaves(item, (*prefix, key))
        else:
            yield prefix
    for address in leaves(actual):
        changed = deepcopy(actual)
        row = changed
        for key in address[:-1]:
            row = row[key]
        row[address[-1]] = {'invented': True}
        with pytest.raises(ArithmeticError):
            compare_receipt(changed, actual)
