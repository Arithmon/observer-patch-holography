"""The shared runtime remains inside the actual native feedback loop."""
from dataclasses import replace
from fractions import Fraction as F
import math

import numpy as np
import pytest

from ..common_source import open_native as original
from ..common_source import open_protocol, protocol, quantum_native
from ..common_source.open_environment import Environment, PhaseNoise
from . import common, native, runtime


@pytest.fixture(scope='module')
def environment():
    return Environment((PhaseNoise(1, 0, F(1, 20), 4), PhaseNoise(6, 1, F(1, 10), 4)))


@pytest.fixture(scope='module', params=('C', 'E'))
def histories(request, environment):
    return tuple(native.execute(environment, request.param, 2, 2, native_mode=mode)
                 for mode in ('protected_qubit', 'passive_qubit'))


def entangled_input():
    initial = np.array([[1, 2j, 0], [1-1j, 0, 3]], complex)
    return initial/np.linalg.norm(initial)


def test_complete_native_work_reference_outputs_and_resources_agree_with_original(histories):
    for history in histories:
        report = native.compatibility_report(history, entangled_input())
        assert report['raw_outcome_count_including_abort'] == (770 if history.law == 'C' else 2566)
        assert report['reference_dimension'] == 3
        assert report['complete_reference_output_trace_distance'] < 2e-12
        assert report['maximum_work_block_trace_distance'] < 2e-12
        assert report['source_resources_agree'] is True
        diagnostics = native.audit(history)
        assert diagnostics['prepared_work_native_input_completeness_error_diagnostic'] < 2e-10
        assert diagnostics['maximally_mixed_native_abort_mass_diagnostic'] < 2e-10
        assert diagnostics['maximally_mixed_native_padding_mass_diagnostic'] < 2e-10
        assert diagnostics['all_native_output_rows'] == 512


def test_every_record_and_raw_zero_abort_output_survives_streaming(histories):
    initial = entangled_input()
    for history in histories:
        expected = original.execute(history.environment, history.law, 2, 2, native_mode=history.native_mode)
        pairs = zip(native.output(history, initial), original.output(expected, initial), strict=True)
        reference = np.zeros((3, 3), complex)
        mass = 0.
        for (record, amplitude, density), (other_record, other_amplitude, other_density) in pairs:
            assert record == other_record
            np.testing.assert_array_equal(amplitude, other_amplitude)
            np.testing.assert_allclose(density, other_density, atol=2e-13, rtol=0)
            weight = float(np.trace(density).real)
            reference += weight*(amplitude.T@amplitude.conj())
            mass += weight*float(np.vdot(amplitude, amplitude).real)
            assert amplitude.shape == (512, 3)
            assert density.shape == (16, 16)
            assert not amplitude.flags.writeable and not density.flags.writeable
        np.testing.assert_allclose(reference, initial.T@initial.conj(), atol=2e-10, rtol=0)
        assert abs(mass-1) < 2e-10


def test_actual_current_changes_work_and_record_selects_actual_final_seam(histories):
    history = histories[1]
    workers = {current: {b.record: b.density for b in worker.branches}
               for current, worker in zip((-1, 0, 1), history.work_histories, strict=True)}
    assert max(np.linalg.norm(workers[-1][record]-workers[1][record]) for record in workers[0]) > 1e-3
    trees = {(c1, raw1): (abort, endings) for c1, raw1, _, abort, _, endings
             in quantum_native._native_tree(history.law, history.native_mode)}
    for branch in history.branches:
        c1, raw1, bits, c2, raw2 = branch.record
        if not bits:
            np.testing.assert_array_equal(branch.work_density, history.work_histories[1].initial_density)
            assert branch.abort
            continue
        current = quantum_native._record_current(history.law, c1, raw1)
        assert branch.work_density is workers[current][bits]
        _, endings = trees[c1, raw1]
        expected = next(matrix for coin, raw, _, matrix in endings[bits[2]] if (coin, raw) == (c2, raw2))
        np.testing.assert_array_equal(branch.native_map, expected)


def test_active_passive_memory_uses_same_mixed_work_and_measured_record_probability(histories):
    active, passive = histories
    witness = native.memory_witness(active, passive)
    probabilities = [math.fsum(branch.probability for branch in worker.branches if branch.record[2] == 0)
                     for worker in active.work_histories]
    qbar = (probabilities[0]+2*probabilities[1]+probabilities[2])/4
    assert abs(witness['second_smooth_plus_probability_diagnostic']-qbar) < 2e-14
    assert abs(witness['pre_read_phase_distances_diagnostic']['active']-1) < 2e-10
    assert abs(witness['pre_read_phase_distances_diagnostic']['passive']-(1-qbar)) < 2e-10
    assert abs(witness['active_encoded_bell_fidelity_diagnostic']-1) < 2e-10
    assert witness['native_storage_or_record_fault_tolerance_certified'] is False


def test_producer_never_uses_old_mixed_or_open_native_trajectory(monkeypatch, environment):
    def forbidden(*args, **kwargs):
        raise AssertionError('older trajectory cannot produce shared-runtime observer')
    monkeypatch.setattr(original, 'execute', forbidden)
    monkeypatch.setattr(open_protocol, 'execute', forbidden)
    history = native.execute(environment, 'C', 2, 2)
    assert len(history.work_histories) == 3
    assert all(type(worker) is runtime.History for worker in history.work_histories)
    assert len(tuple(native.output(history, entangled_input()))) == 770


def test_correlated_native_work_input_and_untyped_inputs_are_rejected(histories):
    history = histories[0]
    invalid = (np.zeros((2, 16, 3)), [[True, 0], [0, 0]],
               np.ma.array(np.eye(2)/math.sqrt(2), mask=False), np.ones((2, 3)),
               np.array([[float('nan')], [0.]]))
    for initial in invalid:
        with pytest.raises(ValueError):
            native.output(history, initial)


def test_resealed_work_replacement_cannot_change_declared_independence_contract(histories):
    history = histories[0]
    program = protocol.continuation_program(2, 2, native_current=-1)
    compiled = common.compile_program(program, history.environment)
    vector = np.zeros(compiled.work_dimension*2)
    vector[0] = 1
    supplied = runtime.execute(compiled, initial=np.outer(vector, vector))
    corrupt = replace(history, work_histories=(supplied, *history.work_histories[1:]), seal='')
    corrupt = replace(corrupt, seal=native._seal(corrupt))
    with pytest.raises(ValueError, match='independently prepared'):
        native.validate_history(corrupt)


def test_resealed_omission_relabelled_record_and_detached_native_map_are_rejected(histories):
    history = histories[0]
    branch = next(row for row in history.branches if row.record[2] and np.linalg.norm(row.native_map) > .01)
    mutated = replace(branch, native_map=branch.native_map*0)
    for branches in (history.branches[:-1], (history.branches[1], history.branches[0], *history.branches[2:]),
                     tuple(mutated if row is branch else row for row in history.branches)):
        corrupt = replace(history, branches=branches, seal='')
        corrupt = replace(corrupt, seal=native._seal(corrupt))
        with pytest.raises(ValueError):
            native.validate_history(corrupt)


def test_lazy_output_is_detached_from_caller_envelopes(histories):
    original_history = histories[0]
    history = replace(original_history, branches=tuple(replace(row) for row in original_history.branches))
    initial = entangled_input()
    stream = native.output(history, initial)
    expected = next(native.output(original_history, initial))
    initial[:] = 0
    history.branches[0].__dict__['work_density'] = np.zeros((16, 16))
    history.branches[0].__dict__['record'] = ('forged',)
    actual = next(stream)
    assert actual[0] == expected[0]
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])


def test_full_carrier_retains_all_78_input_columns_and_512_output_rows(environment):
    history = native.execute(environment, 'C', 2, 1, native_mode='full_carrier')
    assert all(branch.native_map.shape == (512, 78) for branch in history.branches)
    report = native.audit(history)
    assert report['prepared_work_native_input_completeness_error_diagnostic'] < 2e-10
    assert native.compatibility_report(history)['every_native_map_and_raw_record_agrees'] is True
