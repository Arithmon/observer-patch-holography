"""Full-state and continuation controls independent of the screen dual."""
from fractions import Fraction as F
import hashlib
import json

import numpy as np
import pytest

from ..observer import native
from . import oracle as o


EDGES = ((0, 1), (0, 2), (1, 2), (2, 3), (3, 4))
INCIDENT = EDGES[:4]


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('cap', (1, 2, 4))
def test_complete_records_mass_and_reference_coherence(law, cap):
    outputs = o.exact_history(5, EDGES, incident=INCIDENT, law=law, cap=cap)
    assert tuple(outputs) == o.RECORDS and len(outputs) == 21
    q = F(1, 5)**cap
    hit = (1-q)/4
    chi = (1, F(1, 2), F(1, 2), 0) if law == 'C' else (1, 0, 0, 0)
    for (status, word), matrix in outputs.items():
        assert matrix.shape == (8, 8)
        assert all(type(x) is F for x in matrix.flat)
        expected_mass = hit**len(word)*(q if status == 'timeout' else 1)
        coherence = expected_mass/2
        for mask in word:
            coherence *= chi[mask]
        if word and word[0] == 1:
            coherence *= -1
        assert np.trace(matrix) == expected_mass
        assert matrix[2, 5] == matrix[5, 2] == coherence
        assert all(matrix[i, j] == 0 for i in range(8) for j in range(8)
                   if i != j and (i, j) not in ((2, 5), (5, 2)))
    assert sum(np.trace(matrix) for matrix in outputs.values()) == 1
    assert sum(np.trace(v) for (status, _), v in outputs.items() if status == 'timeout') == 2*q-q*q


@pytest.mark.parametrize('law', ('C', 'E'))
def test_remote_wait_moves_the_actual_full_state_and_keeps_its_mass(law):
    edges = ((0, 1), (0, 2), (1, 2), (2, 3))
    # One particle at remote port2; remote attempt moves it to2 or3.
    _, timeout = o._stopped({(8, 8): F(1)}, edges, edges[:3], law, 1)
    assert timeout == {(8, 8): F(1, 8), (16, 16): F(1, 8)}
    assert sum(timeout.values()) == F(1, 4)  # Neither idle replacement nor renormalization.


@pytest.mark.parametrize('law', ('C', 'E'))
def test_mask_zero_attempt_is_a_real_clock_tick_even_when_no_load_moves(law):
    outputs = o.exact_history(5, EDGES, incident=INCIDENT, law=law, cap=1, rho=F(0))
    assert np.trace(outputs['complete', (0, 0)]) == F(1, 25)
    assert np.trace(outputs['timeout', ()]) == F(1, 5)


def test_feedback_changes_phase_without_changing_any_record_probability():
    controlled = o.exact_history(5, EDGES, incident=INCIDENT, cap=2)
    ordinary = o.exact_history(5, EDGES, incident=INCIDENT, cap=2, feedback=False)
    assert controlled['complete', (1, 0)][2, 5] == -F(9, 625)
    for record in o.RECORDS:
        assert np.array_equal(np.diag(controlled[record]), np.diag(ordinary[record]))
        sign = -1 if record[1] and record[1][0] == 1 else 1
        assert controlled[record][2, 5] == sign*ordinary[record][2, 5]


def test_retained_C_reference_block_is_entangled_and_E_control_is_separable():
    for law, witness in (('C', -F(9, 625)), ('E', F(0))):
        matrix = o.exact_history(5, EDGES, incident=INCIDENT, cap=2, law=law)['complete', (1, 0)]
        # Unit witness (|pair1,ref1>+|pair2,ref0>)/sqrt2 after reference transpose.
        # No conditioning or division by the record's probability is performed.
        expectation = (matrix[3, 3]+matrix[4, 4]+matrix[2, 5]+matrix[5, 2])/2
        assert expectation == witness
        if law == 'E':
            assert all(matrix[i, j] == 0 for i in range(8) for j in range(8) if i != j)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_full_runtime_checkpoint_and_all_reduced_reference_outputs(law):
    exact = o.exact_history(5, EDGES, incident=INCIDENT, law=law, cap=2, rho=F(1, 5))
    actual, receipt = o.execute_history(5, EDGES, incident=INCIDENT, law=law, cap=2, rho=F(1, 5))
    assert receipt['full_output_distance_numeric'] < 1e-12
    assert receipt['reduced_output_distance_numeric'] < 1e-12
    assert receipt['checkpoint_distance'] == receipt['zero_padding_nuclear_norm'] == 0
    assert receipt['complete_runtime_records'] == 25 and receipt['comparison_records'] == 21
    assert receipt['checkpoint_sha256'] == hashlib.sha256(receipt['checkpoint_json'].encode()).hexdigest()
    assert receipt['checkpoint_bytes'] == len(receipt['checkpoint_json'].encode())
    assert receipt['maximum_branch_resources'] == {
        'native_attempts_upper': 4, 'stopped_intervals': 2, 'logical_phase_controls': 1}
    assert receipt['expanded_tree_resources'] == {
        'native_attempts_upper': 10, 'stopped_intervals': 5, 'logical_phase_controls': 1}
    assert receipt['expected_native_attempts'] == F(294, 125)
    assert receipt['timeout_probability'] == F(49, 625)
    assert receipt['declared_joint_preparations'] == 1
    assert receipt['physical_pulse_execution'] is False
    for record in o.RECORDS:
        np.testing.assert_allclose(actual[record], np.array(exact[record], complex), atol=1e-15, rtol=0)
    checkpoint = json.loads(receipt['checkpoint_json'])
    assert checkpoint['cursor'] == 1 and checkpoint['reference_dimension'] == 2
    assert len(checkpoint['branches']) == 5


def test_runtime_checkpoint_tampering_is_not_accepted(monkeypatch):
    original = o.runtime.checkpoint_to_json
    def change_record(history):
        value = json.loads(original(history))
        value['branches'][0]['record'] = [1]
        return json.dumps(value)
    monkeypatch.setattr(o.runtime, 'checkpoint_to_json', change_record)
    with pytest.raises(ValueError):
        o.execute_history(2, ((0, 1),), cap=1)


@pytest.mark.parametrize('kwargs', (
    {'n': True}, {'n': 6}, {'edges': list(EDGES)}, {'edges': ((True, 2),)},
    {'edges': ((0, 1), (0, 1))}, {'edges': ((1, 0),)}, {'root': (0, True)},
    {'root': (0,)}, {'incident': ((3, 4),)}, {'incident': ((0, True),)},
    {'cap': True}, {'cap': 0}, {'cap': 5}, {'feedback': 1}, {'law': 1},
    {'law': 'D'}, {'rho': 0.5}, {'rho': F(-1)}, {'rho': F(3, 2)},
))
def test_original_parameter_types_and_clock_geometry_are_checked(kwargs):
    arguments = {'n': 5, 'edges': EDGES, 'incident': INCIDENT, **kwargs}
    with pytest.raises(ValueError):
        o.exact_history(**arguments)


def test_corrupted_fraction_is_rejected_and_returned_fractions_do_not_alias_inputs():
    invalid = F(1, 3)
    invalid._numerator, invalid._denominator = -1, -3
    with pytest.raises(ValueError):
        o.exact_history(5, EDGES, incident=INCIDENT, rho=invalid)
    population = F(1, 5)
    outputs = o.exact_history(5, EDGES, incident=INCIDENT, cap=1, rho=population)
    saved = tuple(tuple((x.numerator, x.denominator) for x in row.flat) for row in outputs.values())
    population._numerator = 2
    assert saved == tuple(tuple((x.numerator, x.denominator) for x in row.flat) for row in outputs.values())


def test_receipt_rechecks_original_source_endpoint_before_claiming_execution(monkeypatch):
    def wrong_native(loads, seam, coin):
        return loads, (coin, 0)
    monkeypatch.setattr(native, 'native', wrong_native)
    with pytest.raises(ArithmeticError, match='original native endpoint'):
        o.run()


def test_run_exports_replayed_structure_without_platform_specific_hashes():
    receipt = o.run()
    assert receipt['original_native_endpoint_checks'] == 8
    for row in receipt['executions']:
        assert row['checkpoint_replayed'] is True
        assert not {'checkpoint_json', 'checkpoint_bytes', 'checkpoint_sha256'} & row.keys()
        assert abs(row['probability_mass']-1) < 1e-12


def test_runtime_rejects_wrong_reduced_reference_state_with_correct_full_output(monkeypatch):
    original = np.trace
    phase = np.array([1., -1.]*4)

    def wrong_reduction(value, *args, **kwargs):
        result = original(value, *args, **kwargs)
        if (getattr(value, 'ndim', 0) == 4 and kwargs.get('axis1') == 0
                and kwargs.get('axis2') == 2):
            # Positive, same-trace output with the opposite reference phase.
            return result*phase[:, None]*phase[None, :]
        return result

    monkeypatch.setattr(o.np, 'trace', wrong_reduction)
    with pytest.raises(ArithmeticError, match='full exact stopped history'):
        o.execute_history(5, EDGES, incident=INCIDENT, cap=2)


@pytest.mark.parametrize('bad', (float('nan'), float('inf'), -1e-10))
@pytest.mark.parametrize('which', (0, 1))
def test_runtime_rejects_invalid_full_or_reduced_distance(monkeypatch, bad, which):
    distances = [0., 0.]
    distances[which] = bad
    monkeypatch.setattr(o, '_distance', lambda *_: distances.pop(0))
    with pytest.raises(ArithmeticError, match='full exact stopped history'):
        o.execute_history(2, ((0, 1),), cap=1)
