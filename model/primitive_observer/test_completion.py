"""Independent finite native-channel and strict-domain completion checks."""
from collections import defaultdict
from fractions import Fraction as F

import numpy as np
import pytest

from . import completion as c


def _endpoint(state, coin):
    """Independent floor/ceiling law; does not call native or the producer."""
    total = state[0]+state[1]
    out = list(state)
    out[:2] = ((total//2, (total+1)//2) if coin == 0 else
               ((total+1)//2, total//2))
    return tuple(out), state[0]-out[0]


def _oracle(matrix, channel='C', reference_dimension=1):
    result = defaultdict(lambda: (F(0), F(0)))
    for (i, j), value in matrix.items():
        a, r = divmod(i, reference_dimension)
        b, s = divmod(j, reference_dimension)
        for coin in (0, 1):
            left, jl = _endpoint(c.STATES[a], coin)
            right, jr = _endpoint(c.STATES[b], coin)
            equal = jl == jr if channel == 'C' else (
                c.STATES[a][0]-c.STATES[a][1] == c.STATES[b][0]-c.STATES[b][1])
            if equal:
                key = (c.INDEX[left]*reference_dimension+r,
                       c.INDEX[right]*reference_dimension+s)
                old = result[key]
                result[key] = (old[0]+value[0]/2, old[1]+value[1]/2)
    return {key: value for key, value in result.items() if value != (0, 0)}


def _dense(matrix, reference_dimension=1):
    result = np.zeros((66*reference_dimension, 66*reference_dimension), complex)
    for key, value in matrix.items():
        result[key] = complex(float(value[0]), float(value[1]))
    return result


def test_exact_attachment_checks_every_matrix_unit_and_keeps_scope():
    result = c.audit()
    assert result['native_matrix_units_checked'] == 4356
    assert result['C_eigenvalue_multiplicities'] == {'1': 2216, '1/2': 1840, '0': 300}
    assert result['retained_algebra'] == 'M46 direct-sum (I2 tensor M10)'
    assert result['scope']['all_repair_records_forgotten']
    assert not result['scope']['retained_record_history_preserved']
    assert not result['scope']['finite_step_law_selected_by_completion']


def test_every_compressed_matrix_unit_matches_independent_two_event_replay():
    for i in range(66):
        for j in range(66):
            matrix = {(i, j): (F(1), F(0))}
            assert c.apply_compressed(matrix, 2) == _oracle(_oracle(matrix))


@pytest.mark.parametrize('iterations', [1, 2, 3, 8])
def test_complex_entangled_reference_is_retained_in_complete_channel(iterations):
    a, u, _ = c.witness_indices()
    # (|a,0> + i|u,1>)/sqrt(2), with an additional unused reference direction.
    p, q = 3*a, 3*u+1
    matrix = {(p, p): (F(1, 2), F(0)), (q, q): (F(1, 2), F(0)),
              (p, q): (F(0), -F(1, 2)), (q, p): (F(0), F(1, 2))}
    wanted = matrix
    for _ in range(iterations):
        wanted = _oracle(wanted, reference_dimension=3)
    actual = c.apply_compressed(matrix, iterations, reference_dimension=3)
    assert actual == wanted
    dense = _dense(actual, 3).reshape(66, 3, 66, 3)
    reduced = np.trace(dense, axis1=0, axis2=2)
    np.testing.assert_array_equal(reduced, np.diag([.5, .5, 0]))
    # Reference-uniform numerical CP adapter must agree on the same entangled state.
    initial = _dense(matrix, 3)
    out = np.zeros_like(initial)
    for k in c.compressed_kraus(iterations):
        extended = np.kron(k, np.eye(3))
        out += extended@initial@extended.conj().T
    np.testing.assert_allclose(out, dense.reshape(198, 198), atol=2e-16, rtol=0)


@pytest.mark.parametrize('channel,count', [('C', 4), ('E', 6)])
def test_floating_kraus_adapter_has_actual_native_branch_support(channel, count):
    kraus = c.binary_kraus(channel)
    assert len(kraus) == count
    expected = {}
    for i, state in enumerate(c.STATES):
        for coin in (0, 1):
            endpoint, current = _endpoint(state, coin)
            record = (coin, current if channel == 'C' else state[0]-state[1])
            expected.setdefault(record, np.zeros((66, 66), complex))[c.INDEX[endpoint], i] = 2**-.5
    for got, key in zip(kraus, sorted(expected), strict=True):
        np.testing.assert_allclose(got, expected[key], atol=2e-16, rtol=0)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(66), atol=3e-16, rtol=0)
    kraus[0][:] = 0
    assert np.any(c.binary_kraus(channel)[0])  # no mutable shared cache


@pytest.mark.parametrize('iterations', [1, 2, 4, 10])
def test_complete_saturating_binary_measurement(iterations):
    a, u, _ = c.witness_indices()
    psi = np.zeros(66, complex)
    psi[[a, u]] = 2**-.5
    state = np.outer(psi, psi.conj())
    e_state = sum(k@state@k.conj().T for k in c.binary_kraus('E'))
    actual = state
    for _ in range(iterations):
        actual = sum(k@actual@k.conj().T for k in c.binary_kraus('C'))
    expected_plus = F(3, 8)+c.half_diamond_distance(iterations)
    plus = float(np.vdot(psi, actual@psi).real)
    e_plus = float(np.vdot(psi, e_state@psi).real)
    assert plus == pytest.approx(float(expected_plus), abs=3e-15)
    assert e_plus == pytest.approx(3/8, abs=3e-16)
    assert np.trace(actual).real == pytest.approx(1, abs=4e-15)
    difference = actual-e_state
    distance = np.sum(np.abs(np.linalg.eigvalsh(difference)))/2
    assert distance == pytest.approx(float(c.half_diamond_distance(iterations)), abs=3e-15)


def test_full_carrier_nonunital_counterexample_is_retained():
    result = c.audit()['full_carrier_counterexample']
    assert result['dimension'] == 78 and result['trace_preserved']
    assert not result['unital'] and not result['binary_tracial_completion_transfers']
    changes = {tuple(row['loads']): F(row['E_of_identity_diagonal'])
               for row in result['changed_diagonal_entries']}
    assert changes[(2, 0)+(0,)*10] == 0
    assert changes[(0, 2)+(0,)*10] == 0
    assert changes[(1, 1)+(0,)*10] == 3


@pytest.mark.parametrize('kwargs', [
    {'retain_records': True}, {'retain_records': 0}, {'retain_records': None},
    {'adaptive_seams': []}, {'adaptive_seams': [(0, 1)]},
    {'adaptive_seams': [(0, 1), (0, 5)]}, {'seam': (0, 5)},
    {'seam': [0, 1]}, {'seam': (False, 1)}, {'seam': None},
])
def test_wrong_record_or_schedule_contract_fails_before_native_execution(kwargs, monkeypatch):
    def unexpected(*args, **kw):
        pytest.fail('invalid channel contract reached the native source')
    monkeypatch.setattr(c.native, 'native', unexpected)
    for function in (lambda: c.apply_compressed({}, 1, **kwargs),
                     lambda: c.binary_kraus(**kwargs),
                     lambda: c.compressed_kraus(1, **kwargs)):
        with pytest.raises(ValueError):
            function()


@pytest.mark.parametrize('bad', [True, False, 0, -1, 1., None, '2'])
def test_invalid_event_count_rejected(bad):
    for function in (lambda: c.apply_compressed({}, bad),
                     lambda: c.compressed_kraus(bad),
                     lambda: c.half_diamond_distance(bad)):
        with pytest.raises(ValueError):
            function()


@pytest.mark.parametrize('bad', [None, [], {(0, 0): True}, {(0, 0): .5},
    {(0, 0): complex(1, 0)}, {(0, 0): float('nan')}, {(0, 0): (F(1), False)},
    {(0, 0): (F(1), F(0), F(0))}, {(False, 0): F(1)}, {(66, 0): F(1)},
    {(-1, 0): F(1)}, {(0,): F(1)}, {'00': F(1)}])
def test_sparse_operator_boundary_requires_exact_declared_entries(bad):
    with pytest.raises(ValueError):
        c.apply_compressed(bad, 1)


@pytest.mark.parametrize('bad', [True, 0, -1, 1., None])
def test_invalid_reference_dimension_rejected(bad):
    with pytest.raises(ValueError):
        c.apply_compressed({}, 1, reference_dimension=bad)


def test_exact_adapter_preserves_tiny_coherence_when_floating_weight_underflows():
    a, u, _ = c.witness_indices()
    exact = c.apply_compressed({(a, u): (1, -2)}, 1500)
    assert exact == {(a, u): (F(1, 2**1500), -F(1, 2**1499))}
    with pytest.raises(ValueError, match='unresolved'):
        c.compressed_kraus(1500)


def test_native_record_mutant_cannot_pass_the_attachment(monkeypatch):
    original = c.native.native
    def changed(state, seam, coin):
        output, _ = original(state, seam, coin)
        return output, (coin, 0)
    monkeypatch.setattr(c.native, 'native', changed)
    with pytest.raises(ArithmeticError, match='decomposition'):
        c.audit()
