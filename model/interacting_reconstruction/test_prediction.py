"""Independent action, tensor-order and disturbed-history prediction tests."""
from itertools import product
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.interacting_reconstruction.execution import coupled_history
from model.interacting_reconstruction.prediction import RecoveredField
from model.relativistic_scalar import experiment as previous
from model.relativistic_scalar.lattice import ScalarLattice, ScalarParameters


def probabilities(rows):
    assert [row['bits'] for row in rows] == [[0, 0], [0, 1], [1, 0], [1, 1]]
    return np.array([row['probability'] for row in rows])


def independent_operator(d, omega):
    extended = d+6
    annihilation = np.zeros((extended, extended))
    for n in range(1, extended):
        annihilation[n-1, n] = math.sqrt(n)
    q = (annihilation+annihilation.T)/math.sqrt(2*omega)
    p = 1j*math.sqrt(omega/2)*(annihilation.T-annihilation)
    return q[:d, :d], (q@q)[:d, :d], (p@p)[:d, :d], np.linalg.matrix_power(q, 4)[:d, :d]


def dense_embed(dimensions, factors):
    # Entrywise occupation indexing independently checks the implementation's
    # sparse Kronecker assembly, including its little-endian tensor convention.
    states = list(product(*(range(d) for d in reversed(dimensions))))
    result = np.zeros((len(states), len(states)), dtype=complex)
    for row, reversed_left in enumerate(states):
        left = tuple(reversed(reversed_left))
        for column, reversed_right in enumerate(states):
            right = tuple(reversed(reversed_right))
            result[row, column] = math.prod(
                factors[site][left[site], right[site]] if site in factors
                else float(left[site] == right[site]) for site in range(len(dimensions)))
    return result


def dense_hamiltonian(dimensions, quadratic, quartic, bonds, omega):
    operators = [independent_operator(d, omega) for d in dimensions]
    h = np.zeros((math.prod(dimensions),)*2, dtype=complex)
    for site, (q, q2, p2, q4) in enumerate(operators):
        h += dense_embed(dimensions, {site: p2/2+quadratic[site]*q2/2+quartic[site]*q4})
    for left in range(len(dimensions)):
        for right in range(left):
            h -= bonds[left, right]*dense_embed(dimensions,
                {left: operators[left][0], right: operators[right][0]})
    return h


def dense_history(dimensions, hamiltonian, duration, read_site, sign, feedback_phase):
    u = expm(-1j*duration*hamiltonian)
    vector = np.zeros(dimensions[read_site], complex)
    vector[0], vector[4] = 1/5, sign*1j*math.sqrt(24)/5
    local_yes = np.outer(vector, vector.conj())
    yes = dense_embed(dimensions, {read_site: local_yes})
    no = np.eye(len(u))-yes
    feedback = np.eye(dimensions[read_site], dtype=complex)
    feedback[4, 4] = np.exp(1j*feedback_phase)
    feedback = dense_embed(dimensions, {read_site: feedback})
    result = []
    for first, first_kraus in enumerate((no, yes)):
        branch = first_kraus@u[:, 0]
        if first:
            branch = feedback@branch
        for last_kraus in (no, yes):
            final = last_kraus@u@branch
            result.append(float(np.vdot(final, final).real))
    return np.array(result)


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('coupling', [0., .7, 1.])
def test_predictor_agrees_with_independent_existing_exact_cube_history(sign, coupling):
    dims = (8,)+(2,)*7
    parameters = ScalarParameters(coupling=coupling)
    lattice = ScalarLattice(dims, parameters)
    bonds = -lattice.stiffness.copy()
    np.fill_diagonal(bonds, 0)
    recovered = RecoveredField(dims, np.diag(lattice.stiffness),
                               np.repeat(parameters.g, 8), bonds)
    assert np.max(abs((recovered.hamiltonian-lattice.hamiltonian).data), initial=0) < 1e-14
    actual = probabilities(recovered.history(float(previous.INTERVAL), sign=sign))
    expected = probabilities(previous.history(lattice, sign=sign, source_execution=False))
    np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)
    assert actual.min() > 0
    assert actual.sum() == pytest.approx(1., abs=3e-13)


@pytest.mark.parametrize('sign,feedback_phase', [(-1, .7), (1, -1.1), (-1, 0.)])
def test_inhomogeneous_non_cube_prediction_matches_independent_dense_circuit(sign, feedback_phase):
    dimensions = (8, 2, 2)
    quadratic, quartic = np.array([2.2, 3.1, 4.7]), np.array([.07, .04, .03])
    bonds = np.array([[0., .4, 0.], [.4, 0., .8], [0., .8, 0.]])
    field = RecoveredField(dimensions, quadratic, quartic, bonds, omega=1.7)
    independent = dense_hamiltonian(dimensions, quadratic, quartic, bonds, 1.7)
    np.testing.assert_allclose(field.hamiltonian.toarray(), independent, atol=2e-14, rtol=0)
    expected = dense_history(dimensions, independent, .17, 0, sign, feedback_phase)
    actual = probabilities(field.history(.17, sign=sign, feedback_phase=feedback_phase))
    np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)


def test_permuted_read_site_and_unequal_occupations_are_covariant():
    dimensions = (8, 2, 3)
    quadratic, quartic = np.array([3.7, 4.2, 2.1]), np.array([.05, .03, .09])
    bonds = np.array([[0., .35, .2], [.35, 0., .8], [.2, .8, 0.]])
    shuffle = np.array([1, 2, 0])
    original = RecoveredField(dimensions, quadratic, quartic, bonds, omega=1.8)
    changed = RecoveredField(tuple(dimensions[x] for x in shuffle), quadratic[shuffle],
                             quartic[shuffle], bonds[np.ix_(shuffle, shuffle)], omega=1.8)
    np.testing.assert_allclose(probabilities(original.history(.19, read_site=0)),
                               probabilities(changed.history(.19, read_site=2)), atol=4e-13, rtol=0)


def test_disturbance_feedback_and_edges_change_joint_records():
    dims = (8, 2)
    quadratic, quartic = [4., 4.], [.3, .2]
    bonds = np.array([[0., 1.7], [1.7, 0.]])
    field = RecoveredField(dims, quadratic, quartic, bonds)
    full = probabilities(field.history(.2)).reshape(2, 2)
    without_feedback = probabilities(field.history(.2, feedback_phase=0)).reshape(2, 2)
    independent_reset = np.outer(full.sum(axis=1), full.sum(axis=1))
    without_edges = probabilities(RecoveredField(dims, quadratic, quartic,
                                                 np.zeros((2, 2))).history(.2)).reshape(2, 2)
    assert .5*np.sum(abs(full-without_feedback)) > .001
    assert .5*np.sum(abs(full-independent_reset)) > .03
    assert .5*np.sum(abs(full-without_edges)) > 1e-5
    assert full.sum() == pytest.approx(1., abs=3e-13)


@pytest.mark.parametrize('read_site,sign', [(0, -1), (3, 1)])
def test_hidden_source_coupled_records_map_to_independent_cube_and_ledger(read_site, sign):
    dimensions = tuple(8 if i == read_site else 2 for i in range(8))
    actual, metadata = coupled_history(dimensions, float(previous.INTERVAL), read_site=read_site,
                                       sign=sign, steps=16)
    # Every cube vertex is equivalent for homogeneous parameters and all
    # spectator cutoffs agree.  This compares the opaque read position with
    # an independent origin computation, without importing its hidden map.
    lattice = ScalarLattice((8,)+(2,)*7)
    expected = previous.history(lattice, sign=sign, source_execution=False)
    assert previous.tv(actual, expected) < 2e-8
    assert metadata['worst_branch_resources'] == previous.source_cost(lattice, 16)
    assert metadata['all_outcomes_retained']
    assert not metadata['reset_between_reads']
    assert metadata['maximum_local_target_residual'] < 2e-10


def test_prediction_does_not_call_source_targets(monkeypatch):
    from model.interacting_reconstruction import execution, source

    def forbidden(*args, **kwargs):
        raise AssertionError('prediction leaked a backend lookup')

    monkeypatch.setattr(source, 'supplied_target', forbidden)
    monkeypatch.setattr(execution, 'target_audit', forbidden)
    monkeypatch.setattr(execution, 'coupled_history', forbidden)
    records = RecoveredField((8, 2), [4., 4.], [.04, .04], [[0., .8], [.8, 0.]]).history(.03)
    assert sum(row['probability'] for row in records) == pytest.approx(1., abs=3e-13)


@pytest.mark.parametrize('kwargs', [
    {'dimensions': (8, True)},
    {'dimensions': (8, 0)},
    {'quadratic': [4., np.nan]},
    {'quartic': [1.]},
    {'bonds': [[0., 1.], [0., 0.]]},
    {'bonds': [[1., 0.], [0., 0.]]},
    {'omega': 0.},
])
def test_invalid_reconstructed_fields_rejected(kwargs):
    arguments = dict(dimensions=(8, 2), quadratic=[4., 4.], quartic=[.1, .1],
                     bonds=[[0., 1.], [1., 0.]])
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        RecoveredField(**arguments)


@pytest.mark.parametrize('kwargs', [
    {'duration': 0.}, {'duration': -.1}, {'duration': np.nan},
    {'read_site': 1}, {'read_site': True}, {'sign': 0}, {'feedback_phase': np.inf},
])
def test_invalid_history_requests_rejected(kwargs):
    field = RecoveredField((8, 2), [4., 4.], [.1, .1], [[0., 1.], [1., 0.]])
    arguments = {'duration': .1, **kwargs}
    with pytest.raises(ValueError):
        field.history(**arguments)
