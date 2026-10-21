"""Independent cold-history and reference-sensitive refutation controls.

Dense finite oracles are numerical cross-checks of the outward residual
certificate, not replacements for its infinite-occupation proof.
"""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import thermal_history as thermal


def independent_operators(d, coupling=1.):
    # Form powers before compression in a padded orthonormal oscillator.
    padded = d+5
    annihilation = np.diag(np.sqrt(np.arange(1, padded)), 1)
    q = (annihilation+annihilation.T)/2
    p = 1j*(annihilation.T-annihilation)
    q2, q4, p2 = q@q, np.linalg.matrix_power(q, 4), p@p
    first = (p2/2+q2/2+coupling*q4/16)[:d, :d]
    second = (p2/2+2*q2+coupling*q4/16)[:d, :d]
    hamiltonian = (np.kron(np.eye(d), first)+np.kron(second, np.eye(d))
                   +coupling*np.kron(q2[:d, :d], q2[:d, :d])/4)
    return hamiltonian, np.kron(q[:d, :d], np.eye(d))


def finite_history(state, sign, hamiltonian, coordinate, *, interval=1/16,
                   return_prefixes=False):
    values, basis = np.linalg.eigh(coordinate)
    reads = tuple((basis*function((values+np.pi/2)/2))@basis.conj().T
                  for function in (np.cos, np.sin))
    evolve = expm(-1j*interval*hamiltonian)
    kick = expm(1j*sign*coordinate/4)
    feedback = expm(1j*coordinate/4)
    prepared = kick@state
    first = [read@evolve@prepared for read in reads]
    first[1] = feedback@first[1]
    final = tuple(read@evolve@branch for branch in first for read in reads)
    return ((prepared,), tuple(first), final) if return_prefixes else final


def pure_distance(left, right):
    # Independent stable rank-two trace norm, sensitive to reference columns.
    x, y = left.reshape(-1), right.reshape(-1)
    xx, yy = np.vdot(x, x).real, np.vdot(y, y).real
    if xx == 0:
        return float(yy/2)
    perpendicular = y-x*(np.vdot(x, y)/xx)
    return float(np.sqrt((xx-yy)**2+4*xx*np.vdot(perpendicular, perpendicular).real)/2)


def complete_distance(left, right):
    return sum(pure_distance(a, b) for a, b in zip(left, right, strict=True))


@pytest.fixture(scope='module')
def cold():
    d, beta = 8, 4
    h, q = independent_operators(d)
    energies, basis = np.linalg.eigh(h)
    weights = np.exp(-beta*(energies-energies[0])); weights /= weights.sum()
    canonical = (basis*np.sqrt(weights))@basis.conj().T
    seed, seed_budget = thermal.finite_gibbs_seed(d, F(beta))
    rows = {sign: thermal.history(seed, sign, occupation=d, cells=8)
            for sign in (-1, 1)}
    return {'d': d, 'h': h, 'q': q, 'basis': basis, 'weights': weights,
            'canonical': canonical, 'seed': seed, 'seed_budget': seed_budget,
            'rows': rows}


def test_cold_seed_matches_independent_canonical_matrix_with_its_vector_bound(cold):
    error = np.linalg.norm(cold['seed']-cold['canonical'])
    assert error <= float(cold['seed_budget']['canonical_seed_vector_error_upper'])+2e-13
    assert cold['seed_budget']['reference_convention'] == 'common frequency-two Hermite basis'
    assert not cold['seed_budget']['eigensolver_accuracy_assumed']
    assert np.linalg.norm(cold['seed']-cold['seed'].conj().T) < 2e-15


@pytest.mark.parametrize('sign', (-1, 1))
def test_cold_complete_four_record_history_matches_independent_finite_oracle(cold, sign):
    actual, budget = cold['rows'][sign]
    expected = finite_history(cold['seed'], sign, cold['h'], cold['q'])
    # The same finite input and finite generator remove occupation error from
    # this comparison. All four vectors include every reference column.
    allowance = budget['finite_residual']+budget['endpoint_jump']
    vector_error = np.sqrt(sum(np.linalg.norm(a-b)**2 for a, b in zip(actual, expected, strict=True)))
    assert vector_error <= allowance+3e-12
    assert budget['all_record_branches'] == 4
    assert budget['retained_reference_columns'] == 64
    assert [state.shape for state in actual] == [(64, 64)]*4
    assert abs(sum(np.linalg.norm(state)**2 for state in actual)-np.linalg.norm(cold['seed'])**2) < 1e-10
    assert budget['occupation_residual'] > 0
    assert budget['complete_output_trace_error_upper'] < .002


def test_wrong_reference_labels_preserve_work_density_but_fail_complete_quantum_comparison(cold):
    wrong = cold['basis']@np.diag(np.sqrt(cold['weights']))
    np.testing.assert_allclose(wrong@wrong.conj().T,
                               cold['canonical']@cold['canonical'].conj().T,
                               atol=2e-15, rtol=0)
    expected = finite_history(cold['canonical'], 1, cold['h'], cold['q'])
    wrong_outputs = finite_history(wrong, 1, cold['h'], cold['q'])
    np.testing.assert_allclose([np.linalg.norm(state)**2 for state in wrong_outputs],
                               [np.linalg.norm(state)**2 for state in expected],
                               atol=2e-15, rtol=0)
    assert complete_distance(expected, wrong_outputs) > .01
    assert complete_distance(expected, wrong_outputs) > 10*cold['rows'][1][1]['complete_output_trace_error_upper']


def test_erasing_a_branch_and_reversing_the_preparation_sign_are_resolved(cold):
    plus, budget = cold['rows'][1]
    minus = cold['rows'][-1][0]
    erased = list(plus); erased[2] = np.zeros_like(erased[2])
    assert complete_distance(plus, erased) > .01
    assert complete_distance(plus, minus) > .05
    assert complete_distance(plus, erased) > 10*budget['complete_output_trace_error_upper']
    assert sum(np.linalg.norm(state)**2 for state in erased) < .9


def test_uniform_occupation_bound_covers_stopped_prefixes_and_classical_timer_mixture():
    # A larger dense truncation is an independent numerical cross-check of
    # the infinite-space certificate, not a replacement for its proof.
    d, larger, beta = 4, 12, 2
    h, q = independent_operators(d)
    larger_h, larger_q = independent_operators(larger)

    def canonical(hamiltonian):
        values, vectors = np.linalg.eigh(hamiltonian)
        weights = np.exp(-beta*(values-values[0])); weights /= weights.sum()
        return (vectors*np.sqrt(weights))@vectors.conj().T

    finite_input, larger_input = canonical(h), canonical(larger_h)
    seed, seed_budget = thermal.finite_gibbs_seed(d, F(beta))
    preparation = thermal.preparation_bracket(d, 12, F(beta))
    ids = [i+larger*j for j in range(d) for i in range(d)]
    weights = (F(25, 36), F(10, 36), F(1, 36))
    assert sum(weights) == 1
    for sign in (-1, 1):
        _, budget = thermal.history(seed, sign, occupation=d)
        finite_error = F(budget['finite_residual'])+F(budget['endpoint_jump'])
        bound = min(F(1), preparation['canonical_preparation_trace_distance_upper']
                    +2*seed_budget['canonical_seed_vector_error_upper']
                    +2*finite_error+F(budget['occupation_residual']))
        finite = finite_history(finite_input, sign, h, q, return_prefixes=True)
        expected = finite_history(larger_input, sign, larger_h, larger_q,
                                  return_prefixes=True)
        assert [len(branches) for branches in finite] == [1, 2, 4]
        distances, weighted_finite, weighted_expected = [], [], []
        for branches, targets, weight in zip(finite, expected, weights, strict=True):
            lifted = []
            for branch in branches:
                # Embed both work and untouched reference in the common
                # Hermite basis; never compare work marginals alone.
                value = np.zeros((larger*larger, larger*larger), complex)
                value[np.ix_(ids, ids)] = branch
                lifted.append(value)
            assert sum(np.linalg.norm(value)**2 for value in lifted) == pytest.approx(1, abs=3e-13)
            assert sum(np.linalg.norm(value)**2 for value in targets) == pytest.approx(1, abs=3e-13)
            distances.append(complete_distance(lifted, targets))
            assert distances[-1] < float(bound)
            amplitude = np.sqrt(float(weight))
            weighted_finite.extend(amplitude*value for value in lifted)
            weighted_expected.extend(amplitude*value for value in targets)
        # Distinct timer and field labels are distinct output blocks. A
        # classical mixture weights complete distances, without conditioning
        # on completing two reads or dividing by the rare two-tick weight.
        mixture = complete_distance(weighted_finite, weighted_expected)
        average = sum(float(weight)*distance for weight, distance in
                      zip(weights, distances, strict=True))
        assert mixture == pytest.approx(average, abs=3e-13)
        assert mixture <= max(distances)+3e-13 < float(bound)


@pytest.mark.parametrize('kind,coefficient', (('coordinate', F(-1, 4)),
                                            ('coordinate', F(1, 4)),
                                            ('hamiltonian', F(-1, 16)),
                                            ('hamiltonian', F(1, 16))))
def test_signed_generators_agree_with_independent_fixed_basis_exponentials(kind, coefficient):
    d = 4
    h, q = independent_operators(d)
    rng = np.random.default_rng(883)
    state = rng.normal(size=(d*d, 3))+1j*rng.normal(size=(d*d, 3))
    state /= np.linalg.norm(state)
    generator = -1j*float(coefficient)*h if kind == 'hamiltonian' else 1j*float(coefficient)*q
    actual, budget = thermal.evolve(state, kind, coefficient, occupation=d, cells=8)
    expected = expm(generator)@state
    assert np.linalg.norm(actual-expected) <= budget['finite_residual']+budget['endpoint_jump']+2e-12
    assert np.linalg.norm(actual-expm(-generator)@state) > .05


@pytest.mark.parametrize('bad', (True, .25, np.float64(.25), '1/4'))
def test_signed_generator_parameters_reject_implicit_scalar_types(bad):
    with pytest.raises(ValueError):
        thermal.generator(4, 'coordinate', bad)
    with pytest.raises(ValueError):
        thermal.evolve(np.eye(16, dtype=complex), 'hamiltonian', bad, occupation=4)


@pytest.mark.parametrize('bad', (True, 0, -1, .25, np.float64(.25), '1/4'))
def test_seed_inverse_temperature_requires_positive_exact_input(bad):
    with pytest.raises(ValueError):
        thermal.finite_gibbs_seed(4, bad)


@pytest.mark.parametrize('mutation', ('nan', 'infinity', 'vector', 'integer', 'empty_reference'))
def test_complete_history_rejects_invalid_amplitude_domains(mutation):
    initial = np.eye(16, dtype=complex)/4
    if mutation == 'nan':
        initial[0, 0] = np.nan
    elif mutation == 'infinity':
        initial[0, 0] = np.inf
    elif mutation == 'vector':
        initial = initial[:, 0]
    elif mutation == 'integer':
        initial = np.eye(16, dtype=int)
    else:
        initial = np.zeros((16, 0), dtype=complex)
    with pytest.raises(ValueError):
        thermal.history(initial, occupation=4)
