"""Independent complete-reference and finite-feedback controls."""
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import native_response_dynamics as m


def test_full_cost_is_the_ray_disturbance_of_all_eighteen_generators():
    p = F(1, 25)
    ray = (m._bell(0)+s.sqrt(24)*m._bell(3))/5
    columns = s.Matrix.hstack(*((s.kronecker_product(g, s.eye(6))
        +s.kronecker_product(s.eye(6), g))*ray for g in m.full_basis()))
    actual = m._clean((columns.conjugate().T*(s.eye(36)-ray*ray.conjugate().T)*columns).applyfunc(s.re))
    assert actual == m.cost_gram(p)
    assert actual.eigenvals() == {0: 7, s.Rational(4, 75): 5,
                                  s.Rational(32, 25): 5, s.Rational(64, 625): 1}


def test_limiting_range_is_closed_but_the_projection_is_not_a_quotient():
    report = m.canonical_certificate()
    assert report['limit_raw_dimension'] == 13
    assert report['limit_observable_dimension'] == 12
    assert report['limit_range_bracket_closed']
    assert not report['projection_is_a_Lie_homomorphism']
    assert report['finite_interior_kernel_observable_dimension'] == 6
    assert not report['ordinary_reference_cooling_alone_implements_feedback']
    assert report['original17_controls_limit_without_work_class_sum_phase']['observable_dimension'] == 11


def test_unbounded_input_amplification_restores_a_suppressed_finite_current():
    delta = 64*s.sqrt(5)
    p = (1-delta/s.sqrt(delta**2+s.Rational(2, 45)**2))/2
    feedback = m.feedback_matrix(p, 16)
    gain = m.feedback_bound(p, 16)['high_current_gain']
    input_coefficients = s.eye(18)[:, 12]/gain
    assert feedback*input_coefficients == s.eye(18)[:, 12]
    assert float(1/gain) > 4e7
    assert m.limit_projection()*s.eye(18)[:, 12] == s.zeros(18, 1)


@pytest.mark.parametrize('p,steps,step', [(F(1, 25), 0, F(1, 2)),
    (F(1, 25), 3, F(1, 2)), (F(1, 4), 4, F(3, 4)), (F(1, 2), 2, F(1, 3)),
    (F(0), 2, F(3, 4)), (F(1), 3, F(1, 2))])
def test_functional_filter_equals_every_literal_feedback_update(p, steps, step):
    literal = (s.eye(18)-s.Rational(step.numerator, step.denominator)*m.cost_gram(p))**steps
    assert m.feedback_matrix(p, steps, step) == literal


@pytest.mark.parametrize('p,steps', [(F(1, 25), 3), (F(1, 4), 5), (F(0), 2), (F(1, 2), 4)])
def test_exact_filter_distance_matches_full_spectrum(p, steps):
    delta = m.feedback_matrix(p, steps)-m.limit_projection()
    exact = max(abs(x) for x in delta.eigenvals())
    bounds = m.feedback_bound(p, steps)
    assert exact == bounds['coefficient_operator_error']
    assert bounds['low_current_loss'] <= bounds['low_current_loss_upper_bound']
    assert bounds['high_current_gain'] <= bounds['high_current_gain_upper_bound']


def test_source_bias_joint_regime_has_finite_small_error_without_a_cutoff():
    delta = 64*s.sqrt(5)
    p = (1-delta/s.sqrt(delta**2+s.Rational(2, 45)**2))/2
    bound = m.feedback_bound(p, 16)
    assert 0 < float(bound['coefficient_operator_error']) < 5.2e-7
    assert float(bound['high_current_gain']) < 2.4e-8
    assert np.linalg.matrix_rank(np.asarray(m.feedback_matrix(p, 16), float)) == 18
    # A fixed finite reference eventually damps both nominally slow sectors.
    zero, first, second, relative = m._spectral_data()
    assert zero.rank() == 7 and m.limit_projection().rank() == 13
    assert first.rank() == 5 and second.rank() == 5 and relative.rank() == 1


def test_continuous_flow_bound_matches_the_complete_matrix_exponential():
    p, time = F(1, 25), F(7, 2)
    actual = expm(-float(time)*np.asarray(m.cost_gram(p), float))
    distance = np.linalg.norm(actual-np.asarray(m.limit_projection(), float), 2)
    bounds = m.continuous_bound(p, time)
    assert abs(distance-float(bounds['coefficient_operator_error'])) < 1e-14
    assert float(bounds['low_current_loss']) <= float(bounds['low_current_loss_upper_bound'])


def test_entire_limiting_group_has_uniform_reference_cost_and_relative_center_attains_it():
    p = F(1, 25)
    ray = np.asarray((m._bell(0)+s.sqrt(24)*m._bell(3))/5, complex)[:, 0]
    center = np.diag([1j, 1j, 1j, 0, 0, 0])
    unitary = expm(np.pi*center/2)
    failure = 1-abs(np.vdot(ray, np.kron(unitary, unitary)@ray))**2
    assert abs(failure-4*float(p)*(1-float(p))) < 1e-14
    rng = np.random.default_rng(982)
    x = np.asarray(m.limit_projection(), float)@rng.normal(size=18)
    X = np.einsum('i,ijk->jk', x, np.asarray(m.full_basis(), complex))
    U = expm(5*X)
    assert 1-abs(np.vdot(ray, np.kron(U, U)@ray))**2 <= 4*float(p)*(1-float(p))+1e-14


@pytest.mark.parametrize('p,expected', [(0, 0), (F(1, 4), F(3, 4)),
    (F(1, 2), 1), (F(3, 4), 1), (1, 1)])
def test_global_limiting_group_bound_keeps_both_weight_halves(p, expected):
    assert m.limiting_group_failure_bound(p) == s.Rational(expected)


def test_upper_half_reference_can_fail_completely_under_a_limiting_group_control():
    # Construct the full reference independently of the producer's Bell helper.
    columns = s.eye(6)
    A = sum((s.kronecker_product(columns[:, j], columns[:, j])
             for j in range(3)), s.zeros(36, 1))/s.sqrt(3)
    B = sum((s.kronecker_product(columns[:, j], columns[:, j])
             for j in range(3, 6)), s.zeros(36, 1))/s.sqrt(3)
    p = s.Rational(3, 4)
    ray = s.sqrt(p)*A+s.sqrt(1-p)*B
    U = s.diag(s.I, s.I, 1, 1, 1, 1)
    overlap = s.simplify((ray.conjugate().T*s.kronecker_product(U, U)*ray)[0])
    failure = 1-overlap*s.conjugate(overlap)
    assert failure == m.limiting_group_failure_bound(p) == 1
    assert failure > 4*p*(1-p)  # The old unqualified certificate claimed 3/4.
    # At p=1 the global bound remains one, rather than the false zero.
    z = (1+s.I*s.sqrt(3))/2
    U = s.diag(1, z, s.conjugate(z), 1, 1, 1)
    assert s.simplify((A.conjugate().T*s.kronecker_product(U, U)*A)[0]) == 0


def test_maximum_step_is_invertible_only_at_interior_weights():
    assert m.feedback_matrix(F(1, 4), 1, F(3, 4)).rank() == 18
    for p in (0, 1):
        assert m.feedback_matrix(p, 0, F(3, 4)).rank() == 18
        assert m.feedback_matrix(p, 1, F(3, 4)).rank() == 13


def test_actual_finite_reference_failure_probability_has_the_computed_hessian():
    p = F(1, 25)
    rng = np.random.default_rng(6501)
    coefficients = rng.normal(size=18)
    coefficients /= np.linalg.norm(coefficients)
    X = sum((a*np.asarray(g, complex) for a, g in zip(coefficients, m.full_basis(), strict=True)),
            np.zeros((6, 6), complex))
    ray = np.asarray((m._bell(0)+s.sqrt(24)*m._bell(3))/5, complex)[:, 0]
    t = 2.**-12
    U = expm(t*X)
    moved = np.kron(U, U)@ray
    failure = np.linalg.norm(moved-ray*np.vdot(ray, moved))**2
    cost = coefficients@np.asarray(m.cost_gram(p), float)@coefficients
    assert abs(failure/t**2-cost) < 2e-7


def test_complete_channel_bound_holds_on_an_entangled_reference_and_both_outputs():
    p, steps, duration = F(1, 25), 9, .3
    rng = np.random.default_rng(913)
    coefficients = rng.normal(size=18)
    coefficients /= np.linalg.norm(coefficients)
    Fm = np.asarray(m.feedback_matrix(p, steps), float)
    P = np.asarray(m.limit_projection(), float)
    basis = np.asarray(m.full_basis(), complex)
    U = expm(duration*np.einsum('i,ijk->jk', Fm@coefficients, basis))
    V = expm(duration*np.einsum('i,ijk->jk', P@coefficients, basis))
    # A maximally entangled work/reference state; retain both detector records.
    state = np.eye(6).reshape(36)/np.sqrt(6)
    actual = np.kron(U, np.eye(6))@state
    target = np.kron(V, np.eye(6))@state
    effects = (np.diag([1, 1, 0, 0, 0, 0]), np.diag([0, 0, 1, 1, 1, 1]))
    distance = 0.
    for effect in effects:
        a, b = np.kron(effect, np.eye(6))@actual, np.kron(effect, np.eye(6))@target
        difference = np.outer(a, a.conj())-np.outer(b, b.conj())
        distance += np.linalg.eigvalsh(difference).__abs__().sum()
    bound = 2*duration*float(m.feedback_bound(p, steps)['coefficient_operator_error'])
    assert distance <= bound+1e-13


def test_finite_current_bracket_leakage_obeys_high_gain_bound():
    rng = np.random.default_rng(781)
    p, steps = F(1, 25), 5
    basis = np.asarray(m.full_basis(), complex)
    Fm, P = np.asarray(m.feedback_matrix(p, steps), float), np.asarray(m.limit_projection(), float)
    x, y = rng.normal(size=(2, 18))
    X, Y = (np.einsum('i,ijk->jk', Fm@row, basis) for row in (x, y))
    coordinates = np.real(np.einsum('ijk,kj->i', basis.conj(), (X@Y-Y@X).T))
    leak = np.linalg.norm((np.eye(18)-P)@coordinates)
    upper = 2*float(m.feedback_bound(p, steps)['high_current_gain'])*np.linalg.norm(x)*np.linalg.norm(y)
    assert leak <= upper+1e-13


@pytest.mark.parametrize('p', [True, False, .1, np.float64(.1), '1/4', s.Symbol('p'), s.oo, -1, 2])
def test_invalid_reference_weights_are_not_silently_coerced(p):
    with pytest.raises(ValueError):
        m.cost_gram(p)


@pytest.mark.parametrize('steps,step', [(True, F(1, 2)), (-1, F(1, 2)),
    (1., F(1, 2)), (1, True), (1, 0), (1, 1), (1, .5)])
def test_invalid_feedback_controls_are_refused(steps, step):
    with pytest.raises(ValueError):
        m.feedback_matrix(F(1, 25), steps, step)


def test_frames_and_outputs_cannot_be_mutated():
    with pytest.raises(TypeError):
        m.full_basis()[0][0, 0] = 1
    with pytest.raises(TypeError):
        m.feedback_matrix(F(1, 25), 3)[0, 0] = 1
    record = m.feedback_bound(F(1, 25), 3)
    record['high_current_gain'] = 10
    assert m.feedback_bound(F(1, 25), 3)['high_current_gain'] < 1


@pytest.fixture(scope='module')
def source_certificate():
    return m.source_certificate()


def test_complete_native_source_and_actual_biased_ground_bind_the_cost(source_certificate):
    assert source_certificate['complete_native18_bidirectional_span_checked']
    assert source_certificate['actual_full36_ground_disturbance_Gram_checked']
    assert source_certificate['limit_common_phase_quotient_matches_declared12_as_validation_only']
    assert json.loads(json.dumps(source_certificate, allow_nan=False)) == source_certificate


def test_a_false_ground_weight_cannot_replace_the_actual_full_reference(source_certificate, monkeypatch):
    from . import native_response_biased as biased
    original = biased.exact_data(F(16))
    changed = dict(original)
    changed['first_triplet_weight'] = s.Rational(1, 4)
    monkeypatch.setattr(biased, 'exact_data', lambda bias: changed)
    with pytest.raises(ArithmeticError, match='disturbance Gram'):
        m.source_certificate()
