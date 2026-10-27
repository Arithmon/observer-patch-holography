"""Independent exact action, invariant and operational controls for17versus12."""
import json

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import native_response_comparison as m


@pytest.fixture(scope='module')
def data():
    return m.exact_data()


@pytest.fixture(scope='module')
def certificate():
    return m.exact_certificate()


def rank(matrices):
    columns = s.Matrix.hstack(*(s.Matrix([s.re(x).expand() for x in a]+[s.im(x).expand() for x in a])
                               for a in matrices))
    return columns.to_DM(extension=True).rank()


def zero(matrix):
    return all(s.simplify(x) == 0 for x in matrix)


def test_actual_native_lifts_reconstructed_from_full_native_columns(data):
    # The independent oracle forms the78 Grams directly from complete
    # 10626-row source columns rather than reusing response_matrices().
    _, columns = m.response.native_port_columns(0, 1)
    native = m.response.geometric_coordinates()
    R, A = native['frame'], data['native_chart']
    assert zero(R*A-data['geometric_frame'])
    assert zero(A.T*native['metric']*A-s.eye(6))
    grams = tuple(s.Matrix((columns[p].T@columns[q]).toarray())/11664
                  for p, q in m.response.PORT_PAIRS)
    for coefficients, h in zip(data['native_source_coefficients'], data['native_hermitian'], strict=True):
        lift = sum((c*g for c, g in zip(coefficients, grams, strict=True)), s.zeros(36))
        assert zero(lift*data['geometric_frame']-data['geometric_frame']*h)
    assert [sum(x != 0 for x in row) for row in data['native_source_coefficients']] == [1, 11, 11, 11, 11, 18]


def test_same_object_intersection_union_and_complete_observable_dimensions(data, certificate):
    G, D = data['native_basis'], data['declared_basis']
    assert (rank(G), rank(D), rank(G+D)) == (17, 12, 18)
    assert rank(G+data['intersection_basis']) == 17
    assert rank(D+data['intersection_basis']) == 12
    assert rank(data['intersection_basis']) == 11
    assert rank(tuple(g-s.trace(g)*s.eye(6)/6 for g in G)) == 16
    assert rank(tuple(g-s.trace(g)*s.eye(6)/6 for g in D)) == 12
    assert certificate['algebras']['union_Lie_closure_dimension'] == 18
    assert not certificate['algebras']['declared_algebra_is_a_subalgebra_of_native']


def test_exact_change_of_span_preserves_each_original_entry(data):
    for actual, canonical, forward, backward in (
        ('native_basis', 'canonical_native_basis', 'native_from_canonical', 'canonical_from_native'),
        ('declared_basis', 'canonical_declared_basis', 'declared_from_canonical', 'canonical_from_declared')):
        for j, generator in enumerate(data[actual]):
            assert zero(generator-sum((data[forward][i, j]*g for i, g in enumerate(data[canonical])), s.zeros(6)))
        assert zero(data[forward]*data[backward]-s.eye(len(data[actual])))


def test_center_mismatch_is_not_just_a_dimension_mismatch(data):
    G, D = data['native_basis'], data['declared_basis']
    identity, relative = s.I*s.eye(6), s.I*s.diag(1, 1, 1, 0, 0, 0)
    assert (rank(G+(identity,)), rank(G+(relative,))) == (17, 18)
    assert (rank(D+(identity,)), rank(D+(relative,))) == (13, 12)
    assert all(s.simplify(s.trace(g[:3, :3])-s.trace(g[3:, 3:])) == 0 for g in G)
    assert s.trace(relative[:3, :3])-s.trace(relative[3:, 3:]) == 3*s.I


def test_conserved_projectors_give_only_four_invariant_state_subspaces(data, certificate):
    p, q = data['commutant_basis']
    assert p+q == s.eye(6) and p*q == s.zeros(6)
    assert all(zero(g*p-p*g) and zero(g*q-q*g) for g in data['native_hermitian'])
    # A reducing projector in the exact two-dimensional commutant is aP+bQ;
    # idempotence forces a,b in{0,1}. Its ranks are0,3,3,6, never a carrier
    # carrying a faithful twelve-dimensional restricted current action.
    assert [int(s.trace(a*p+b*q)) for a, b in ((0, 0), (1, 0), (0, 1), (1, 1))] == [0, 3, 3, 6]
    assert rank(tuple(g[:3, :3] for g in data['native_basis'])) == 9
    assert rank(tuple(g[3:, 3:] for g in data['native_basis'])) == 9
    assert certificate['invariant_restrictions']['complex_commutant_dimension'] == 2


def test_fixed_pointer_menu_retains_all_outcomes_and_uses_no_response_eigenvectors():
    for probe in m.fixed_probes():
        rho, effects = probe['density'], probe['effects']
        assert s.trace(rho) == 1 and rho*rho == rho and rho == s.conjugate(rho)
        assert sum(effects, s.zeros(6)) == s.eye(6)
        assert [s.trace(e) for e in effects] == [1, 1, 4]
        assert all(e == e.conjugate().T and e*e == e for e in effects)
        assert all(effects[i]*effects[j] == s.zeros(6) for i in range(3) for j in range(i))


def test_five_extra_phases_have_independent_fixed_contrasts_already_before_closure(data, certificate):
    jacobian = m.response_jacobian(tuple(s.I*h for h in data['extra_hermitian_basis']))[:5, :]
    assert jacobian == s.Matrix([[-2, 1, 0, 0, 0], [1, -2, 0, 0, 0],
                                [0, 0, 2, 0, 0], [0, 0, 0, 2, 0], [0, 0, 0, 0, 2]])
    assert jacobian.det() == 24
    assert m.response_jacobian(tuple(s.I*h for h in data['native_hermitian']))[:5, :].rank() == 5
    assert m.response_jacobian(data['declared_basis'])[:5, :] == s.zeros(5, 12)
    assert certificate['fixed_records']['extra_phase_rank_already_in_six_source_Hermitian_responses'] == 5


def test_finite_extra_phase_has_exact_records_unavailable_to_real_triplet_rotations(certificate):
    witness = s.Matrix(certificate['fixed_records']['finite_extra_phase_unitary'])
    probe = m.fixed_probes()[0]
    law = [s.simplify(s.trace(e*witness*probe['density']*witness.conjugate().T)) for e in probe['effects']]
    assert law == [0, 1, 0]
    # Independent rational proper rotations keep every input/output real;
    # the imaginary antisymmetric Y contrast is consequently identically0.
    R = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5), 0],
                  [s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    rotation = s.diag(s.eye(3), R)
    for row in m.fixed_probes()[:5]:
        assert s.trace(row['contrast']*rotation*row['density']*rotation.T) == 0


def test_complete_central_tangent_functional_and_single_marker_false_positive(data):
    derivatives = m.response_jacobian(data['native_basis'])
    assert sum((derivatives[i, :] for i in range(5, 8)), s.zeros(1, 17)) == s.zeros(1, 17)
    central = m.response_jacobian((data['relative_center'],))
    assert central[5:, :] == -s.ones(3, 1)
    # A single cross-triplet marker is insufficient: native su3 can mimic
    # the center on that prepared state. The other coordinates reject it.
    mimic = s.I*s.diag(1, -1, 0, 0, 0, 0)
    assert rank(data['native_basis']+(mimic,)) == 17
    mimic_values = m.response_jacobian((mimic,))
    assert mimic_values[5] == central[5]
    assert mimic_values[5:, :] != central[5:, :]
    assert sum(mimic_values[5:, 0]) == 0


def test_finite_determinant_ratio_is_global_phase_invariant_and_has_discrete_exceptions(certificate):
    central = s.Matrix(certificate['finite_group_invariant']['finite_central_witness'])
    assert central[:3, :3].det()/central[3:, 3:].det() == -s.I
    for phase in (1, -1, s.I, (3+4*s.I)/5):
        U = phase*central
        assert s.simplify(U[:3, :3].det()/U[3:, 3:].det()) == -s.I
    omega = (-1+s.I*s.sqrt(3))/2
    allowed = s.diag(omega, omega, omega, 1, 1, 1)
    assert s.simplify(allowed[:3, :3].det()/allowed[3:, 3:].det()) == 1
    assert certificate['finite_group_invariant']['uniform_declared_phase_native_iff'] == 't belongs to (2 pi/3) Z'


def test_exact_derivatives_agree_with_independent_finite_complex_born_experiment(data):
    generator = np.asarray(data['native_basis'][7], complex)
    derivatives = np.asarray(m.response_jacobian((data['native_basis'][7],)), float)[:, 0]
    t = 2.**-15
    U, V = expm(t*generator), expm(-t*generator)
    for row, derivative in zip(m.fixed_probes(), derivatives, strict=True):
        rho, effect = np.asarray(row['density'], complex), np.asarray(row['contrast'], complex)
        actual = np.trace(effect@(U@rho@U.conj().T-V@rho@V.conj().T)).real/(2*t)
        assert abs(actual-derivative) < 1e-7


@pytest.mark.parametrize('bad', [[], (), (s.eye(6),), (s.zeros(5),), (s.I*s.eye(6)*s.Float(.5),),
                               (s.I*s.eye(6)*s.Symbol('x', real=True),), (s.diag(s.oo, 0, 0, 0, 0, 0),)])
def test_response_input_refuses_wrong_or_inexact_domains(bad):
    with pytest.raises(ValueError):
        m.response_jacobian(bad)


def test_owned_outputs_and_scope_cannot_be_mutated(certificate, data):
    certificate = m.exact_certificate()
    certificate['algebras']['native_real_dimension'] = 12
    assert m.exact_certificate()['algebras']['native_real_dimension'] == 17
    data = m.exact_data()
    data['native_basis'] = ()
    assert len(m.exact_data()['native_basis']) == 17
    with pytest.raises(TypeError):
        m.exact_data()['native_basis'][0][0, 0] = 12
    report = m.exact_certificate()
    assert json.loads(json.dumps(report)) == report
    assert not report['scope']['signed_native_port_Gram_schedules_source_executed_here']
    assert not report['scope']['fixed_pointer_source_words_executed_here']
    assert not report['scope']['native_C_E_repairs_select_either_current_grammar']


def test_subclass_cannot_hide_wrong_entries_with_overridden_iteration():
    class HiddenFloat(s.MutableDenseMatrix):
        def __iter__(self):
            return iter(s.zeros(6))
    wrong = HiddenFloat(s.I*s.Float(.5)*s.eye(6))
    with pytest.raises(ValueError, match='concrete trusted'):
        m.response_jacobian((wrong,))


def test_fixed_probe_containers_are_detached():
    probes = m.fixed_probes()
    probes[0]['density'] = s.zeros(6)
    probes[1]['effects'] = ()
    assert s.trace(m.fixed_probes()[0]['density']) == 1
    assert len(m.fixed_probes()[1]['effects']) == 3


def test_wrong_hermitian_response_with_unchanged_native_lift_is_rejected(data, monkeypatch):
    original = m.response.preserving_responses(0, 1)
    changed = dict(original)
    changed['basis'] = (2*original['basis'][0], *original['basis'][1:])
    monkeypatch.setattr(m.response, 'preserving_responses', lambda *args: changed)
    m._exact_data.cache_clear()
    try:
        with pytest.raises(ArithmeticError, match='complete-output source lift'):
            m.exact_data()
    finally:
        m._exact_data.cache_clear()
