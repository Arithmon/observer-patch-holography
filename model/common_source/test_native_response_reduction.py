"""Exact controls for channel and quotient reductions of the native response."""
from fractions import Fraction
import json

import numpy as np
import pytest
import sympy as s

from . import native_response_reduction as r


def test_source_target_centers_and_quotients_are_not_conflated():
    report = r.canonical_certificate()
    rows = report['canonical_algebras']
    assert (rows['source_real_dimension'], rows['source_observable_dimension']) == (17, 16)
    assert (rows['target_real_dimension'], rows['target_observable_dimension']) == (12, 12)
    assert rows['intersection_real_dimension'] == 11
    assert not rows['target_relative_center_in_source']
    assert 12 not in report['exact_unitary_or_quotient']['surjective_real_Lie_quotient_dimensions']


@pytest.mark.parametrize('a,b', [(1, 1), (1, 0), (0, 1), (-1, 1),
                                (Fraction(2, 3), 0), (0, 0), (Fraction(-1, 2), Fraction(1, 10))])
def test_full_choi_spectrum_matches_exact_three_sector_formula(a, b):
    expected = {}
    for value, multiplicity in r.multipole_eigenvalues(a, b):
        expected[value] = expected.get(value, 0)+multiplicity
    assert r.multipole_choi(a, b).eigenvals() == expected


@pytest.mark.parametrize('a,b', [(1, 0), (0, 1), (-1, 1), (Fraction(2, 3)+Fraction(1, 1000), 0)])
def test_unphysical_projection_and_too_good_erasure_are_rejected(a, b):
    with pytest.raises(ValueError, match='not completely positive'):
        r.covariant_channel(a, b)


@pytest.mark.parametrize('value', [True, False, 0.5, np.float64(.5), np.int64(1), '1/2', None])
def test_exact_coefficients_reject_coercions(value):
    with pytest.raises(ValueError, match='original exact rational'):
        r.covariant_channel(value, 0)


def test_exact_fixing_rotations_reconstructs_all_matrix_units():
    js = r.spin_generators()
    assert sum((j*j for j in js), s.zeros(3)) == 2*s.eye(3)
    for i in range(3):
        for j in range(3):
            assert int(i == j)*s.eye(3)-js[i]*js[j] == r._unit(j, i)


def test_projector_hypothesis_is_needed_on_the_full_carrier():
    p = s.diag(0, 0, 0, 1, 1, 1)
    kraus = (p, *(r._unit(3, i, 6) for i in range(3)))
    assert sum((k.conjugate().T*k for k in kraus), s.zeros(6)) == s.eye(6)
    apply = lambda x: sum((k.conjugate().T*x*k for k in kraus), s.zeros(6))
    for j in r.spin_generators():
        x = s.diag(s.zeros(3), j)
        assert apply(x) == x
    assert apply(p) == s.eye(6)
    q = s.diag(0, 0, 0, 1, -1, 0)
    assert apply(q) != q


def test_physical_quintet_erasure_loses_one_third_of_every_rotation():
    channel = r.covariant_channel(Fraction(2, 3), 0)
    for j in r.spin_generators():
        assert channel*j.reshape(9, 1) == s.Rational(2, 3)*j.reshape(9, 1)
    for q in (s.diag(1, -1, 0), r._unit(0, 1)+r._unit(1, 0)):
        assert channel*q.reshape(9, 1) == s.zeros(9, 1)
    assert channel*s.eye(3).reshape(9, 1) == s.eye(3).reshape(9, 1)


def test_twirl_erases_both_nontrivial_sectors():
    channel = r.covariant_channel(0, 0)
    for x in (*r.spin_generators(), s.diag(1, -1, 0)):
        assert channel*x.reshape(9, 1) == s.zeros(9, 1)


def test_basis_snapshots_are_immutable():
    basis = r.canonical_bases()
    with pytest.raises(TypeError):
        basis['source'][0][0, 0] = 1
    basis['source'] = ()
    assert len(r.canonical_bases()['source']) == 17


def test_wrong_spin_normalization_cannot_support_the_channel_obstruction(monkeypatch):
    js = r.spin_generators()
    monkeypatch.setattr(r, 'spin_generators', lambda: (2*js[0], *js[1:]))
    with pytest.raises(ArithmeticError, match='Casimir'):
        r.canonical_certificate()


def test_false_choi_spectrum_is_not_an_unchecked_report_claim(monkeypatch):
    original = r.multipole_eigenvalues

    def shifted(a, b):
        rows = original(a, b)
        return ((rows[0][0]+1, rows[0][1]), *rows[1:])

    monkeypatch.setattr(r, 'multipole_eigenvalues', shifted)
    with pytest.raises(ArithmeticError, match='affine Choi'):
        r.canonical_certificate()


@pytest.fixture(scope='module')
def source_report():
    return r.report()


def test_actual_source_matrices_and_full_native_covariance_bind_the_theorem(source_report):
    assert json.loads(json.dumps(source_report, allow_nan=False)) == source_report
    actual = source_report['actual_source_binding']
    assert actual['both_exact_real_spans_independently_checked']
    assert actual['right_block_is_actual_native_three_prime']
    assert actual['actual_proper_native_rotations_checked'] == 60
    assert actual['right_EndC3_irreducible_dimensions'] == [1, 3, 5]
    assert actual['exact_character_inner_product_matrix'] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


@pytest.mark.parametrize('mutation', ['lost_generator', 'relative_center', 'wrong_triplet'])
def test_false_native_binding_is_rejected(source_report, monkeypatch, mutation):
    from . import native_response_comparison as comparison
    data = comparison.exact_data()
    if mutation == 'lost_generator':
        data['native_basis'] = data['native_basis'][:-1]
    elif mutation == 'relative_center':
        data['native_basis'] = (*data['native_basis'][:-1], r.canonical_bases()['target_center'])
    else:
        frame = data['geometric_frame']
        data['geometric_frame'] = s.ImmutableMatrix.hstack(frame[:, 3:], frame[:, :3])
    monkeypatch.setattr(comparison, 'exact_data', lambda: data)
    with pytest.raises(ArithmeticError, match='span|three-prime'):
        r.report()
