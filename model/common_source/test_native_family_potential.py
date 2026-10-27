"""Exact native potential invariants, stationary completeness and scope."""
import pytest
import sympy as s

from . import native_family_potential as p


@pytest.fixture(scope='module')
def data():
    return p.architecture()


def test_complete_native_character_average_and_scalar_moments(data):
    assert sum(data['trace_counts'].values()) == 60
    assert data['dimensions'] == (1, 0, 1, 0, 1, 0, 2)
    x = data['variables']; norm2 = (x.T*x)[0]
    assert p._scalar(data['moments'][2]-4*norm2) == 0
    assert p._scalar(data['moments'][4]-s.Rational(12, 5)*norm2**2) == 0


def test_first_anisotropy_and_all_stationary_strata_are_distinct(data):
    assert len(data['interior_points']) == 4 and len(data['boundary_roots']) == 2
    assert {name: row['size'] for name, row in data['classes'].items()} == {'vertex': 12, 'face': 20, 'edge': 30}
    assert data['classes']['vertex']['value'] == s.Rational(52, 25)
    assert data['classes']['face']['value'] == s.Rational(68, 45)
    assert data['classes']['edge']['value'] == s.Rational(8, 5)
    assert sum(row['size'] for row in data['classes'].values()) == 62


def test_actual_stationary_curvatures_have_the_correct_morse_types(data):
    for name in ('vertex', 'face'):
        a, b = data['classes'][name]['tangent_eigenvalues']
        assert a == b and (a.is_negative if name == 'vertex' else a.is_positive)
    a, b = data['classes']['edge']['tangent_eigenvalues']
    assert a.is_negative and b.is_positive


def test_potential_coordinates_are_the_same_selected_native_triplet(data):
    F = data['native_family_isometry']; selected = p.selector.geometry()
    assert F.shape == (36, 3)
    assert p._zero(F.T*F-s.eye(3)) and p._zero(F*F.T-selected['P'])


def test_native_positive_gram_derivative_has_exact_selected_doublets(data):
    x = data['variables']
    assert p._zero(30*data['source_Gram']-s.hessian(data['moments'][6], x))
    assert data['classes']['face']['source_Gram_eigenvalues'] == (s.Rational(4, 9), s.Rational(4, 9), s.Rational(68, 45))
    assert data['classes']['vertex']['source_Gram_eigenvalues'] == (s.Rational(4, 25), s.Rational(4, 25), s.Rational(52, 25))
    assert data['nonstationary_Gram_discriminant'].is_positive


def test_general_complex_yukawa_counterexample_prevents_false_degeneracy_claim():
    # Any native unit axis is orthogonally conjugate to this one. This matrix
    # commutes with every rotation about that axis, including C3 and C5.
    u = s.Matrix([0, 0, 1]); hat = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    Y = 3*s.eye(3)+2*u*u.T+s.I*hat
    assert Y.H == Y and Y*hat == hat*Y
    assert Y.eigenvals() == {s.Integer(2): 1, s.Integer(4): 1, s.Integer(5): 1}


def test_report_retains_physical_and_representation_domains(data):
    r = p.report()
    assert r['empirical_inputs'] == []
    assert r['first_anisotropic_degree'] == 6
    assert not r['fixed_radius_selection']['radial_stability_or_nonzero_vacuum_selected']
    b = r['coupling_boundary']
    assert not b['tangent_scalar_curvature_is_fermion_Yukawa_spectrum']
    assert not b['C3_or_C5_alone_forces_general_complex_Yukawa_degeneracy']
