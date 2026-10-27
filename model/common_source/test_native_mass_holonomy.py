"""Native geometric phase lattice and explicitly conditional mass map."""
from fractions import Fraction as F

import pytest
import sympy as s

from . import native_mass_holonomy as m


@pytest.fixture(scope='module')
def data():
    return m.geometry()


def test_all_actual_faces_have_the_exact_spin_coherent_phase(data):
    assert len(data['faces']) == len(data['products']) == 20
    assert m._scalar(data['phase']**5-s.I) == 0
    assert m._scalar(data['phase']**20-1) == 0
    for a, b, c in data['faces']:
        reverse = s.trace(data['projectors'][c]*data['projectors'][b]*data['projectors'][a])
        forward = s.trace(data['projectors'][a]*data['projectors'][b]*data['projectors'][c])
        assert m._scalar(reverse-s.conjugate(forward)) == 0


def test_integer_cycle_certificate_is_unimodular_and_covers_large_windings(data):
    assert abs(data['reduced_determinant']) == 1
    assert data['boundary1']*data['boundary2'] == s.zeros(12, 20)
    coefficients = s.Matrix([(-1)**i*(i+1)**3 for i in range(20)])
    chain = data['boundary2']*coefficients
    decoded = m.integral_face_decomposition([int(x) for x in chain])
    assert data['boundary2']*s.Matrix(decoded) == chain
    assert (sum(decoded)-sum(coefficients)) % 20 == 0
    for index in range(20):
        decoded = m.integral_face_decomposition([int(x) for x in data['boundary2'][:, index]])
        assert sum(decoded) % 20 == 1


def test_spin_one_symmetric_square_doubles_the_same_source_phase(data):
    basis = s.Matrix([[1, 0, 0], [0, 1/s.sqrt(2), 0], [0, 1/s.sqrt(2), 0], [0, 0, 1]])
    assert basis.T*basis == s.eye(3)
    a, b, c = data['faces'][0]
    ps = [basis.T*s.kronecker_product(data['projectors'][i], data['projectors'][i])*basis
          for i in (a, b, c)]
    for p in ps: assert m._zero(p*p-p) and m._scalar(s.trace(p)-1) == 0
    first = s.trace(data['projectors'][a]*data['projectors'][b]*data['projectors'][c])
    assert m._scalar(s.trace(ps[0]*ps[1]*ps[2])-first**2) == 0


def test_cycle_decoder_rejects_nonclosed_and_noninteger_inputs():
    with pytest.raises(ValueError): m.integral_face_decomposition([1]+[0]*29)
    with pytest.raises(ValueError): m.integral_face_decomposition([0.]*30)
    with pytest.raises(ValueError): m.integral_face_decomposition([False]*30)
    with pytest.raises(ValueError): m.integral_face_decomposition([0]*29)


def test_native_family_action_really_is_regular_c3(data):
    R = data['C3_family_action']
    assert m._zero(R**3-s.eye(3)) and s.trace(R) == 0
    assert m._zero(R.T*R-s.eye(3))
    # All three C3 characters occur once, rather than an assumed cyclic label.
    assert R.charpoly().as_expr() == R.charpoly().gen**3-1


def test_actual_family_circulant_has_the_claimed_conditional_mass_invariants(data):
    x, y = s.symbols('x y', real=True)
    R = data['C3_family_action']
    C = s.eye(3)+(x+s.I*y)*R+(x-s.I*y)*R**2
    assert m._zero(C.conjugate().T-C)
    assert m._scalar(s.trace(C)-3) == 0
    e2 = (s.trace(C)**2-s.trace(C*C))/2
    assert m._scalar(e2-(3-3*(x*x+y*y))) == 0
    assert m._scalar(C.det()-(1-3*(x*x+y*y)+2*(x**3-3*x*y*y))) == 0


def test_conditional_positive_root_certificates_and_phase_classes():
    rows = [m.conditional_spectrum(k) for k in (0, 1, 2)]
    for row in rows:
        intervals = [list(map(F, pair)) for pair in row['positive_root_intervals_exact']]
        assert all(0 < lo <= hi for lo, hi in intervals)
        assert sum(pair[0] for pair in intervals) <= 3 <= sum(pair[1] for pair in intervals)
    assert rows[0]['conditional_mass_ratios_after_one_exact'][0] == ['1', '1']
    for bad in (-1, 3, True):
        with pytest.raises(ValueError): m.conditional_spectrum(bad)


def test_both_fixed_rules_and_the_full_phase_lattice_fail_the_heldout_windows():
    report = m.report()
    contract = report['conditional_mass_contract']
    assert contract['supplied_balance'] == '|b/a|=1/sqrt(2)'
    assert not contract['actual_native_GNS_orientation_attachment_executed_here']
    assert not contract['source_selected_physical_square_root_mass_map']
    assert contract['positive_spin_half_lattice_residues_mod20'] == [0, 1, 2, 18, 19]
    comparison = report['separate_heldout_comparison']
    assert comparison['complete_balanced_geometric_loop_class_excluded']
    assert all(row['window_intersections'] == [False, False] for row in comparison['conditional_phase_classes'])
    assert not comparison['continuously_chosen_or_refined_connection_phases_excluded']


def test_source_phase_and_conditional_spectrum_do_not_read_mass_data(monkeypatch):
    def unavailable(*args, **kwargs): raise AssertionError('empirical data entered the phase producer')
    monkeypatch.setattr(m.strict_json, 'load', unavailable)
    m.geometry()
    for k in (0, 1, 2): m.conditional_spectrum(k)


def test_rational_square_root_enclosures_are_outward():
    for bounds in ((F(2), F(2)), (F(5), F(5)), (F(7, 13), F(9, 11))):
        lo, hi = m._sqrt_interval(bounds)
        assert lo*lo <= bounds[0] <= bounds[1] < hi*hi
