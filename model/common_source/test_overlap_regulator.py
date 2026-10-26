"""Independent Wilson, chiral, gauge and complete native-kernel controls."""
from itertools import product
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import block_diag, expm

from . import overlap_regulator as r


def test_actual_exterior_module_preserves_nonabelian_group_product_and_charges():
    a = np.zeros((5, 5), complex)
    a[0, 1] = a[1, 0] = .21
    a[3, 4] = -.13j
    a[4, 3] = .13j
    b = np.diag([.17, -.11, -.06, .07, -.07])
    u, v = expm(1j*a), expm(1j*b)
    rho_u, rho_v = r.exterior_link(u), r.exterior_link(v)
    assert np.linalg.norm(r.exterior_link(u@v)-rho_u@rho_v, 'fro') < 1e-12
    assert np.linalg.norm(rho_u.conj().T@rho_u-np.eye(15), 'fro') < 1e-12
    assert [len(s) for s in r.matter_basis()].count(2) == 10
    assert [len(s) for s in r.matter_basis()].count(4) == 5


@pytest.mark.parametrize('bad', [np.zeros((5, 5)), np.eye(5)*1j, np.ma.array(np.eye(5)), np.full((5, 5), np.nan)])
def test_nonunitary_wrong_determinant_and_masked_links_are_rejected(bad):
    with pytest.raises((ValueError, ArithmeticError)):
        r.exterior_link(bad)


def test_all_free_lattice_input_columns_match_independent_fourier_symbol():
    # Fourier modes are enumerated independently of the Wilson assembly.
    links = np.ones((16, 4, 1, 1), complex)
    wilson, g5 = r.wilson_operator(links)
    d, hat, _ = r.overlap(wilson, g5)
    sites = tuple(product((0, 1), repeat=4))
    fourier = np.array([[(-1)**sum(a*b for a, b in zip(x, p, strict=True))/4
                         for p in sites] for x in sites])
    transformed = np.kron(fourier.T, np.eye(4))@d@np.kron(fourier, np.eye(4))
    expected = block_diag(*(np.zeros((4, 4)) if not any(p) else 2*np.eye(4) for p in sites))
    assert np.linalg.norm(transformed-expected, 'fro') < 1e-12
    assert np.count_nonzero(np.linalg.eigvalsh(d.conj().T@d) < 1e-10) == 4
    # A left Weyl field has two spinor components, not a second Dirac species.
    assert round(np.trace((np.eye(64)-hat)/2).real) == 32
    assert np.linalg.norm(g5@d+d@g5-d@g5@d, 'fro') < 1e-12


def test_nonconstant_gauge_transform_is_checked_on_actual_wilson_links():
    links = np.ones((16, 4, 1, 1), complex)
    gauges = np.exp(1j*np.arange(16)*.037).reshape(16, 1, 1)
    changed, whole = r.gauge_transform(links, gauges)
    w, g5 = r.wilson_operator(links)
    new_w, _ = r.wilson_operator(changed)
    d, hat, _ = r.overlap(w, g5)
    new_d, new_hat, _ = r.overlap(new_w, g5)
    assert np.linalg.norm(new_w-whole@w@whole.conj().T, 'fro') < 1e-12
    assert np.linalg.norm(new_d-whole@d@whole.conj().T, 'fro') < 1e-12
    assert np.linalg.norm(new_hat-whole@hat@whole.conj().T, 'fro') < 1e-12
    wrong = changed.copy()
    wrong[0, 1] = changed[0, 1].conj()
    wrong_w, _ = r.wilson_operator(wrong)
    assert np.linalg.norm(wrong_w-whole@w@whole.conj().T, 'fro') > .01
    with pytest.raises(ArithmeticError, match='gauge unitarity'):
        r.gauge_transform(links, np.ones((16, 1, 1))*2)
    real_links = np.ones((16, 4, 1, 1), float)
    complex_result, _ = r.gauge_transform(real_links, gauges)
    assert complex_result.dtype.kind == 'c'
    assert abs(complex_result[0, 0, 0, 0].imag) > .01
    with pytest.raises(ValueError, match='gauge matrix'):
        r.gauge_transform(real_links, np.ma.array(gauges))


def test_modified_projector_is_needed_for_the_gauged_weyl_map():
    links = np.ones((16, 4, 1, 1), complex)
    links[0, 0, 0, 0] = np.exp(.07j)
    w, g5 = r.wilson_operator(links)
    d, hat, _ = r.overlap(w, g5)
    left = (np.eye(64)-hat)/2
    right = (np.eye(64)+g5)/2
    assert np.linalg.norm(d@left-right@d, 'fro') < 1e-12
    assert np.linalg.norm(d@(np.eye(64)-g5)/2-right@d, 'fro') > 1


@pytest.mark.parametrize('p', [(1/3, 0., 0., 0.), (.17, -.29, .11, -.13)])
def test_free_continuum_error_is_controlled_under_one_fixed_momentum(p):
    gamma, _ = r.gamma_matrices()
    target = 1j*sum((x*g for x, g in zip(p, gamma, strict=True)), np.zeros((4, 4), complex))
    errors = []
    for a in (.125, .0625, .03125):
        d, _ = r.free_symbol(p, a)
        bound = r.free_continuum_bound(p, a)
        errors.append(np.linalg.norm(d-target, 2))
        assert errors[-1] <= bound
    assert errors[2] < errors[1] < errors[0]
    with pytest.raises(ValueError, match='comparison band'):
        r.free_continuum_bound(p, 10.)


def test_complete_native_query_keeps_failure_phase_and_arbitrary_reference_columns():
    momentum = (1/3, 1/5, 1/7, 1/11)
    report, actual = r.native_free_query(momentum)
    d, v = r.free_symbol(momentum, .125)
    target = np.vstack((.125*d/2, (np.eye(4)-v)/2))
    assert np.linalg.norm(actual-target, 2) < 3e-6
    assert actual.shape == (8, 4)
    assert np.linalg.norm(actual.conj().T@actual-np.eye(4), 'fro') < 1e-7
    assert report['source_native_pulses'] > 10000
    # One untouched maximally entangled reference retains every input column.
    entangled = actual/2
    changed = actual.copy()
    changed[4:] *= 1j
    assert np.allclose(abs(changed), abs(actual))
    assert np.linalg.norm(changed/2-entangled, 'fro') > .1


def test_zero_success_branch_remains_a_valid_complete_native_query():
    report, actual = r.native_free_query((0., 0., 0., 0.))
    assert np.linalg.norm(actual[:4], 'fro') < 1e-12
    assert np.linalg.norm(actual[4:]-np.eye(4), 'fro') < 1e-12
    assert report['ideal_trotter_operator_bound_exact'] == '0'


def test_wrong_current_primitive_wire_is_not_hidden_by_probability_comparison(monkeypatch):
    original = r.compiler.execute
    def changed(word, qubits):
        result = original(word, qubits)
        # Preserve entry magnitudes but corrupt an actual source phase.
        return 1j*result
    monkeypatch.setattr(r.compiler, 'execute', changed)
    with pytest.raises(ArithmeticError):
        r.native_free_query((1/3, 1/5, 1/7, 1/11))


def test_wilson_sign_gap_and_original_numeric_domains_are_enforced():
    with pytest.raises(ValueError, match='spectral gap'):
        r.overlap(np.eye(4), np.eye(4))
    with pytest.raises(ValueError):
        r.free_symbol([True]*4)
    with pytest.raises(ValueError):
        r.free_symbol([0]*4, True)


def test_local_polynomial_bound_uses_the_actual_gapped_wilson_operator():
    links = np.ones((16, 4, 1, 1), complex)
    links[0, 0, 0, 0] = np.exp(.013j)
    w, g5 = r.wilson_operator(links)
    exact, _, _ = r.overlap(w, g5)
    error = []
    for degree in (64, 128, 192):
        approximation, certificate = r.local_overlap_polynomial(w, g5, degree, links=links)
        bound = F(certificate['uniform_operator_tail_exact'])
        assert bound == 16*F(15, 17)**(degree+1)
        assert certificate['maximum_range_in_lattice_links'] == 2*degree+1
        defect = g5@approximation+approximation@g5-approximation@g5@approximation
        error.append(np.linalg.norm(approximation-exact, 2))
        assert error[-1] <= float(bound)+1e-11
        assert np.linalg.norm(defect, 2) <= float(F(certificate['Ginsparg_Wilson_defect_bound_exact']))+1e-11
    assert error[-1] < error[0]
    with pytest.raises(ValueError, match='not built from these local links'):
        r.local_overlap_polynomial(np.eye(4), np.eye(4), links=links)
    with pytest.raises(ValueError, match='on-site chirality'):
        r.local_overlap_polynomial(w, .5*g5, links=links)
    with pytest.raises(ValueError):
        r.local_overlap_polynomial(w, g5, True, links=links)
    with pytest.raises(TypeError):
        r.local_overlap_polynomial(w, g5)
    # A dense, gapped sign kernel alone does not imply finite lattice range.
    all_sites = np.ones(16)/4
    dense_h = np.kron(np.eye(16)-2*np.outer(all_sites, all_sites), np.eye(4))
    forged_w = np.eye(64)+g5@dense_h
    with pytest.raises(ValueError, match='nearest-neighbor Wilson assembly'):
        r.local_overlap_polynomial(forged_w, g5, 0, links=links)
    tiny_nonlocal = w.copy()
    tiny_nonlocal[0:4, 60:64] += 1e-12*np.eye(4)
    tiny_nonlocal[60:64, 0:4] += 1e-12*np.eye(4)
    with pytest.raises(ValueError, match='nearest-neighbor Wilson assembly'):
        r.local_overlap_polynomial(tiny_nonlocal, g5, 0, links=links)


def test_receipt_comparison_keeps_phases_zeros_exact_contracts_and_types():
    report, _ = r.native_free_query((1/3, 1/5, 1/7, 1/11))
    r.compare(deepcopy(report), report)
    bad = deepcopy(report)
    column = bad['complete_success_and_failure_reference_columns'][0][0]
    column[0], column[1] = -column[1], column[0]
    with pytest.raises(ValueError, match='complex native output'):
        r.compare(bad, report)
    bad = deepcopy(report)
    bad['source_numerical_operator_admission'] *= 1.000001
    with pytest.raises(ValueError, match='contract changed'):
        r.compare(bad, report)
    bad = deepcopy(report)
    bad['both_quantum_outcomes_and_all_input_columns_checked'] = 1
    with pytest.raises(ValueError, match='type mismatch'):
        r.compare(bad, report)


def test_background_satisfies_the_declared_representation_uniform_domain():
    assert r.GENERATOR_BOUND == F(7, 5000)
    assert r.PLAQUETTE_BOUND == F(14, 625) < F(1, 40)
    links = r.sample_links()
    for matrix in links.reshape(-1, 15, 15):
        assert np.linalg.norm(matrix-np.eye(15), 2) <= float(4*r.GENERATOR_BOUND)
    assert 0 < r.plaquette_defect(links) < 1/40
