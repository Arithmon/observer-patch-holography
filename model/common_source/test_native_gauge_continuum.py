from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import native_gauge_continuum as g
from . import native_gauge_matter_transport as transport
from . import native_staggered_field as field


@pytest.fixture(scope='module')
def arrays():
    return g.execute()


def test_actual_current_lift_preserves_selected_and_full_complement():
    a = g.architecture(); e = a['encoder']; q = s.eye(32)-e*e.T
    for index, A in zip(g.CURRENT_INDICES, a['currents'], strict=True):
        assert A+A.H == s.zeros(32)
        assert g.matter.zero(q*A*e)
        assert g.matter.zero(e.T*A*q)
        assert g.matter.zero(e.T*A*e-transport.exact_currents()[index])
    assert any(not g.matter.zero(x) for x in a['commutators'])


def test_complete_source_links_and_report(arrays):
    report = g.report_from_arrays(arrays)
    assert report['selected_Dirac_fields_with_family_factor'] == 90
    assert report['retained_auxiliary_complement_dimension'] == 17
    assert report['physical_chiral_SM_or_dynamical_gauge_continuum_claimed'] is False
    assert report['fixed_chirality_site_gauge_leakage_Frobenius_diagnostic'] > 0.1
    assert report['omitted_field_strength_operator_norm_diagnostic'] > 0.1
    for row in report['source_link_refinements']:
        assert row['symbol_discretization_error_diagnostic'] <= (
            row['sine_bound_diagnostic']+row['fixed_source_link_error_upper_conditional'])


def test_noncommuting_uniform_sine_bound_and_finite_time_duhamel(arrays):
    A = arrays['currents']; p = (0.4, -0.1, 0.3); c = float(field.native_parameters()['coordinate_speed'])
    hd = g.symbol(p, F(1), A, continuum=True)
    previous = float('inf')
    for a in g.SPACINGS:
        h = g.symbol(p, a, A)
        bound = c*float(a)**2*sum(np.linalg.norm(p[j]*np.eye(32)-1j*A[j], 2)**3 for j in range(3))/6
        error = np.linalg.norm(h-hd, 2)
        assert error <= bound+1e-13
        assert error < previous/3.9
        assert np.linalg.norm(expm(-1j*h/3)-expm(-1j*hd/3), 2) <= error/3+1e-13
        previous = error


def test_full_real_space_action_intertwines_every_corner_block(arrays):
    # A factor proof covers all32 internal columns; this direct complete
    # 64x32x(8x32) application independently checks the link orientation.
    d = 32; W = arrays['site_fourier']; a = g.SPACINGS[0]
    for block, b in enumerate(product((0, 1), repeat=3)):
        wave = np.kron(W[:, 8*block:8*block+8], np.eye(d)).reshape(64,d,8*d)
        p = -np.pi*np.asarray(b)/(2*float(a))
        h = g.source_symbol(p, a, arrays['links'][0])
        actual = g.apply_lattice(wave, arrays['links'][0], spacing=a)
        expected = np.einsum('xij,jk->xik', wave, h)
        assert np.max(np.abs(actual-expected)) < 2e-12


def test_local_gauge_covariance_on_all_real_space_link_columns(arrays):
    links = arrays['links'][0]; gauges = np.repeat(np.eye(32, dtype=complex)[None], 64, axis=0)
    gauges[0] = arrays['site_gauge']; gauges[21] = arrays['site_gauge'].conj().T
    changed = g.transformed_links(links, gauges)
    # This checks every nonzero Hamiltonian block/all internal columns,
    # not a single prepared state or a constant-background-only formula.
    for j, shift in enumerate(g.site_factors()[0]):
        for x, y in zip(*np.nonzero(shift), strict=True):
            assert np.linalg.norm(changed[x,j]-gauges[x]@links[j]@gauges[y].conj().T, 2) < 1e-14
            assert np.linalg.norm(changed[x,j]@gauges[y]-gauges[x]@links[j], 2) < 1e-11
    rng = np.random.default_rng(831)
    state = rng.standard_normal((64,32,2))+1j*rng.standard_normal((64,32,2))
    transformed = np.einsum('xij,xja->xia', gauges, state)
    left = g.apply_lattice(transformed, changed, spacing=g.SPACINGS[0])
    right = np.einsum('xij,xja->xia', gauges, g.apply_lattice(state, links, spacing=g.SPACINGS[0]))
    assert np.linalg.norm(left-right) < 2e-10
    assert np.linalg.norm(changed[0]-changed[10]) > 0.1


def test_both_tastes_and_chiralities_retained_for_constant_background(arrays):
    I, Z, Y, X = np.eye(2), np.diag([1,-1]), np.array([[0,-1j],[1j,0]]), np.array([[0,1],[1,0]])
    tastes = [np.kron(np.kron(I,Z),Y), np.kron(np.kron(Z,X),Y), np.kron(np.kron(Z,Y),I)]
    h = g.symbol((0.3,0.2,-0.1), F(1,4), arrays['currents'])
    for t in tastes:
        T = np.kron(t, np.eye(32))
        assert np.linalg.norm(T@h-h@T) < 1e-12
    alpha = [np.asarray(x, complex) for x in field.clifford_matrices()]
    chi = 1j*alpha[0]@alpha[1]@alpha[2]
    assert np.count_nonzero(np.linalg.eigvalsh(chi) > 0) == 4
    h0 = g.symbol((0.3,0.2,-0.1), F(1,4), arrays['currents'], mass=0)
    Chi = np.kron(chi, np.eye(32))
    assert np.linalg.norm(Chi@h0-h0@Chi) < 1e-12


def test_same_target_free_family_mass_commutes_with_actual_gauge_action(arrays):
    Y = arrays['family_mass']; U = arrays['links'][0,0]
    mass = np.kron(np.eye(32), Y); gauge = np.kron(U, np.eye(3))
    assert np.linalg.norm(gauge@mass-mass@gauge) < 1e-13
    assert np.linalg.eigvalsh(Y)[0] > 0
    assert np.linalg.norm(Y, 2) <= 1
    p = (0.2,-0.3,0.1); a = F(1,4)
    delta = g.family_symbol(p,a,arrays['currents'],Y)-g.family_symbol(p,a,arrays['currents'],Y,continuum=True)
    old = g.symbol(p,a,arrays['currents'])-g.symbol(p,a,arrays['currents'],continuum=True)
    assert np.linalg.norm(delta-np.kron(old,np.eye(3))) < 1e-12
    with pytest.raises(ValueError): g.family_symbol(p,a,arrays['currents'],Y[:2,:2])


def test_fixed_chirality_does_not_survive_one_actual_site_rotation(arrays):
    alpha = [np.asarray(x, complex) for x in field.clifford_matrices()]
    chi = 1j*alpha[0]@alpha[1]@alpha[2]; W = arrays['site_fourier']
    P = W@np.kron(np.eye(8),(np.eye(8)+chi)/2)@W.conj().T
    assert np.linalg.norm(P@P-P) < 1e-13
    # Exact factorization of (I-PtensorI) Gsite (PtensorI) removes the
    # identity term. Every unselected row remains in the leakage operator.
    E = np.zeros((64,64)); E[0,0] = 1
    left = (np.eye(64)-P)@E@P
    assert np.linalg.norm(left) > 0.1
    assert np.linalg.norm(arrays['site_gauge']-np.eye(32)) > 0.1


def test_refinement_contract_requires_improving_source_accuracy():
    kwargs = dict(duration=F(1), momentum_cutoff=F(1), particles=2)
    previous = None
    for q in (4,8,16):
        a = F(1,q); n = q**4; count = 7*(4*q)**3*n
        result = g.refinement_contract(spacing=a, splitting_steps=n,
            source_macro_count=count, source_macro_error=a*a/(2*count), **kwargs)
        assert result['source_precision_scaling_implemented'] is False
        errors = tuple(F(result[k]) for k in ('one_body_discretization_error_upper',
            'one_body_product_formula_error_upper', 'complete_Fock_source_error_upper_conditional'))
        if previous is not None:
            assert all(x < y/3.5 for x,y in zip(errors, previous, strict=True))
        previous = errors
    with pytest.raises(ValueError, match='accumulated source'):
        g.refinement_contract(spacing=F(1,16), splitting_steps=16**4,
            source_macro_count=7*64**3*16**4, source_macro_error=F(1,10**10), **kwargs)


def test_nonunitary_source_norm_budget_is_not_clipped_at_two():
    result = g.refinement_contract(spacing=F(1,4), duration=F(1), momentum_cutoff=F(1),
        splitting_steps=1, source_macro_count=1, source_macro_error=F(3,4))
    assert F(result['complete_Fock_source_error_upper_conditional']) == 3
    assert F(result['fixed_particle_propagator_error_upper_conditional']) == 5


@pytest.mark.parametrize('change', [dict(spacing=0), dict(spacing=0.25),
    dict(splitting_steps=True), dict(source_macro_count=-1), dict(momentum_cutoff=100)])
def test_refinement_domains_fail_closed(change):
    kwargs = dict(spacing=F(1,4), duration=F(1), momentum_cutoff=F(1),
        splitting_steps=64, source_macro_count=100, source_macro_error=F(1,10**12))
    kwargs.update(change)
    with pytest.raises(ValueError): g.refinement_contract(**kwargs)


def test_retained_factor_menu_shape_current_and_nan_are_checked(arrays):
    for bad in ({k:v for k,v in arrays.items() if k != 'site_gauge'}, dict(arrays, extra=np.zeros(1))):
        with pytest.raises(ValueError): g.validate_arrays(bad)
    for key in arrays:
        bad = {k:v.copy() for k,v in arrays.items()}; bad[key].flat[0] = np.nan
        with pytest.raises(ValueError): g.validate_arrays(bad)
    bad = {k:v.copy() for k,v in arrays.items()}; bad['currents'][0,0,0] += 0.01j
    with pytest.raises(ValueError): g.validate_arrays(bad)


def test_own_array_report_binding_rejects_claim_and_bound_tampering(arrays):
    import copy
    import json
    report = g.report_from_arrays(arrays)
    assert g.verify_arrays(report, arrays)
    assert g.verify_arrays(json.loads(json.dumps(report, allow_nan=False)), arrays)
    bad = copy.deepcopy(report); bad['physical_chiral_SM_or_dynamical_gauge_continuum_claimed'] = True
    with pytest.raises(ValueError): g.verify_arrays(bad, arrays)
    bad = copy.deepcopy(report)
    bad['source_link_refinements'][0]['fixed_source_link_error_upper_conditional'] *= 1.00001
    with pytest.raises(ValueError): g.verify_arrays(bad, arrays)
    bad = copy.deepcopy(report); del bad['refinement_contract_example']['source_precision_scaling_implemented']
    with pytest.raises(ValueError): g.verify_arrays(bad, arrays)


def test_source_symbol_keeps_nonunitary_link_leakage_in_hermitian_action(arrays):
    links = arrays['links'][0].copy(); links[0,0,0] += 0.01
    h = g.source_symbol((0.2,-0.1,0.3), F(1,4), links)
    assert np.linalg.norm(h-h.conj().T) < 1e-13
    assert np.linalg.norm(h-g.source_symbol((0.2,-0.1,0.3), F(1,4), arrays['links'][0])) > 0.001
