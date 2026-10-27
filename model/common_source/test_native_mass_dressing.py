"""Source-current covariance, full auxiliary resolvent and retained spectra."""
from copy import deepcopy

import numpy as np
import pytest
import sympy as s
from scipy import sparse
from scipy.sparse.linalg import spsolve

from . import native_mass_dressing as m


@pytest.fixture(scope='module')
def execution():
    return m.execute()


def test_current_map_and_adjoint_auxiliary_indices_are_complete(execution):
    d = m.current_data()
    K = s.Matrix.vstack(*d['matter_frame'])
    assert m._zero(K.H*K-d['casimir15'])
    assert len(d['port_coefficients']) == 12
    assert np.max(execution[1]['dressing_finite_covariance_errors']) < 1e-12
    # Omitting the adjoint-index action is a physically different coupling.
    native = np.asarray(d['representation']['matter_currents'][0], complex)
    Kn = np.asarray(K, complex)
    assert np.linalg.norm(np.kron(np.eye(12), native)@Kn-Kn@native) > .1


def test_full_domain_counts_keep_bright_dark_and_neutral_states(execution):
    r = execution[0]['complete_domain']
    assert r['total_dimension'] == 90*(1+2*12)
    assert r['bright_dimension']+sum(x['multiplicity'] for x in r['all_dark_poles']) == 2250
    assert 14*18+2*9 == r['current_image_reduction_dimension']
    assert r['family_dark_dimension_inside_reduction']+r['current_orthogonal_auxiliary_dimension'] == 2040
    assert r['full_operator_zero_multiplicity'] == 6


def test_complete_finite_resolvent_matches_all_poles_and_matrix_residues(execution):
    _, arrays = execution
    z = .125+.03125j
    light = list(range(3))+list(range(9, 12))
    for H, poles, rs in zip(arrays['dressing_charged_hamiltonians'], arrays['dressing_charged_poles'],
                            arrays['dressing_charged_residues'], strict=True):
        actual = np.linalg.inv(z*np.eye(18)-H)[np.ix_(light, light)]
        spectral = sum(R/(z-p) for p, R in zip(poles, rs, strict=True))
        np.testing.assert_allclose(actual, spectral, atol=2e-12, rtol=0)


def test_feshbach_self_energy_keeps_frequency_dependence(execution):
    a = m.architecture(); Q = np.asarray(a['Q'], float); Y = np.asarray(a['Y'], float)
    z = .125+.03125j
    light = list(range(3))+list(range(9, 12))
    for H, row in zip(execution[1]['dressing_charged_hamiltonians'], a['channels'], strict=True):
        bare = np.block([[np.zeros((3, 3)), Y], [Y, np.zeros((3, 3))]])
        C = np.block([[float(row['left'])*Q, np.zeros((3, 3))],
                      [np.zeros((3, 3)), float(row['right'])*Q]])
        effective = np.linalg.inv(z*np.eye(6)-bare-C*z/(z*z-1))
        actual = np.linalg.inv(z*np.eye(18)-H)[np.ix_(light, light)]
        np.testing.assert_allclose(actual, effective, atol=2e-12, rtol=0)
        approximate = np.linalg.inv(z*(np.eye(6)+C)-bare)
        assert np.linalg.norm(actual-approximate) > 1e-4


def test_literal2250_operator_resolvent_keeps_arbitrary_light_reference_columns(execution):
    bare, K = m.full_factors(); full = m.full_sparse_operator()
    assert full.shape == (2250, 2250)
    assert np.max(abs((full-full.conj().T).data), initial=0) < 1e-14
    rng = np.random.default_rng(81)
    inputs = rng.normal(size=(90, 3))+1j*rng.normal(size=(90, 3))
    rhs = np.zeros((2250, 3), complex); rhs[:90] = inputs
    z = .125+.25j
    literal = spsolve(z*sparse.eye(2250, format='csc')-full.tocsc(), rhs)[:90]
    effective = np.linalg.solve(z*np.eye(90)-bare-(K.conj().T@K)*z/(z*z-1), inputs)
    np.testing.assert_allclose(literal, effective, atol=5e-13, rtol=0)


def test_complete_operator_sign_symmetry_is_kept_distinct_from_physical_spin(execution):
    full = m.full_sparse_operator()
    light = sparse.diags([1]*45+[-1]*45, format='csr', dtype=float)
    aux = sparse.kron(sparse.eye(12), light, format='csr')
    gamma = sparse.bmat([[light, None, None], [None, None, -aux], [None, -aux, None]], format='csr')
    residual = gamma@full@ gamma+full
    assert np.max(abs(residual.data), initial=0) < 1e-14
    assert not execution[0]['scope']['full_physical_spin_particle_hole_history_derived']
    assert not execution[0]['scope']['auxiliary_Weyl_kinetics_and_physical_spin_attached']


def test_decomposition_exhausts_full2250_resolvent_trace_and_all_dark_poles(execution):
    _, arrays = execution; bare, K = m.full_factors(); gram = K.conj().T@K
    for z in (.125+.25j, -.5+.375j, 1.5+.5j):
        light = np.linalg.inv(z*np.eye(90)-bare-gram*z/(z*z-1))
        literal_trace = (np.trace(light)+1080*(1/(z-1)+1/(z+1))
            +.5*(1/(z-1)**2+1/(z+1)**2)*np.trace(gram@light))
        decomposed = sum(row['multiplicity']*np.sum(1/(z-poles)) for row, poles in
            zip(m.architecture()['channels'], arrays['dressing_charged_poles'], strict=True))
        decomposed += 2*np.sum(1/(z-arrays['dressing_neutral_poles']))
        decomposed += 990*(1/(z-1)+1/(z+1))
        assert abs(literal_trace-decomposed) < 3e-11


def test_sector_splitting_is_from_one_fixed_rule_not_three_fitted_strengths(execution):
    report, _ = execution
    assert report['common_coupling_exact'] == report['common_auxiliary_hopping_exact'] == '1'
    assert [(str(r['left']), str(r['right'])) for r in m.architecture()['channels']] == [
        ('25/8', '4'), ('25/8', '3'), ('9/8', '3')]
    poles = [r['light_poles'] for r in report['charged_channels']]
    assert len({tuple(x) for x in poles}) == 3
    assert not report['scope']['poles_inferred_from_detector_records_here']
    assert not report['scope']['complete2250_operator_compiled_into_native_source_words']


def test_source_rank_one_decomposition_bounds_exact_pole_ratio_without_tuning(execution):
    bound = execution[0]['exact_no_tuning_shape_bound']
    assert s.Rational(bound['middle_over_smallest_ratio_upper_exact']) < 32
    assert not bound['parameter_scan_or_low_energy_kinetic_approximation_used']
    assert not bound['three_poles_below_threshold_asserted_for_every_Delta']
    assert bound['subthreshold_poles_have_nonzero_light_component']
    assert not bound['above_threshold_reassignment_excluded']
    for row in execution[0]['charged_channels']:
        assert row['light_poles'][1]/row['light_poles'][0] < 32
    a = m.architecture(); P = s.eye(3)-a['Q']; Y = a['Y']
    S = (Y-Y*P*Y/s.trace(P*Y)).applyfunc(s.simplify)
    assert m._zero(S*P) and S.rank() == 2 and (Y-S).rank() == 1


@pytest.mark.parametrize('change', ['missing', 'extra', 'nonfinite', 'dark_pole', 'residue'])
def test_complete_spectral_custody_fails_closed(execution, change):
    arrays = deepcopy(execution[1])
    if change == 'missing': arrays.pop('dressing_neutral_residues')
    elif change == 'extra': arrays['untracked'] = np.zeros(1)
    elif change == 'nonfinite': arrays['dressing_charged_hamiltonians'][0, 0, 0] = np.nan
    elif change == 'dark_pole': arrays['dressing_charged_poles'][0, 4] += .01
    else: arrays['dressing_charged_residues'][0, 0] *= .9
    with pytest.raises(ValueError): m.validate_arrays(arrays)
