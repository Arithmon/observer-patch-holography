"""Complete native incidence, coherent lift, covariance and all bath poles."""
from copy import deepcopy

import numpy as np
import pytest
import sympy as s
from scipy import sparse
from scipy.sparse.linalg import spsolve

from . import native_mass_bath as m


@pytest.fixture(scope='module')
def execution(): return m.execute()


def test_exact_native_gram_retains_all_inputs_and_timeout(execution):
    a = m.source()
    assert a['F'].T*a['F'] == a['D'] == s.eye(12)-a['T']
    assert a['D'] == a['D'].T and a['D'].rank() == 9
    assert a['timeout'] == s.Rational(5, 32)
    assert a['F'][20, :] == s.zeros(1, 12)
    assert execution[0]['source']['complete_reduced_one_attempt_records'] == 2706
    assert not execution[0]['source']['timeout_outcomes_conditioned_away']


def test_source_anisotropy_is_exact_and_kept_outside_closed_compression(execution):
    a = m.source()
    assert a['commutator_squared'] == (s.Integer(4065763817)-1618423471*s.sqrt(5))/52399067238400
    assert a['commutator_squared'].is_positive
    assert not execution[0]['fixed_rule']['family_domain_compression_used_as_closed_bath']
    assert execution[0]['fixed_rule']['coherent_Hermitian_lift_is_an_added_repair_service_choice']


def test_full_domain_zero_modes_and_padding_are_not_removed(execution):
    report, arrays = execution; d = report['complete_domain']
    assert d['total_dimension'] == 90+12*30*64 == 23130
    assert d['current_image_reduction_dimension'] == 14*134+2*67
    assert d['current_orthogonal_dimension'] == 330*64
    assert d['full_zero_multiplicity'] == 14*88+2*45+330*46
    assert d['full_zero_multiplicity']+2*d['full_positive_pole_multiplicity'] == d['total_dimension']
    assert np.count_nonzero(arrays['bath_hamiltonian'][33:]) == 0
    assert np.count_nonzero(arrays['bath_family_embedding'][12:]) == 0
    assert np.max(abs(arrays['bath_charged_residues'][:, 23])) < 1e-18
    assert np.trace(arrays['bath_neutral_residues'][11]) > .02


def test_exact_current_covariance_is_the_same_full_auxiliary_intertwiner(execution):
    data = m.dressing.current_data(); K = np.vstack(execution[1]['bath_current_frame'])
    np.testing.assert_allclose(K.conj().T@K, np.asarray(data['casimir15'], complex), atol=1e-13, rtol=0)
    for rho, adjoint in zip(data['representation']['matter_currents'], data['adjoints'], strict=True):
        rho = np.asarray(rho, complex); adjoint = np.asarray(adjoint, complex)
        residual = (np.kron(adjoint, np.eye(15))+np.kron(np.eye(12), rho))@K-K@rho
        assert np.linalg.norm(residual, 2) < 1e-12


def test_full_charged_resolvent_equals_frequency_dependent_source_schur(execution):
    _, arrays = execution; z = .0625+.125j
    W, D = arrays['bath_family_embedding'], arrays['bath_dirichlet']
    Rb = np.linalg.inv(z*np.eye(64)-arrays['bath_hamiltonian'])
    S = W.T@Rb@W
    np.testing.assert_allclose(S, W[:12].T@np.linalg.solve(z*z*np.eye(12)-D, z*W[:12]), atol=1e-13, rtol=0)
    Y = np.asarray(m.feedback.architecture()['Y'], float)
    bare = np.block([[np.zeros((3, 3)), Y], [Y, np.zeros((3, 3))]])
    for i, row in enumerate(m.dressing.architecture()['channels']):
        self_energy = np.block([[float(row['left'])*S, np.zeros((3, 3))],
                                [np.zeros((3, 3)), float(row['right'])*S]])
        expected = np.linalg.inv(z*np.eye(6)-bare-self_energy)
        actual = np.linalg.inv(z*np.eye(134)-arrays['bath_charged_hamiltonians'][i])[:6, :6]
        np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=0)
        grouped_poles = [float(np.mean(arrays['bath_charged_poles'][i, list(g)])) for g in m.groups(23, 88)]
        spectral = sum(R/(z-e) for R, e in zip(arrays['bath_charged_residues'][i], grouped_poles, strict=True))
        np.testing.assert_allclose(actual, spectral, atol=1e-13, rtol=0)


def test_literal23130_operator_retains_unknown_light_reference_inputs(execution):
    arrays = execution[1]; H = m.full_sparse_operator(arrays)
    assert H.shape == (23130, 23130)
    assert np.max(abs((H-H.conj().T).data), initial=0) < 1e-14
    rng = np.random.default_rng(18); inputs = rng.normal(size=(90, 2))+1j*rng.normal(size=(90, 2))
    rhs = np.zeros((23130, 2), complex); rhs[:90] = inputs; z = .125+.25j
    literal = spsolve(z*sparse.eye(23130, format='csc')-H.tocsc(), rhs)[:90]
    W = arrays['bath_family_embedding']; S = W.T@np.linalg.solve(z*np.eye(64)-arrays['bath_hamiltonian'], W)
    C = np.kron(np.eye(2), np.asarray(m.dressing.current_data()['casimir15'], float))
    expected = np.linalg.solve(z*np.eye(90)-arrays['bath_bare_response']-np.kron(C, S), inputs)
    np.testing.assert_allclose(literal, expected, atol=2e-12, rtol=0)


def test_complete_trace_resolves_every_light_and_dark_pole(execution):
    arrays = execution[1]; W = arrays['bath_family_embedding']
    C = np.kron(np.eye(2), np.asarray(m.dressing.current_data()['casimir15'], float))
    for z in (.0625+.125j, .5+.375j):
        Rb = np.linalg.inv(z*np.eye(64)-arrays['bath_hamiltonian'])
        light = np.linalg.inv(z*np.eye(90)-arrays['bath_bare_response']-np.kron(C, W.T@Rb@W))
        literal = np.trace(light)+360*np.trace(Rb)+np.trace(light@np.kron(C, W.T@Rb@Rb@W))
        decomposed = sum(row['multiplicity']*np.sum(1/(z-poles)) for row, poles in zip(
            m.dressing.architecture()['channels'], arrays['bath_charged_poles'], strict=True))
        decomposed += 2*np.sum(1/(z-arrays['bath_neutral_poles']))+330*np.sum(1/(z-arrays['bath_poles']))
        assert abs(literal-decomposed) < 1e-9


def test_predetermined_window_keeps_all_five_charged_modes(execution):
    report = execution[0]
    assert not report['comparison_window']['three_charged_modes_selected_by_this_window']
    assert [row['frozen_window_pole_count'] for row in report['charged_channels']] == [5, 5, 5]
    assert all(min(row['frozen_window_light_residue_traces']) > 0 for row in report['charged_channels'])
    assert not report['scope']['rule_or_window_tuned_after_spectrum']
    assert not report['scope']['coherent_bath_words_or_detector_records_executed']


@pytest.mark.parametrize('change', ['missing', 'extra', 'nan', 'dark_pole', 'residue'])
def test_full_bath_custody_fails_closed(execution, change):
    arrays = deepcopy(execution[1])
    if change == 'missing': arrays.pop('bath_incidence')
    elif change == 'extra': arrays['hidden'] = np.zeros(1)
    elif change == 'nan': arrays['bath_hamiltonian'][0, 0] = np.nan
    elif change == 'dark_pole': arrays['bath_poles'][30] += .01
    else: arrays['bath_neutral_residues'][11] *= .9
    with pytest.raises(ValueError): m.validate_arrays(arrays)
