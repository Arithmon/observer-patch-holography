"""Actual-history observer reconstruction and discriminating modular controls."""
import numpy as np
import pytest
from scipy.linalg import expm

from . import experiment as e
from . import operational as o


@pytest.fixture(scope='module', params=(2, 4))
def history(request):
    return e.history(request.param, 16)


def test_gns_reconstructs_positive_functional_and_left_representation():
    rng = np.random.default_rng(891)
    for rank in (1, 2, 3):
        source = rng.normal(size=(3, rank))+1j*rng.normal(size=(3, rank))
        rho = source@source.conj().T
        rho /= np.trace(rho)
        gns = o.MatrixGNS.from_density(rho)
        a, b = (rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)) for _ in range(2))
        assert gns.dimension == 3*rank
        assert gns.inner(a, b) == pytest.approx(np.trace(rho@a.conj().T@b), abs=1e-13)
        np.testing.assert_allclose(gns.lift(a@b), a@gns.lift(b), atol=1e-13)
        assert gns.expectation(a) == pytest.approx(np.trace(rho@a), abs=1e-13)


def test_pure_state_gns_is_not_just_its_one_dimensional_support_corner():
    rho = np.diag([1., 0.])
    gns = o.MatrixGNS.from_density(rho)
    raising = np.array([[0., 0.], [1., 0.]])
    assert gns.dimension == 2
    assert gns.inner(raising, raising) == pytest.approx(1.)
    with pytest.raises(ValueError, match='faithful'):
        o.modular_flow(rho, raising, .5)


def test_same_history_readout_marginals_and_explicit_access_choices(history):
    receipt = o.observation_receipt(history)
    d = history['dimension']
    assert receipt['all_pass']
    assert receipt['gns_hilbert_dimension'] == d*d
    assert receipt['record_and_work_block_algebra_gns_dimension'] == 16*d*d
    assert receipt['records_only_gns_dimension'] == 4
    assert receipt['records_only_modular_flow'] == 'identity'
    assert not receipt['clock_calibration_established']
    rho, _ = o.scalar_state(history)
    weighted = np.zeros_like(rho)
    for row in history['branches']:
        state, probability = o.scalar_state(history, record=row['record'])
        weighted += probability*state
        assert o.MatrixGNS.from_density(state).rank <= d//2
        with pytest.raises(ValueError, match='faithful'):
            o.modular_hamiltonian(state)
    np.testing.assert_allclose(weighted, rho, atol=1e-13)


def test_partial_trace_pairing_against_full_working_observables(history):
    d = history['dimension']
    rng = np.random.default_rng(893+d)
    probe = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
    vectors = o.terminal_vectors(history)
    for site in (0, 1):
        rho, _ = o.scalar_state(history, site)
        embedded = (np.kron(np.kron(np.eye(d), probe), np.eye(4)) if site == 0
                    else np.kron(np.kron(probe, np.eye(d)), np.eye(4)))
        direct = sum(np.vdot(v, embedded@v) for v in vectors.values())
        assert np.trace(rho@probe) == pytest.approx(direct, abs=1e-12)


def test_both_observer_gns_spaces_embed_into_one_history_representation(history):
    from .spatial import scalar_operator
    d = history['dimension']
    rng = np.random.default_rng(901+d)
    x = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
    for site in (0, 1):
        inclusion, local, common_cyclic = o.regional_gns_embedding(history, site)
        np.testing.assert_allclose(inclusion.conj().T@inclusion, np.eye(d*d), atol=2e-11)
        np.testing.assert_allclose(inclusion@local.root_columns.ravel(), common_cyclic, atol=2e-12)
        actual = np.concatenate([scalar_operator(d, site, x)@block
                                 for block in inclusion.reshape(4, 4*d*d, d*d)], axis=0)
        np.testing.assert_allclose(actual, inclusion@np.kron(x, np.eye(d)), atol=2e-11)


def test_modular_group_state_and_kms_identity_on_actual_state(history):
    rho, _ = o.scalar_state(history)
    d = len(rho)
    rng = np.random.default_rng(897+d)
    a, b = (rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)) for _ in range(2))
    flow = lambda x, s: o.modular_flow(rho, x, s)
    np.testing.assert_allclose(flow(flow(a, .3), -.8), flow(a, -.5), atol=1e-12)
    np.testing.assert_allclose(flow(a@b, .3), flow(a, .3)@flow(b, .3), atol=1e-12)
    np.testing.assert_allclose(flow(a.conj().T, .3), flow(a, .3).conj().T, atol=1e-12)
    assert np.trace(rho@flow(a, .3)) == pytest.approx(np.trace(rho@a), abs=1e-12)
    # Analytic continuation sigma_-i(B)=rho B rho^-1 fixes the sign convention.
    sigma_minus_i_b = np.linalg.solve(rho.T, (rho@b).T).T
    assert np.trace(rho@a@sigma_minus_i_b) == pytest.approx(np.trace(rho@b@a), abs=1e-11)


def test_record_aware_scalar1_modular_observer_keeps_all_committed_bits(history):
    result = o.record_aware_scalar_view(history)
    assert result['gns_hilbert_dimension'] == 4*history['dimension']**2
    assert {tuple(row['record']) for row in result['blocks']} == {(0, 0), (0, 1), (1, 0), (1, 1)}
    assert sum(row['probability'] for row in result['blocks']) == pytest.approx(1., abs=1e-12)
    assert all(row['number_probe_moves'] for row in result['blocks'])


def test_gibbs_time_correspondence_requires_the_extra_thermal_state_hypothesis():
    h = np.array([[.4, .7j], [-.7j, 1.3]])
    beta, s = .8, .6
    rho = expm(-beta*h)
    rho /= np.trace(rho)
    probe = np.diag([1., -1.])
    model_time = -beta*s
    physical = expm(1j*model_time*h)@probe@expm(-1j*model_time*h)
    np.testing.assert_allclose(o.modular_flow(rho, probe, s), physical, atol=1e-13)
    # A faithful tracial state supplies no nontrivial modular clock.
    np.testing.assert_allclose(o.modular_flow(np.eye(2)/2, probe, 123.), probe, atol=1e-13)


@pytest.mark.parametrize('rho', (np.diag([1.2, -.2]), np.eye(2), [[1, 1], [0, 0]]))
def test_rejects_nonstates_instead_of_repairing_them(rho):
    with pytest.raises(ValueError):
        o.MatrixGNS.from_density(rho)
