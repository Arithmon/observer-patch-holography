"""Independent unsplit dynamics checks using synthetic preparations only."""
import ast
import inspect
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import reference, source


SETTINGS = [{'u': u, 'v': v, 'theta': phase}
            for u, v in ((0., 0.), (0., .37), (.23, -.41))
            for phase in (0., -math.pi/2)]


def _observables(records):
    return np.asarray([value for record in records
                       for value in [*record['spin_probabilities'],
                                     record['selected_position_mean'],
                                     *[read['selected_mean'] for read in record['reads']]]])


def test_reference_does_not_import_source_or_compiler():
    tree = ast.parse(inspect.getsource(reference))
    imports = [node for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))]
    for node in imports:
        names = [alias.name for alias in node.names]
        if isinstance(node, ast.ImportFrom):
            assert node.level == 0
            names.append(node.module or '')
        assert not any('source' in name or 'compiler' in name for name in names)


def test_unsplit_reference_matches_dense_liouvillian():
    dimension = 3
    rho = reference.preparation(dimension, -.48, .21)
    a, number, spin_x = reference.operators(dimension)
    joint_number = np.kron(np.eye(2), number)
    hamiltonian = 7.2*joint_number-3.1/2*np.kron(spin_x, a+a.T)
    collapse = math.sqrt(.8)*joint_number
    loss = collapse.T @ collapse
    identity = np.eye(2*dimension)
    # Column-major vectorization is independent of the ODE's row-major RHS.
    liouvillian = (-1j*(np.kron(identity, hamiltonian)-np.kron(hamiltonian.T, identity))
                   +np.kron(collapse.conj(), collapse)
                   -(np.kron(identity, loss)+np.kron(loss.T, identity))/2)
    expected = (expm(.017*liouvillian) @ rho.ravel(order='F')).reshape(rho.shape, order='F')
    actual = reference.lindblad_pulse(rho, delta=7.2, omega=3.1, duration=.017,
                                     sign=1, dephasing_rate=.8)
    np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0)


def test_diffusion_rate_is_half_squared_number_difference():
    rng = np.random.default_rng(741)
    vector = rng.normal(size=8)+1j*rng.normal(size=8)
    vector /= np.linalg.norm(vector)
    rho = np.outer(vector, vector.conj())
    numbers = np.tile(np.arange(4), 2)
    differences = numbers[:, None]-numbers[None, :]
    expected = rho*np.exp(-.7*.09*differences**2/2)
    unsplit = reference.lindblad_pulse(rho, delta=0., omega=0., duration=.09,
                                      sign=1, dephasing_rate=.7)
    np.testing.assert_allclose(unsplit, expected, atol=2e-12, rtol=0)
    np.testing.assert_allclose(source.phase_diffusion(rho, 4, .7, .09), expected,
                               atol=2e-16, rtol=0)
    # Spin coherence at equal occupation is preserved by motional dephasing.
    assert unsplit[0, 4] == pytest.approx(rho[0, 4], abs=2e-13)


@pytest.mark.parametrize('dimension', [4, 8])
def test_source_unitary_path_matches_unsplit_zero_noise_reference(dimension):
    arguments = (dimension, -.61, .14, [0, 1, 3], SETTINGS)
    actual = source.run(*arguments, dephasing_rate=0., subdivisions=1)
    expected = reference.run(*arguments, dephasing_rate=0.)
    np.testing.assert_allclose(_observables(actual['records']), _observables(expected),
                               atol=3e-10, rtol=0)
    for got, want in zip(actual['records'], expected, strict=True):
        for got_read, want_read in zip(got['reads'], want['reads'], strict=True):
            np.testing.assert_allclose(got_read['spin_detector_joint_probabilities'],
                                       want_read['spin_detector_joint_probabilities'],
                                       atol=3e-10, rtol=0)


@pytest.mark.parametrize('dimension', [4, 8])
def test_symmetric_diffusion_split_converges_to_unsplit_reference(dimension):
    arguments = (dimension, -.61, .14, [1, 3], SETTINGS)
    # Amplify the known supplied noise channel to resolve splitting convergence.
    expected = _observables(reference.run(*arguments, dephasing_rate=400.))
    errors = []
    for subdivisions in (1, 2, 4):
        actual = source.run(*arguments, dephasing_rate=400., subdivisions=subdivisions)
        errors.append(np.max(abs(_observables(actual['records'])-expected)))
    assert errors[0] > 1e-6
    assert errors[1] < .30*errors[0]
    assert errors[2] < .30*errors[1]


def test_actual_published_noise_setting_matches_refined_source():
    arguments = (4, -.77, .09, [2, 10], SETTINGS)
    expected = _observables(reference.run(*arguments))
    actual = _observables(source.run(*arguments, subdivisions=8)['records'])
    refined = _observables(source.run(*arguments, subdivisions=16)['records'])
    # Eight slices are an approximation; the acceptance mass is the largest
    # difference here. Check its declared numerical refinement, not identity.
    error = np.max(abs(actual-expected))
    refined_error = np.max(abs(refined-expected))
    assert error < 3e-5
    assert refined_error < .30*error


def test_final_selection_differs_from_repeated_spin_projection():
    arguments = (4, -.61, .14, [3], SETTINGS)
    retained = reference.run(*arguments, dephasing_rate=0.)[0]
    reset = reference.preparation(4, -.61, .14)
    for _ in range(3):
        reset = reference.step(reset, dephasing_rate=0.)
        selected = reset[:4, :4].copy()
        reset[:] = 0
        reset[:4, :4] = selected/np.trace(selected)
    repeated = reference.readout(reset, SETTINGS)
    assert abs(retained['selected_position_mean']-repeated['selected_position_mean']) > 1e-3
    assert retained['spin_probabilities'][1] > 1e-3


def test_small_rotation_has_positive_cosine_potential_on_down_zero():
    dimension = 8
    a, _, spin_x = reference.operators(dimension)
    alpha = .3
    theta = 1e-5
    spin_z = np.diag([1., -1.])
    displacement = expm(1j*alpha*np.kron(spin_x, a+a.T))
    rotation = np.kron(expm(-.5j*theta*spin_z), np.eye(dimension))
    gate = displacement.conj().T @ rotation @ displacement @ displacement @ rotation @ displacement.conj().T
    expected = (expm(2j*alpha*(a+a.T))+expm(-2j*alpha*(a+a.T)))/2
    np.testing.assert_allclose((gate[:dimension, :dimension]-np.eye(dimension))/theta,
                               -1j*expected, atol=theta, rtol=0)
