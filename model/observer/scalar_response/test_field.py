"""Independent dense action and finite-oscillator checks of the field adapter."""
from decimal import Decimal, localcontext
from math import comb, factorial, sqrt

import numpy as np
import pytest
from scipy.linalg import eigh, expm

from model.observer.scalar_response.field import CoherentSum, ScalarField


@pytest.fixture(scope="module")
def dense_q13():
    """Assemble all 1728 modes directly from edge energies and compact profiles.

    This uses neither the adapter's coefficient recurrence nor its tensor
    eigensolver. A closed polynomial profile supplies the same smearing.
    """
    q = 13
    with localcontext() as context:
        context.prec = 70
        phi = (1 + Decimal(5).sqrt()) / 2
        nodes = sorted((b * phi) % 1 for b in range(q)) + [Decimal(1)]
        gaps = np.array([float(b-a) for a, b in zip(nodes, nodes[1:])])
        x = np.array([float(a) for a in nodes[1:-1]])
    mass = (gaps[:-1] + gaps[1:]) / 2
    n = len(x)
    boundary = np.zeros((q, n))
    for edge in range(q):
        if edge:
            boundary[edge, edge-1] = -1
        if edge < n:
            boundary[edge, edge] = 1
    stiffness = boundary.T @ np.diag(1/gaps) @ boundary
    whitened = stiffness / np.sqrt(np.outer(mass, mass))
    one = np.eye(n)
    action = (np.eye(n**3) + np.kron(whitened, np.eye(n*n))
              + np.kron(one, np.kron(whitened, one))
              + np.kron(np.eye(n*n), whitened))
    values, vectors = eigh(action, check_finite=False, driver="evd")
    omega = np.sqrt(values)
    # Integral_0^1 [sin(pi*x)(1+cos(pi*x))^8]^2 dx, evaluated by beta integral.
    normalization = sqrt(comb(34, 17)/(18*2**17))
    profile = lambda z: np.sin(np.pi*z)*(1+np.cos(np.pi*z))**8/normalization
    prep = np.sqrt(mass)*profile(x)*(x < .45)
    detector = np.sqrt(mass)*profile(1-x)*(x > .55)
    transverse = np.sqrt(2*mass)*np.sin(np.pi*x)
    f = vectors.T @ np.kron(prep, np.kron(transverse, transverse))
    g = vectors.T @ np.kron(detector, np.kron(transverse, transverse))
    return omega, f, g


@pytest.mark.parametrize("time", (19/20, 39/40, 1.))
def test_all_mode_tensor_response_matches_independent_dense_action(dense_q13, time):
    omega, f, g = dense_q13
    mean = 4*np.dot(f*g, np.sin(time*omega)/omega)
    variance = np.dot(g*g, 1/(2*omega))
    plus = (1+np.exp(-variance/2)*np.sin(mean))/2
    field = ScalarField(13)
    actual = field.response(time)
    assert actual["mode_count"] == len(omega) == 1728
    assert actual["mean"] == pytest.approx(mean, rel=0, abs=2e-12)
    assert actual["vacuum_variance"] == pytest.approx(variance, rel=0, abs=2e-13)
    assert actual["plus_probability"] == pytest.approx(plus, rel=0, abs=1e-12)
    assert field.response(time, current=0)["plus_probability"] == .5
    assert field.response(time, current=-1)["plus_probability"] == pytest.approx(1-plus, abs=1e-12)


def test_terminal_branch_gram_matches_explicit_coherent_instrument():
    field = ScalarField(13)
    modes = field.mode_data()
    for current in (-1, 0, 1):
        initial = field.prepared_state(39/40, current=current)
        probabilities = []
        for outcome in (-1, 1):
            probability, branch = initial.centered_instrument(modes["field_vector"], outcome)
            descriptor = field.branch_descriptor(39/40, outcome, current=current)
            raw = np.asarray(descriptor["gram"])
            gram = raw[:, :, 0] + 1j*raw[:, :, 1]
            raw = np.asarray(descriptor["coefficients"])
            coefficients = raw[:, 0] + 1j*raw[:, 1]
            assert probability == pytest.approx(descriptor["probability"], abs=3e-15)
            assert coefficients.conj() @ gram @ coefficients == pytest.approx(1, abs=3e-15)
            assert branch.norm_squared() == pytest.approx(1, abs=3e-15)
            probabilities.append(probability)
        assert sum(probabilities) == pytest.approx(1, abs=3e-15)


def test_zero_preparation_and_zero_pointer_coupling_controls():
    field = ScalarField(13)
    assert field.response(39/40, gamma=0)["plus_probability"] == .5
    source = CoherentSum([[.2+.3j], [-.5+.1j]], [1., .4-.2j]).normalized()
    for outcome in (-1, 1):
        probability, branch = source.centered_instrument([0.], outcome)
        assert probability == pytest.approx(.5, abs=2e-15)
        # Zero coupling must leave the complete input ray unchanged.
        assert abs(source.inner(branch))**2 == pytest.approx(1, abs=3e-15)


def fock_vector(alpha, cutoff=40):
    return np.array([np.exp(-abs(alpha)**2/2)*alpha**n/sqrt(factorial(n))
                     for n in range(cutoff)])


def materialize(state, cutoff=40):
    assert state.amplitudes.shape[1] == 1
    return sum(coefficient*fock_vector(alpha[0], cutoff)
               for coefficient, alpha in zip(state.coefficients, state.amplitudes))


@pytest.mark.parametrize("outcome", (-1, 1))
def test_coherent_backaction_against_independent_fock_kraus_matrices(outcome):
    # A non-Gaussian input makes an incorrect "update the Gaussian mean"
    # implementation fail even if its first-readout probability is right.
    cutoff = 40
    source = CoherentSum([[.2+.3j], [-.5+.1j]], [1., .4-.2j]).normalized()
    annihilation = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    beta = .7
    position = beta*(annihilation+annihilation.T)
    kraus = (expm(-.5j*position)-1j*outcome*expm(.5j*position))/2
    expected = kraus @ materialize(source, cutoff)
    probability = np.vdot(expected, expected).real
    actual_probability, actual_state = source.centered_instrument([beta], outcome)
    assert actual_probability == pytest.approx(probability, rel=0, abs=3e-14)
    np.testing.assert_allclose(materialize(actual_state, cutoff), expected/np.sqrt(probability),
                               rtol=0, atol=3e-14)

    # Continue the disturbed branch through free evolution and another
    # noncommuting readout. This checks the retained relative Weyl phases.
    omega, time = 1.3, .47
    unitary = np.diag(np.exp(-1j*omega*time*np.arange(cutoff)))
    continued = unitary @ (expected/np.sqrt(probability))
    second = (expm(-.2j*position)-1j*expm(.2j*position))/2
    expected_second = second @ continued
    branch = actual_state.evolve([omega], time)
    second_probability, _ = branch.centered_instrument([.4*beta], 1)
    assert second_probability == pytest.approx(np.vdot(expected_second, expected_second).real,
                                               rel=0, abs=4e-14)
