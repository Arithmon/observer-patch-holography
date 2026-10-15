"""Independent finite-matrix checks of inverse and history error bounds."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.interacting_reconstruction import bounds
from model.relativistic_scalar.lattice import local_operators


def _coherence(hamiltonian, time, j, k):
    state = expm(-1j*time*hamiltonian)[:, j]
    probabilities = []
    for sign in (-1, 1):
        vector = np.zeros(len(hamiltonian), complex)
        vector[j], vector[k] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
        probabilities.append(abs(np.vdot(vector, state))**2)
    return (probabilities[0]-probabilities[1])/2


@pytest.mark.parametrize("dimension", [2, 4, 8])
def test_richardson_bound_for_independent_real_generators(dimension):
    random = np.random.default_rng(913+dimension)
    matrix = random.normal(size=(dimension, dimension))
    matrix = (matrix+matrix.T)/2
    eigenvalues = np.linalg.eigvalsh(matrix)
    centered = matrix-np.eye(dimension)*(eigenvalues[0]+eigenvalues[-1])/2
    norm_cap = np.linalg.norm(centered, 2)+1e-12
    duration = .2
    estimate = (8*_coherence(matrix, duration/2, 0, 1)
                - _coherence(matrix, duration, 0, 1))/(3*duration)
    assert abs(estimate-matrix[1, 0]) <= bounds.coefficient_radius(duration, norm_cap)
    # A scalar generator shift cannot change the measured coherence.
    assert _coherence(matrix+7*np.eye(dimension), duration, 0, 1) == pytest.approx(
        _coherence(matrix, duration, 0, 1), abs=1e-14)


def test_coefficient_bound_is_exact_rational_then_outward_float():
    t, b, delta = F(1, 256), F(16), F(1, 10**8)
    exact = b**5*t**4/18+3*delta/t
    actual = bounds.coefficient_radius(t, b, delta)
    assert F(actual) >= exact
    assert F(math.nextafter(actual, -math.inf)) < exact
    assert bounds.coefficient_radius(t, 0) == 0


def test_generalized_richardson_radius_reduces_to_dyadic_bound():
    t, b, delta = F(1, 256), F(16), F(1, 10**8)
    certificate = bounds.richardson_radius(t/2, t, b, delta)
    assert certificate["radius"] == bounds.coefficient_radius(t, b, delta)
    assert F(certificate["finite_time_radius_exact"]) == b**5*t**4/18
    assert F(certificate["probability_radius_exact"]) == 3*delta/t


def test_generalized_richardson_controls_non_dyadic_measurements():
    matrix = np.asarray([[1, .2, -.3], [.2, 2, .4], [-.3, .4, -1]])
    u, v = .07, .19
    estimate = (v*v*_coherence(matrix, u, 0, 1)/u
                - u*u*_coherence(matrix, v, 0, 1)/v)/(v*v-u*u)
    certificate = bounds.richardson_radius(u, v, 3)
    assert abs(estimate-matrix[1, 0]) <= certificate["radius"]
    with pytest.raises(ValueError):
        bounds.richardson_radius(v, u, 3)


def test_simultaneous_probability_perturbations_attain_sensitivity():
    # Worst-sign perturbations of the four individually bounded laws.
    t, delta = F(1, 256), F(1, 10**5)
    changes = [delta, -delta, -delta, delta]  # minus/plus at half/full time
    change = (4*(changes[0]-changes[1])-(changes[2]-changes[3])/2)/(3*t)
    assert change == 3*delta/t
    assert float(change) <= bounds.coefficient_radius(t, 0, delta)


@pytest.mark.parametrize("dimension", [1, 2, 4, 8, 16])
@pytest.mark.parametrize("omega", [F(1, 2), F(2), F(7, 2)])
def test_polynomial_norm_caps_use_true_infinite_compressions(dimension, omega):
    operators = local_operators(dimension, float(omega))
    upper = bounds.oscillator_norm_bounds(dimension, omega)
    for name in ("q", "q2", "q4"):
        actual = np.linalg.norm(operators[name], 2)
        assert actual <= float(upper[name])+1e-12
    if dimension == 1:
        assert operators["q"][0, 0] == 0
        assert operators["q4"][0, 0] > 0


def test_hamiltonian_uncertainty_controls_independent_dense_difference():
    dims = (4, 8)
    quadratic, quartic = [F(1, 100), F(1, 200)], [F(1, 300), F(1, 400)]
    bond = F(1, 500)
    result = bounds.parameter_uncertainty_bound(dims, quadratic, quartic,
                                                {(0, 1): bond}, F(1, 8))
    left, right = [local_operators(d, 2) for d in dims]
    delta = (np.kron(np.eye(dims[1]), float(quadratic[0])/2*left["q2"]
                     - float(quartic[0])*left["q4"])
             + np.kron(-float(quadratic[1])/2*right["q2"]
                       + float(quartic[1])*right["q4"], np.eye(dims[0]))
             - float(bond)*np.kron(right["q"], left["q"]))
    assert np.linalg.norm(delta, 2) <= result["hamiltonian_norm_error_upper"]
    assert result["complete_history_tv_upper"] <= result["hamiltonian_norm_error_upper"]/8+1e-15
    assert F(result["hamiltonian_norm_error_upper"]) >= F(result["exact_hamiltonian_norm_bound"])


def test_full_disturbed_history_obeys_parameter_budget():
    # Independent two-site, two-read history, retaining both outcomes and
    # conditional unitary feedback.  It does not assume a vacuum eigenstate.
    left, right = local_operators(4), local_operators(2)
    dimension = 8
    onsite = .5*left["p2"]+2*left["q2"]+left["q4"]/24
    h0 = np.kron(np.eye(2), onsite)+np.kron(right["p2"]/2+right["q2"], np.eye(4))
    h0 -= np.kron(right["q"], left["q"])
    delta_g, delta_k = .002, .003
    h1 = h0+delta_g*np.kron(np.eye(2), left["q4"])-delta_k*np.kron(right["q"], left["q"])
    e = np.zeros(dimension, complex)
    e[0], e[2] = 1/math.sqrt(2), 1j/math.sqrt(2)
    yes = np.outer(e, e.conj())
    projectors = [np.eye(dimension)-yes, yes]
    feedback = np.diag([1, 1, 1j, 1, 1, 1, 1, 1])
    durations = (.13, .17)

    def history(h):
        initial = np.eye(dimension, dtype=complex)[:, 0]
        first = expm(-1j*durations[0]*h)@initial
        second_u = expm(-1j*durations[1]*h)
        probabilities = []
        for bit, projector in enumerate(projectors):
            selected = projector@first
            if bit:
                selected = feedback@selected
            continued = second_u@selected
            for final in projectors:
                probabilities.append(float(np.linalg.norm(final@continued)**2))
        return np.asarray(probabilities)

    p0, p1 = history(h0), history(h1)
    assert sum(p0) == pytest.approx(1)
    assert sum(p1) == pytest.approx(1)
    certificate = bounds.parameter_uncertainty_bound(
        (4, 2), [0, 0], [delta_g, 0], [(0, 1, delta_k)], sum(durations))
    assert .5*np.sum(abs(p0-p1)) <= certificate["complete_history_tv_upper"]


def test_zero_uncertainty_and_probability_cap():
    zero = bounds.parameter_uncertainty_bound((4, 2), [0, 0], [0, 0], [], 1)
    assert zero["hamiltonian_norm_error_upper"] == 0
    assert zero["complete_history_tv_upper"] == 0
    large = bounds.parameter_uncertainty_bound((4, 2), [10, 10], [1, 1], [], 1)
    assert large["complete_history_tv_upper"] == 1


def test_hoeffding_budget_meets_union_bound():
    delta, laws, alpha = F(1, 100), 176, F(1, 20)
    result = bounds.simultaneous_shot_count(delta, laws, alpha)
    n = result["trials_per_law"]
    assert 2*laws*math.exp(-2*n*float(delta)**2) <= float(alpha)
    assert result["total_fresh_trials"] == laws*n
    assert result["sampling_executed"] is False


@pytest.mark.parametrize("time,norm,error", [(0, 1, 0), (-1, 1, 0), (1, -1, 0),
                                             (1, 1, -1), (True, 1, 0), (1, math.nan, 0)])
def test_invalid_derivative_bound_inputs(time, norm, error):
    with pytest.raises(ValueError):
        bounds.coefficient_radius(time, norm, error)


@pytest.mark.parametrize("pairs", [[(0, 0, 1)], [(1, 0, 1)], [(0, 2, 1)],
                                    [(0, 1, -1)], [(0, 1, 1), (0, 1, 1)]])
def test_invalid_uncertainty_edges(pairs):
    with pytest.raises(ValueError):
        bounds.parameter_uncertainty_bound((4, 4), [0, 0], [0, 0], pairs, 1)


@pytest.mark.parametrize("margin,laws,failure", [(0, 1, .05), (1, 1, .05),
                                                 (.1, 0, .05), (.1, True, .05),
                                                 (.1, 1, 0), (.1, 1, 1)])
def test_invalid_trial_budget(margin, laws, failure):
    with pytest.raises(ValueError):
        bounds.simultaneous_shot_count(margin, laws, failure)
