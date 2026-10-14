"""Independent coherent histories and hostile inputs for the cutoff bound."""
from fractions import Fraction as F
from math import factorial, sqrt

import numpy as np
import pytest
from scipy.linalg import expm

from model.observer_reconstruction.bounds import coherent_history_cutoff


def test_declared_eight_mode_budgets_are_exact_and_decrease():
    # |v_xj|=1/sqrt(8), omega_j>=5/4, gamma+eta=7/10.
    a2 = F(49, 2000)
    budgets = [coherent_history_cutoff(d, amplitude_squared_upper=a2)
               for d in (8, 16)]
    for d, report in zip((8, 16), budgets):
        low = F(d, factorial(d))*a2**(d+1)
        tail = a2**d/(factorial(d)*(1-a2/F(d+1)))
        square = report["single_mode_vector_error_squared_upper"]
        root = report["single_mode_vector_error_upper"]
        assert square == low + tail
        assert type(square) is F and type(root) is F
        assert root**2 >= square
        assert (root-F(1, 10**36))**2 < square
        assert report["joint_vector_error_upper"] == 8*root
        assert report["pointer_mean_error_upper"] == 16*root
        assert not report["scope"]["numerical_source_execution_error_included"]
        assert not report["scope"]["inverse_reconstruction_error_included"]
    assert budgets[0]["pointer_probability_error_upper"] < F(2, 10**8)
    assert budgets[1]["pointer_probability_error_upper"] < F(3, 10**19)
    assert budgets[1]["pointer_probability_error_upper"] < budgets[0]["pointer_probability_error_upper"]


def test_zero_displacement_is_exact_and_large_bounds_clip_only_trace_distance():
    zero = coherent_history_cutoff(8, amplitude_squared_upper=F(0))
    assert zero["joint_vector_error_upper"] == 0
    assert zero["pointer_probability_error_upper"] == 0
    large = coherent_history_cutoff(2, amplitude_squared_upper=F(2), modes=8)
    assert large["joint_vector_error_upper"] > 1
    assert large["complete_output_trace_distance_upper"] == 1
    assert large["pointer_mean_error_upper"] == 2


@pytest.mark.parametrize("dimension", [2, 4, 8])
@pytest.mark.parametrize("theta", [0.0, 0.8, 2.2])
def test_bound_contains_full_controlled_finite_vs_infinite_history(dimension, theta):
    # A deliberately larger amplitude than the retained benchmark keeps
    # truncation visible.  Both coherent pointer arms retain their phases.
    kick, probe = 0.35j, 0.45j
    lowering = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    raising = lowering.T
    rotation = np.diag(np.exp(-1j*theta*np.arange(dimension)))
    vacuum = np.eye(dimension, dtype=complex)[:, 0]
    prepared = rotation @ expm(kick*raising-kick.conjugate()*lowering) @ vacuum
    finite = np.concatenate((prepared,
        expm(probe*raising-probe.conjugate()*lowering) @ prepared))/sqrt(2)

    # The infinite comparator is a coherent-state series, independently
    # computed without a larger truncated Hamiltonian or exponentiation.
    first_alpha = kick*np.exp(-1j*theta)
    amplitudes = (first_alpha, first_alpha+probe)
    phases = (1, np.exp((probe*first_alpha.conjugate()-probe.conjugate()*first_alpha)/2))
    low = []
    retained_mass = 0
    for alpha, phase in zip(amplitudes, phases):
        coefficients = phase*np.exp(-abs(alpha)**2/2)*np.array(
            [alpha**n/sqrt(factorial(n)) for n in range(dimension)])/sqrt(2)
        low.append(coefficients)
        retained_mass += np.vdot(coefficients, coefficients).real
    projected = np.concatenate(low)
    squared_error = np.linalg.norm(finite-projected)**2 + max(0.0, 1-retained_mass)
    report = coherent_history_cutoff(dimension, amplitude_squared_upper=F(16, 25), modes=1)
    assert squared_error <= float(report["single_mode_vector_error_squared_upper"])+2e-15

    # Ramsey X and Y readouts include interference between the retained arms.
    finite_characteristic = np.vdot(finite[:dimension], finite[dimension:])*2
    infinite_characteristic = np.exp(-abs(probe)**2/2 +
        probe*first_alpha.conjugate()-probe.conjugate()*first_alpha)
    probability_budget = float(report["pointer_probability_error_upper"])
    assert abs(finite_characteristic.real-infinite_characteristic.real)/2 <= probability_budget+2e-15
    assert abs(finite_characteristic.imag-infinite_characteristic.imag)/2 <= probability_budget+2e-15


@pytest.mark.parametrize("dimension", [True, 1, 0, -2, 8.0, "8"])
def test_rejects_invalid_dimensions(dimension):
    with pytest.raises(ValueError, match="dimension"):
        coherent_history_cutoff(dimension, amplitude_squared_upper=F(1, 10))


@pytest.mark.parametrize("amplitude", [True, 0, 0.1, "1/10", F(-1), F(9), F(10)])
def test_rejects_invalid_amplitude_squared_budgets(amplitude):
    with pytest.raises(ValueError, match="amplitude"):
        coherent_history_cutoff(8, amplitude_squared_upper=amplitude)


@pytest.mark.parametrize("modes", [True, 0, -2, 8.0, "8"])
def test_rejects_invalid_mode_counts(modes):
    with pytest.raises(ValueError, match="modes"):
        coherent_history_cutoff(8, amplitude_squared_upper=F(1, 10), modes=modes)


@pytest.mark.parametrize("digits", [True, 0, -2, 8.0, "8"])
def test_rejects_invalid_precision(digits):
    with pytest.raises(ValueError, match="precision"):
        coherent_history_cutoff(8, amplitude_squared_upper=F(1, 10), digits=digits)
