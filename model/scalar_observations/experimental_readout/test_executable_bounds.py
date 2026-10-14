"""Independent checks of the displaced-thermal occupation certificate."""
from fractions import Fraction as F
import cmath
import math

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.special import eval_laguerre

from .executable_bounds import cutoff_certificate, thermal_moment, thermal_tail


def test_capped_preparation_moves_exact_geometric_tail():
    dimension, ratio = 12, F(3, 17)
    full_inside = [(1-ratio)*ratio**n for n in range(dimension)]
    capped = full_inside[:-1]+[ratio**(dimension-1)]
    assert sum(capped) == 1
    outside = ratio**dimension
    assert (sum(abs(left-right) for left, right in zip(capped, full_inside))
            + outside)/2 == outside


@pytest.mark.parametrize("thermal,radius,q", [(F(0), F(9, 8), F(5, 2)),
                                             (F(1, 8), F(3, 2), F(2)),
                                             (F(1, 25), F(8), F(3, 2))])
def test_number_moment_bounds_independent_photon_distribution(thermal, radius, q):
    nbar, square, weight = float(thermal), float(radius), float(q)
    probabilities = []
    for number in range(100):
        if nbar == 0:
            probability = math.exp(-square)*square**number/math.factorial(number)
        else:
            ratio = nbar/(1+nbar)
            probability = (math.exp(-square/(1+nbar))/(1+nbar)*ratio**number
                           *eval_laguerre(number, -square/(nbar*(1+nbar))))
        probabilities.append(probability)
    moment = float(thermal_moment(q, thermal, radius))
    assert sum(p*weight**n for n, p in enumerate(probabilities)) <= moment*(1+3e-14)
    assert sum(probabilities) == pytest.approx(1, abs=3e-14)
    for dimension in (8, 16, 32):
        certificate = thermal_tail(dimension, thermal, radius)
        assert sum(probabilities[dimension:]) <= float(certificate["occupation_tail_upper"])+1e-15


@pytest.mark.parametrize("alpha,beta", [(1.1+.4j, -3+2j), (-.7, 4+5j), (0, 5j)])
def test_pointer_radius_covers_whole_two_arm_path(alpha, beta):
    certificate = cutoff_certificate(
        32, thermal=F(1, 25), alpha_real=F(float(alpha.real)),
        alpha_imag=F(float(alpha.imag)), beta_real=F(float(beta.real)),
        beta_imag=F(float(beta.imag)))
    radius = float(certificate["pointer_tail_certificate"]["squared_displacement_upper"])
    for arm in (-1, 1):
        for position in np.linspace(0, 1, 101):
            assert abs(alpha+arm*position*beta/2)**2 <= radius+2e-14


@pytest.mark.parametrize("dimension,alpha,beta,nbar", [
    (8, .7+.2j, .5+1.2j, .04),
    (16, -.95, 2+1j, .125),
    (32, 1.05, 3+2j, .04),
    (64, .96341309, 4+4j, .12777123),
    (64, .76974579, 5+5j, 0.),
])
def test_finite_preparation_and_readout_fit_unrestricted_bound(dimension, alpha, beta, nbar):
    # Independent dense oscillator construction. The finite pointer's signed
    # mean equals Tr(rho_d exp(beta*a_d†-beta* a_d)); no infinite effect is
    # compressed and no Gaussian formula supplies the finite prediction.
    annihilator = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    displacement = expm(alpha*annihilator.T-np.conj(alpha)*annihilator)
    ratio = nbar/(1+nbar)
    weights = (1-ratio)*ratio**np.arange(dimension)
    weights[-1] = ratio**(dimension-1)
    rho = (displacement*weights)@displacement.conj().T
    characteristic = expm(beta*annihilator.T-np.conj(beta)*annihilator)
    observed = np.trace(rho@characteristic)
    reference = cmath.exp(-(nbar+.5)*abs(beta)**2
                         +beta*np.conj(alpha)-np.conj(beta)*alpha)
    certificate = cutoff_certificate(
        dimension, thermal=F(float(nbar)), alpha_real=F(float(np.real(alpha))),
        alpha_imag=F(float(np.imag(alpha))), beta_real=F(float(np.real(beta))),
        beta_imag=F(float(np.imag(beta))))
    probability_error = float(certificate["probability_error_upper"])
    assert abs(observed.real-reference.real)/2 <= probability_error+2e-14
    assert abs(observed.imag-reference.imag)/2 <= probability_error+2e-14
    if dimension == 64:
        assert probability_error < 4e-7


def test_zero_coordinate_has_exact_probability_identity():
    certificate = cutoff_certificate(4, thermal=F(1), alpha_real=F(4),
                                     beta_real=F(0), beta_imag=F(0))
    assert certificate["probability_error_upper"] == 0
    assert certificate["complete_instrument_trace_distance_upper"] > 0


def test_complete_retained_pointer_and_field_state_bound():
    alpha, beta, nbar, phase = .6, 1+.6j, .1, -.7

    def output(dimension):
        a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
        displaced = expm(alpha*(a.T-a))
        ratio = nbar/(1+nbar)
        weights = (1-ratio)*ratio**np.arange(dimension)
        weights[-1] = ratio**(dimension-1)
        rho = (displaced*weights)@displaced.conj().T
        eigenvalues, eigenvectors = np.linalg.eigh(-1j*(beta*a.T-beta.conjugate()*a))
        branches = []
        for function in (np.cos, np.sin):
            k = (eigenvectors*function((eigenvalues+phase)/2))@eigenvectors.conj().T
            branches.append(k@rho@k.conj().T)
        return branches

    finite, reference = output(8), output(64)
    distance = 0
    for low, high in zip(finite, reference):
        embedded = np.zeros_like(high)
        embedded[:8, :8] = low
        distance += sum(abs(np.linalg.eigvalsh(embedded-high)))/2
    arguments = dict(thermal=F(nbar), alpha_real=F(alpha),
                     beta_real=F(beta.real), beta_imag=F(beta.imag))
    low_bound = float(cutoff_certificate(8, **arguments)["complete_instrument_trace_distance_upper"])
    high_bound = float(cutoff_certificate(64, **arguments)["complete_instrument_trace_distance_upper"])
    assert 1e-6 < distance < low_bound < .2
    assert distance <= low_bound+high_bound+2e-14


def test_geometric_zero_displacement_bound_has_no_rounding_floor():
    for dimension in (8, 64):
        assert thermal_tail(dimension, F(0), F(0))["occupation_tail_upper"] == 0
        assert thermal_tail(dimension, F(1, 25), F(0))["occupation_tail_upper"] == F(1, 26)**dimension


@pytest.mark.parametrize("replacement", [
    {"dimension": True}, {"dimension": 1}, {"thermal": -.1},
    {"thermal": F(-1)}, {"alpha_real": .1}, {"beta_real": 2},
    {"beta_imag": float("nan")}, {"digits": 2}, {"digits": True},
])
def test_certificate_rejects_non_exact_or_invalid_inputs(replacement):
    arguments = dict(dimension=16, thermal=F(1, 25), alpha_real=F(1),
                     beta_real=F(1), beta_imag=F(2))
    arguments.update(replacement)
    with pytest.raises(ValueError):
        cutoff_certificate(**arguments)


@pytest.mark.parametrize("q,nbar", [(F(1), F(0)), (F(3), F(1)),
                                   (F(2), F(-1)), (1.5, F(0))])
def test_moment_requires_integrable_positive_weight(q, nbar):
    with pytest.raises(ValueError):
        thermal_moment(q, nbar, F(1))
