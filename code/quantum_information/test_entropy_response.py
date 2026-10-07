"""Independent scalar spectra, complete tangent bases and hostile inputs."""

import math

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm, logm

from quantum_information.entropy_response import (
    central_normalization_diagnostic, entropy_tangent, finite_entropy_balance,
    first_law_diagnostic, gibbs_entropy_response,
)

X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])
I = np.eye(2)


def mp_matrix(a):
    return mp.matrix([[mp.mpc(float(z.real), float(z.imag)) for z in row] for row in a])


def qubit_entropy(a):
    # Closed characteristic polynomial, independent of production eigensolves.
    radius = mp.sqrt((a[0, 0]-a[1, 1])**2+4*a[0, 1]*a[1, 0])
    values = [(a[0, 0]+a[1, 1]+radius)/2, (a[0, 0]+a[1, 1]-radius)/2]
    return mp.re(-sum(p*mp.log(p) for p in values if p))


@pytest.mark.parametrize("epsilon", [1e-4, 1e-8, 1e-12, 1e-16, 1e-100, 1e-150])
def test_entropy_loss_below_epsilon_is_not_a_zero_first_law_error(epsilon):
    rho = I/2 + epsilon*X
    with mp.workdps(330):
        expected = float(qubit_entropy(mp_matrix(rho))-mp.log(2))
    report = finite_entropy_balance(rho, I/2)
    assert report["modular_change"] == 0
    assert report["relative_entropy"] > 0
    assert report["entropy_change"] == pytest.approx(expected, rel=2e-14, abs=0)


@pytest.mark.parametrize("seed", range(8))
def test_noncommuting_finite_balance_matches_independent_entropy_spectra(seed):
    rng = np.random.default_rng(seed)
    a = rng.uniform(-.15, .15, 6)
    rho = I/2+a[0]*X+a[1]*Y+a[2]*Z
    sigma = I/2+a[3]*X+a[4]*Y+a[5]*Z
    with mp.workdps(90):
        expected = float(qubit_entropy(mp_matrix(rho))-qubit_entropy(mp_matrix(sigma)))
    report = finite_entropy_balance(rho, sigma)
    assert report["entropy_change"] == pytest.approx(expected, abs=2e-15)
    # Independent full matrix logarithm, including the sign of the remainder.
    linear = -np.trace(logm(sigma) @ (rho-sigma)).real
    assert report["modular_change"] == pytest.approx(linear, abs=2e-15)
    assert report["relative_entropy"] == pytest.approx(linear-expected, abs=2e-15)


@pytest.mark.parametrize("direction", [X, Y, Z, X+2*Y-3*Z])
def test_tangent_matches_high_precision_derivative_of_characteristic_polynomial(direction):
    sigma = np.array([[.625, .0625+.125j], [.0625-.125j, .375]])
    with mp.workdps(90):
        a, d = mp_matrix(sigma), mp_matrix(direction)
        expected = float(mp.diff(lambda x: qubit_entropy(a+x*d), 0))
    assert entropy_tangent(sigma, direction) == pytest.approx(expected, abs=3e-15)


@pytest.mark.parametrize("dimension", [2, 3, 4])
def test_complete_diagnostic_equals_sum_over_an_independent_tangent_basis(dimension):
    rng = np.random.default_rng(dimension)
    m = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
    sigma = m @ m.conj().T + np.eye(dimension)
    sigma /= np.trace(sigma)
    k = m+m.conj().T  # Independently supplied wrong generator.
    report = first_law_diagnostic(sigma, k)
    basis = []
    for i in range(dimension):
        for j in range(i):
            e = np.zeros_like(m)
            e[i, j], e[j, i] = 1/math.sqrt(2), 1/math.sqrt(2)
            basis.append(e)
            e = np.zeros_like(m)
            e[i, j], e[j, i] = 1j/math.sqrt(2), -1j/math.sqrt(2)
            basis.append(e)
    for n in range(1, dimension):
        diagonal = [1.]*n+[-float(n)]+[0.]*(dimension-n-1)
        basis.append(np.diag(diagonal)/math.sqrt(n*(n+1)))
    log_sigma = logm(sigma)
    errors = [np.trace((k+log_sigma) @ b).real for b in basis]
    assert report["all_tangent_defect"] == pytest.approx(math.hypot(*errors), rel=3e-15)
    witness = report["witness"]
    attained = np.trace(k @ witness).real-entropy_tangent(sigma, witness)
    assert attained == pytest.approx(report["all_tangent_defect"], rel=3e-15)
    assert report["witness_response"] == pytest.approx(attained, rel=3e-15)
    assert report["witness_trace_defect"] < 1e-15
    assert report["witness_norm_defect"] < 1e-15


@pytest.mark.parametrize("amplitude", [1., 1e-100, 1e-250, 1e-310])
@pytest.mark.parametrize("origin", [0., 1e20, -1e100])
def test_wrong_coherent_generator_survives_tiny_units_and_huge_scalar_origin(amplitude, origin):
    report = first_law_diagnostic(I/2, origin*I+amplitude*Y)
    assert report["all_tangent_defect"] == pytest.approx(math.sqrt(2)*amplitude, rel=2e-13, abs=0)
    np.testing.assert_allclose(report["witness"], Y/math.sqrt(2), atol=3e-14, rtol=0)


def test_equilibrium_and_scalar_generator_have_no_false_witness():
    report = first_law_diagnostic(I/2, 1e100*I)
    assert report["all_tangent_defect"] == 0
    assert not report["witness"].any()
    assert finite_entropy_balance(I/2, I/2)["entropy_change"] == 0


def test_complete_central_test_and_explicit_probability_transfer():
    report = central_normalization_diagnostic([2, 3, 5], [0., 0., 0.])
    assert report["all_sector_transfer_defect"] == pytest.approx(math.log(5/2), abs=1e-15)
    np.testing.assert_array_equal(report["witness"], [1., 0., -1.])
    # First variations cannot fix a common origin of the central observable.
    common = central_normalization_diagnostic([2, 3, 5], np.log([2., 3., 5.])+7)
    assert common["all_sector_transfer_defect"] < 2e-15
    assert central_normalization_diagnostic([4], [-12.])["all_sector_transfer_defect"] == 0


@pytest.mark.parametrize("lam", [-4., -.2, 0., 1e-20, .7, 10., 300.])
def test_gibbs_derivatives_match_separate_high_precision_partition_function(lam):
    energies = [0., 1., 2.]
    with mp.workdps(160):
        def moments(x):
            weights = [mp.exp(-x*e) for e in energies]
            p = [w/sum(weights) for w in weights]
            return sum(v*e for v, e in zip(p, energies)), -sum(v*mp.log(v) for v in p)
        dt = float(mp.diff(lambda x: moments(x)[0], mp.mpf(lam)))
        ds = float(mp.diff(lambda x: moments(x)[1], mp.mpf(lam)))
    report = gibbs_entropy_response(np.diag(energies), lam)
    assert report["dt_dlambda"] == pytest.approx(dt, rel=5e-14, abs=0)
    if lam:
        assert report["ds_dlambda"] == pytest.approx(ds, rel=5e-14, abs=0)
        assert report["ds_dt"] == pytest.approx(lam, rel=3e-16, abs=0)
    else:
        assert report["ds_dlambda"] == report["ds_dt"] == 0


@pytest.mark.parametrize("scale", [1e-100, 1., 1e100])
def test_gibbs_response_respects_energy_units_and_origin(scale):
    beta = .3/scale
    report = gibbs_entropy_response(1e200*I+scale*Y, beta)
    np.testing.assert_allclose(report["state"], (I-math.tanh(.3)*Y)/2, atol=2e-15, rtol=0)
    assert report["energy_variance"] == pytest.approx(scale**2/math.cosh(.3)**2, rel=3e-15, abs=0)
    assert report["ds_dlambda"] == pytest.approx(-.3*scale/math.cosh(.3)**2, rel=3e-15, abs=0)


def test_finite_source_can_be_pure_but_reference_must_be_faithful():
    pure = np.diag([1., 0.])
    assert finite_entropy_balance(pure, I/2)["entropy_change"] == pytest.approx(-math.log(2))
    with pytest.raises(ValueError, match="faithful"):
        finite_entropy_balance(I/2, pure)
    with pytest.raises(ValueError, match="support"):
        finite_entropy_balance(np.diag([1., -1e-14]), I/2)


def test_trace_roundoff_is_reported_not_silently_normalized():
    rho = np.diag([.5+2**-45, .5])
    with mp.workdps(90):
        expected = float(qubit_entropy(mp_matrix(rho))-mp.log(2))
    report = finite_entropy_balance(rho, I/2)
    assert report["trace_difference"] == 2**-45
    assert report["bregman_remainder"] > 0
    assert report["entropy_change"] == pytest.approx(expected, rel=2e-15, abs=0)


@pytest.mark.parametrize("direction", [1e-200*I, np.array([[0., 1e-200], [0., 0.]]),
                                        np.zeros((3, 3)), [[False, 1], [1, 0]]])
def test_invalid_tangent_is_rejected_at_its_own_scale(direction):
    with pytest.raises(ValueError):
        entropy_tangent(I/2, direction)


@pytest.mark.parametrize("bad", [np.nan*I, [[True, 0], [0, False]], [[".5", "0"], ["0", ".5"]],
                                np.ma.array(I/2, mask=[[False, True], [False, False]])])
@pytest.mark.parametrize("function", [lambda a: entropy_tangent(a, X),
                                     lambda a: finite_entropy_balance(a, I/2),
                                     lambda a: first_law_diagnostic(I/2, a),
                                     lambda a: gibbs_entropy_response(a, .5)])
def test_invalid_numeric_inputs_are_not_converted_into_valid_evidence(bad, function):
    with pytest.raises(ValueError):
        function(bad)


@pytest.mark.parametrize("call", [lambda: finite_entropy_balance(I/2+1e-200*X, I/2),
                                  lambda: gibbs_entropy_response(I, .4),
                                  lambda: gibbs_entropy_response(1e-200*X, 1.),
                                  lambda: gibbs_entropy_response(Z, 1000.),
                                  lambda: central_normalization_diagnostic([True], [0.]),
                                  lambda: central_normalization_diagnostic([2, 3], [0.])])
def test_undefined_or_unresolved_responses_raise_instead_of_zero_or_nan(call):
    with pytest.raises(ValueError):
        call()


def test_gibbs_state_matches_independent_matrix_exponential():
    observable = .25*X+.375*Y+.5*Z
    thermal = expm(-.7*observable)
    expected = thermal/np.trace(thermal)
    np.testing.assert_allclose(gibbs_entropy_response(observable, .7)["state"], expected, atol=1e-15, rtol=0)
