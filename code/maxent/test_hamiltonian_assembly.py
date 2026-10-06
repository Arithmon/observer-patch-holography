"""Independent controls for the Hamiltonian supplied to finite A3 numerics."""

import itertools
from fractions import Fraction

import numpy as np
import pytest
from scipy.linalg import expm

from maxent.information_projection import (
    constrained_hamiltonian, duhamel_covariance, gibbs_state,
    hamiltonian_assembly_diagnostics,
)


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])
HEISENBERG = sum(np.kron(a, a) for a in (X, Y, Z))


@pytest.mark.parametrize('order', tuple(itertools.permutations(range(3))))
def test_bookkeeping_cannot_erase_an_entangled_thermal_state(order):
    coefficients = [1e20, 1., -1e20]
    multipliers = [coefficients[j] for j in order]
    # This exact binary64 sum is one in every order.
    assert sum(map(Fraction, multipliers)) == 1
    actual_h = constrained_hamiltonian([HEISENBERG]*3, multipliers)
    np.testing.assert_array_equal(actual_h, HEISENBERG)
    rho, log_z = gibbs_state([HEISENBERG]*3, multipliers)
    # H=I-4|singlet><singlet|: diagonalize the known projectors, not H.
    singlet = np.array([0., 1., -1., 0.])/np.sqrt(2)
    projector = np.outer(singlet, singlet)
    probability = 1/(1+3*np.exp(-4))
    expected = probability*projector+(1-probability)*(np.eye(4)-projector)/3
    np.testing.assert_allclose(rho, expected, atol=2e-15, rtol=0)
    assert log_z == pytest.approx(3+np.log1p(3*np.exp(-4)), abs=2e-15)
    partial_transpose = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    assert np.linalg.eigvalsh(partial_transpose)[0] == pytest.approx(.5-probability, abs=2e-15)


@pytest.mark.parametrize('observable', (X, Y, Z, HEISENBERG))
def test_products_are_accumulated_before_their_roundoff_can_cancel(observable):
    # (1+2^-27)(1-2^-27)-1 = -2^-54 exactly. Rounding the product
    # first gives zero even with a compensated sum of rounded products.
    a, b = 1+2.**-27, 1-2.**-27
    expected = -2.**-54*observable
    actual = constrained_hamiltonian([b*observable, observable], [a, -1.])
    np.testing.assert_array_equal(actual, expected)


def test_scalar_partition_energy_survives_product_cancellation():
    a, b = 1+2.**-27, 1-2.**-27
    # Three terms have exact scalar energy -1 after the two large products
    # cancel. This affects log Z even though it cannot change the state.
    rho, log_z = gibbs_state([b*I, I, X], [a*2.**54, -2.**54, .25])
    np.testing.assert_allclose(rho, (I-np.tanh(.25)*X)/2, atol=2e-15, rtol=0)
    assert log_z == pytest.approx(1+np.log(2*np.cosh(.25)), abs=2e-15)


def test_large_cancelled_terms_do_not_hide_a_noncommuting_remainder():
    small = .3*X+.2*Y+.4*Z
    rho, _ = gibbs_state([X, Y, small, X, Y], [1e20, -1e20, 1., -1e20, 1e20])
    exponential = expm(-small)
    np.testing.assert_allclose(rho, exponential/np.trace(exponential), atol=2e-15, rtol=0)


def test_product_roundoff_can_erase_an_order_one_thermal_interaction():
    a, b = 1+2.**-27, 1-2.**-27
    # -2^54 ((1+2^-27)(1-2^-27)-1) H = H exactly.
    rho, _ = gibbs_state([b*HEISENBERG, HEISENBERG], [-a*2.**54, 2.**54])
    exponential = expm(-HEISENBERG)
    np.testing.assert_allclose(rho, exponential/np.trace(exponential), atol=2e-15, rtol=0)


@pytest.mark.parametrize('coefficient', (1e20, 1e200))
def test_intermediate_products_can_cancel_outside_binary64_range(coefficient):
    # Each input is finite; even products outside binary64 can cancel to X.
    actual = constrained_hamiltonian([1e200*X, X, 1e200*X],
                                     [coefficient, 1., -coefficient])
    np.testing.assert_array_equal(actual, X)
    rho, _ = gibbs_state([1e200*X, X, 1e200*X], [coefficient, 1., -coefficient])
    np.testing.assert_allclose(rho, (I-np.tanh(1.)*X)/2, atol=2e-15, rtol=0)


def test_subnormal_products_accumulate_before_conversion():
    tiny = np.nextafter(0., 1.)
    actual = constrained_hamiltonian([tiny*X, tiny*X], [.5, .5])
    np.testing.assert_array_equal(actual, tiny*X)
    with pytest.raises(ValueError, match='underflow'):
        constrained_hamiltonian([tiny*X], [.5])
    with pytest.raises(ValueError, match='range'):
        constrained_hamiltonian([1e200*X], [1e200])


def test_reconstruction_of_the_full_matrix_does_not_cancel_a_diagonal_entry():
    # Centering and adding the scalar back in binary64 would erase the 2.
    h = np.diag([1e20, 2.])
    np.testing.assert_array_equal(constrained_hamiltonian([h], [1.]), h)


@pytest.mark.parametrize('observable', (Z, HEISENBERG))
def test_scalar_origin_is_removed_before_diagonal_variation_is_rounded(observable):
    rho, _ = gibbs_state([np.eye(len(observable)), observable], [1e20, 1.])
    exponential = expm(-observable)
    np.testing.assert_allclose(rho, exponential/np.trace(exponential), atol=2e-15, rtol=0)


def test_duhamel_response_uses_the_retained_interaction():
    probability = 1/(1+3*np.exp(-4))
    expected_variance = 16*probability*(1-probability)
    covariance = duhamel_covariance([HEISENBERG]*3, [1e20, 1., -1e20])
    np.testing.assert_allclose(covariance, np.full((3, 3), expected_variance),
                               atol=2e-14, rtol=0)


@pytest.mark.parametrize('bad', ([[True, 0], [0, 1]],
                               np.ma.array(X, mask=[[False, True], [False, False]]),
                               np.array([[np.nan, 0], [0, 1.]]),
                               np.array([[0., 1.], [0., 0.]])))
def test_assembly_diagnostics_cannot_certify_invalid_inputs(bad):
    with pytest.raises(ValueError):
        hamiltonian_assembly_diagnostics([bad], [1.])


def test_certificate_does_not_silently_replace_a_nonhermitian_input():
    almost = X.copy(); almost[0, 1] += 2.**-50
    # Gibbs numerics tolerate Hermiticity roundoff, but the exact-input
    # certificate must not certify only a symmetrized replacement.
    gibbs_state([almost], [1.])
    with pytest.raises(ValueError, match='exact Hermitian'):
        hamiltonian_assembly_diagnostics([almost], [1.])


def test_rounding_certificate_matches_independent_exact_entry_errors():
    import sympy as sp
    # The scalar 1/3 and centered thirds force nonzero conversion errors.
    h = np.diag([1., 0., 0.]).astype(complex)
    h[0, 1], h[1, 0] = .3j, -.3j
    report = hamiltonian_assembly_diagnostics([h], [.7])
    exact = sp.Rational(.7)*sp.Matrix(h).applyfunc(
        lambda x: sp.Rational(float(sp.re(x)))+sp.I*sp.Rational(float(sp.im(x))))
    mean = sp.trace(exact)/3
    rounded = sp.Matrix(report.matrix).applyfunc(
        lambda x: sp.Rational(float(sp.re(x)))+sp.I*sp.Rational(float(sp.im(x))))
    error = rounded-(exact-mean*sp.eye(3))
    row_bound = max(sum(abs(sp.re(error[i,j]))+abs(sp.im(error[i,j]))
                        for j in range(3)) for i in range(3))
    scalar_bound = abs(sp.Rational(report.scalar)-mean)
    assert row_bound > 0 and scalar_bound > 0
    for value, bound in ((row_bound, report.operator_rounding_bound),
                         (scalar_bound, report.scalar_rounding_bound),
                         (row_bound+scalar_bound, report.ideal_log_partition_error_bound)):
        assert sp.Rational(bound) >= value
        assert sp.Rational(float(np.nextafter(bound, 0.))) < value
    assert report.ideal_gibbs_trace_distance_bound == report.operator_rounding_bound


def test_exact_assembly_certificate_distinguishes_zero_from_tiny_error():
    report = hamiltonian_assembly_diagnostics([HEISENBERG]*3, [1e20, 1., -1e20])
    assert report.operator_rounding_bound == report.scalar_rounding_bound == 0
    # An error below the smallest float must get an upward bound, not zero.
    from maxent.hamiltonian_assembly import _upper_float
    exact = Fraction(1, 2**2000)
    assert Fraction(_upper_float(exact)) >= exact > 0


def test_ideal_gibbs_bounds_cover_independent_high_precision_functions():
    import mpmath as mp
    ctx = mp.mp.clone(); ctx.dps = 90
    h = np.diag([1., 0., 0.]).astype(complex)
    h[0,1], h[1,0] = .3j, -.3j
    report = hamiltonian_assembly_diagnostics([h], [.7])
    def exact_float_matrix(matrix):
        return ctx.matrix([[ctx.mpc(float(x.real), float(x.imag)) for x in row]
                           for row in matrix])
    source = ctx.mpf(.7)*exact_float_matrix(h)
    assembled = exact_float_matrix(report.matrix)+ctx.mpf(report.scalar)*ctx.eye(3)
    source_exp, assembled_exp = ctx.expm(-source), ctx.expm(-assembled)
    z_source, z_assembled = (sum(e[i,i] for i in range(3)) for e in (source_exp, assembled_exp))
    difference = source_exp/z_source-assembled_exp/z_assembled
    distance = sum(abs(x) for x in ctx.eighe(difference, eigvals_only=True))
    assert 0 < distance <= ctx.mpf(report.ideal_gibbs_trace_distance_bound)
    assert abs(ctx.log(z_source)-ctx.log(z_assembled)) <= ctx.mpf(report.ideal_log_partition_error_bound)
