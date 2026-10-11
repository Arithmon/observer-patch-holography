"""Original-entry controls for finite operator actions and residual norms."""

from fractions import Fraction

import numpy as np
import pytest
import sympy as sp

from quantum_information.algebras import FiniteAlgebra, matrix_units
from quantum_information.expectations import (
    apply_map, repair_generator, state_preserving_expectation,
)


def _exact(values):
    return sp.Matrix([
        [sp.Rational(float(z.real)) + sp.I*sp.Rational(float(z.imag))
         for z in row] for row in np.asarray(values)
    ])


def _rounded(values, shape):
    return np.array([complex(z) for z in values], dtype=complex).reshape(shape)


@pytest.mark.parametrize("consumer", ("projection", "expectation", "repair"))
@pytest.mark.parametrize("phase", (1., 1j))
def test_cancelled_trace_retains_ordinary_size_observable(consumer, phase):
    identity = np.eye(4)
    algebra = FiniteAlgebra([identity])
    values = [2.**56, -2.**56, 1., 0.]
    x = phase*np.diag(values)
    # The trace is exactly one in the original inputs, independently of
    # the summation order used by either public matrix action.
    average = sum((Fraction(v) for v in values), Fraction())/4
    assert average == Fraction(1, 4)
    expected = phase*float(average)*identity
    if consumer == "projection":
        actual = algebra.project(x)
    else:
        channel = state_preserving_expectation(algebra, identity/4)
        if consumer == "expectation":
            actual = apply_map(channel, x)
        else:
            generator, receipt = repair_generator(
                [channel], [algebra], identity/4, [1.])
            assert receipt["intersection_dimension"] == 1
            assert receipt["gap"] == pytest.approx(1.)
            actual = apply_map(generator, x)
            expected = _rounded(_exact(expected)-_exact(x), x.shape)
            assert actual[2, 2] == -.75*phase
    assert actual[3, 3] == .25*phase
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("imaginary", (False, True))
def test_scalar_projection_needs_no_representable_intermediate_coefficient(imaginary):
    maximum = np.finfo(float).max
    x = np.eye(4, dtype=complex)*maximum
    if imaginary:
        x.imag = np.eye(4)*maximum
    algebra = FiniteAlgebra([np.eye(4)])
    # Q=I/2 gives a coefficient 2*maximum, but Q times that coefficient
    # is exactly the finite supplied scalar matrix, in each real component.
    with np.errstate(over="ignore", invalid="ignore"):
        actual = algebra.project(x)
    assert np.all(np.isfinite(actual))
    np.testing.assert_array_equal(actual, x)


@pytest.mark.parametrize("consumer", ("projection", "expectation"))
def test_subnormal_products_can_sum_to_a_representable_fixed_observable(consumer):
    tiny = np.nextafter(0., 1.)
    x = np.eye(4, dtype=complex)*tiny
    x.imag = np.eye(4)*tiny
    algebra = FiniteAlgebra([np.eye(4)])
    # Every product of a diagonal entry with 1/2 or 1/4 rounds to zero,
    # but the exact scalar projection fixes both components of this matrix.
    with np.errstate(under="ignore"):
        if consumer == "projection":
            actual = algebra.project(x)
        else:
            channel = state_preserving_expectation(algebra, np.eye(4)/4)
            actual = apply_map(channel, x)
    np.testing.assert_array_equal(actual, x)


@pytest.mark.parametrize("exponent", (-1074, -600, 600))
def test_distance_retains_representable_single_component_norm(exponent):
    algebra = FiniteAlgebra([np.eye(4)])
    x = np.zeros((4, 4))
    x[0, 1] = 2.**exponent
    # E01 is exactly orthogonal to the scalar algebra, so its HS distance
    # is its single nonzero entry; no norm implementation supplies the oracle.
    with np.errstate(over="ignore", under="ignore"):
        actual = algebra.distance(x)
    assert actual == 2.**exponent


def test_small_residual_survives_a_large_conserved_component():
    algebra = FiniteAlgebra([np.eye(4)])
    x = 2.**500*np.eye(4)
    x[0, 1] = 2.**-600
    assert algebra.distance(x) == 2.**-600


@pytest.mark.parametrize("full", (False, True))
@pytest.mark.parametrize("consumer", ("projection", "expectation"))
def test_mixed_scale_fixed_observables_keep_each_original_coordinate(full, consumer):
    basis = matrix_units(4) if full else [np.diag(row) for row in np.eye(4)]
    algebra = FiniteAlgebra(basis)
    x = np.diag([np.finfo(float).max, np.nextafter(0., 1.), 2.**-600, 0.])
    if full:
        x = x.astype(complex)
        x[0, 1], x[1, 2] = 2.**600*(1.+1j), 1j*2.**-600
    if consumer == "projection":
        actual = algebra.project(x)
    else:
        channel = state_preserving_expectation(algebra, np.eye(4)/4)
        actual = apply_map(channel, x)
    # A common rescaling of the whole input can silently erase the last
    # two nonzero entries even though this action is simply the identity.
    np.testing.assert_array_equal(actual, x)
    assert algebra.distance(x) == 0.


@pytest.mark.parametrize("dual", (False, True))
def test_complex_map_action_uses_all_original_terms_and_correct_adjoint(dual):
    target = np.array([[1., -1., 1., 0.], [1j, -1j, 0., 1.],
                       [0., 0., 1.+1j, 0.], [1., -1., 0., 1j]])
    supplied = target.conj().T if dual else target
    x = np.array([[2.**56, 2.**56], [1., 2.+3j]])
    exact_map = _exact(supplied)
    if dual:
        exact_map = exact_map.H
    expected = _rounded(exact_map*_exact(x).reshape(4, 1), x.shape)
    assert expected[0, 0] == 1.
    np.testing.assert_array_equal(apply_map(supplied, x, dual=dual), expected)


def test_map_products_can_overflow_before_a_finite_cancelled_sum():
    maximum = np.finfo(float).max
    supplied = np.array([[maximum, maximum, -maximum, -maximum],
                         [maximum, -maximum, 1., 0.],
                         [0., 0., .5, 0.], [0., 0., 0., -.5]])
    x = np.full((2, 2), 2.)
    expected = _rounded(_exact(supplied)*_exact(x).reshape(4, 1), x.shape)
    np.testing.assert_array_equal(expected, [[0., 2.], [1., -1.]])
    with np.errstate(over="ignore", invalid="ignore"):
        actual = apply_map(supplied, x)
    np.testing.assert_array_equal(actual, expected)


def test_nontrivial_projection_uses_the_original_operator_span():
    n = np.kron(np.array([[0., 1.-2j], [1.+2j, 0.]]), np.eye(2))
    supplied = [np.eye(4), n]
    algebra = FiniteAlgebra(supplied)
    x = np.diag([2.**56, -2.**56, 1., 0.]).astype(complex)
    x[0, 2] = 1.+3j
    columns = _exact(np.column_stack([b.reshape(-1) for b in supplied]))
    # Solve the original rational Gram system, before any numerical chart
    # or operator normalization has changed the supplied entries.
    projector = columns*(columns.H*columns).inv()*columns.H
    expected = _rounded(projector*_exact(x).reshape(16, 1), x.shape)
    actual = algebra.project(x)
    assert expected[3, 3] == .25
    np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=0.)


@pytest.mark.parametrize("phase", (1., 1j))
def test_nonzero_final_map_component_cannot_silently_underflow(phase):
    x = np.eye(2, dtype=complex)*np.nextafter(0., 1.)*phase
    # Each nonzero exact output is half a minimum subnormal, which rounds
    # to zero. This differs from several tiny products summing to a finite
    # representable result and requires the API's explicit refusal.
    with np.errstate(under="ignore"):
        with pytest.raises(ValueError, match="underflow|represent|range"):
            apply_map(.5*np.eye(4), x)


@pytest.mark.parametrize("consumer", ("map", "projection", "distance"))
def test_genuinely_unrepresentable_final_results_raise(consumer):
    maximum = np.finfo(float).max
    if consumer == "projection":
        supplied = [np.eye(2), np.array([[3., 4.], [4., -3.]])]
        algebra = FiniteAlgebra(supplied)
        x = np.full((2, 2), maximum)
        columns = _exact(np.column_stack([b.reshape(-1) for b in supplied]))
        projector = columns*(columns.H*columns).inv()*columns.H
        exact_result = projector*_exact(x).reshape(4, 1)
        assert exact_result[0] > sp.Rational(maximum)
    with np.errstate(over="ignore", invalid="ignore"):
        with pytest.raises(ValueError, match="finite|range|overflow|represent"):
            if consumer == "map":
                apply_map(2*np.eye(4), maximum*np.eye(2))
            elif consumer == "projection":
                algebra.project(x)
            else:
                x = np.zeros((4, 4))
                x[0, 1] = x[1, 0] = maximum
                FiniteAlgebra([np.eye(4)]).distance(x)
