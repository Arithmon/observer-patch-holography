"""Original-entry controls for numerical charts of exact operator spans."""

from fractions import Fraction

import numpy as np
import pytest
import sympy as sp

from geometry.null_net_standardness import separating_modulus
from quantum_information.algebras import (
    FiniteAlgebra, algebra_intersection, matrix_units,
    orthonormal_operator_basis, separation_modulus,
)
from quantum_information.expectations import (
    expectation_diagnostics, repair_generator, state_preserving_expectation,
)


def _exact(matrix):
    """Convert the supplied entries themselves, before any normalization."""
    return sp.Matrix([
        [sp.Rational(float(x.real)) + sp.I*sp.Rational(float(x.imag))
         for x in row] for row in np.asarray(matrix)
    ])


def _modulus_from_original_entries(basis, omega):
    # Solve det(H-lambda G)=0 with exact rational Gram components. This
    # neither uses the production chart nor subtracts rounded near-equals.
    matrices = [_exact(a) for a in basis]
    vector = _exact(np.asarray(omega).reshape(-1, 1))
    gram = sp.Matrix([[(a.H*b).trace() for b in matrices]
                      for a in matrices]).applyfunc(sp.expand)
    action = sp.Matrix.hstack(*(a*vector for a in matrices))
    h = (action.H*action).applyfunc(sp.expand)
    assert gram.det() > 0
    coefficient = sp.expand(gram[0, 0]*h[1, 1] + gram[1, 1]*h[0, 0]
                            - gram[0, 1]*h[1, 0] - gram[1, 0]*h[0, 1])
    if h.det() == 0:
        return 0.
    squared = 2*h.det()/(coefficient + sp.sqrt(
        coefficient**2 - 4*gram.det()*h.det()))
    return float(sp.N(sp.sqrt(squared), 80))


def _projection_from_original_entries(basis):
    columns = sp.Matrix.hstack(*[
        _exact(a).reshape(np.asarray(a).size, 1) for a in basis
    ])
    projector = columns*(columns.H*columns).inv()*columns.H
    return np.array(projector.evalf(40)).astype(complex)


@pytest.mark.parametrize("power", (20, 40))
@pytest.mark.parametrize("presentation", ("original", "reversed", "complex_shear"))
@pytest.mark.parametrize("slot", (0, 1))
def test_near_dependent_chart_preserves_zero_and_positive_separation(power, presentation, slot):
    delta = 2.**-power
    a = np.array([[1., 2.], [3., 4.]])
    b = np.array([[1., 2.+delta], [3., 4.-delta]])
    basis = [a, b]
    if presentation == "reversed":
        basis.reverse()
    elif presentation == "complex_shear":
        basis = [a+2j*b, b]
    omega = np.eye(2)[slot]
    expected = _modulus_from_original_entries(basis, omega)
    # The nonzero supplied B-A annihilates e0 exactly. The same span acting
    # on e1 has modulus sqrt(9/14), independently excluding a zero floor.
    assert expected == (0. if slot == 0 else pytest.approx(np.sqrt(9/14)))
    for entry in (separation_modulus, separating_modulus):
        actual = entry(basis, omega)
        if expected == 0:
            assert actual == 0.
        else:
            assert actual == pytest.approx(expected, rel=3e-13, abs=0)


@pytest.mark.parametrize("exponent", (-1000, 0, 1000))
@pytest.mark.parametrize("phase", (1., 1j))
def test_rescaling_keeps_original_input_positive_generalized_gram_value(exponent, phase):
    a = np.array([[1., 2.], [3., 4.]])
    v = np.array([[4., 1.], [-2., 3.]])
    basis = [phase*2.**exponent*a, 2.**(-exponent)*(a+2.**-42*v)]
    omega = np.array([1., 0.])
    expected = _modulus_from_original_entries(basis, omega)
    assert expected == pytest.approx(np.sqrt((79-5*np.sqrt(85))/126))
    assert separation_modulus(basis, omega) == pytest.approx(expected, rel=3e-13, abs=0)


def test_complex_span_and_reference_are_not_replaced_by_their_conjugates():
    unitary = np.diag([1., 1j])
    a = np.array([[1., 2.], [3., 4.]])
    v = np.array([[4., 1.], [-2., 3.]])
    basis = [unitary@b@unitary.conj().T for b in (a, a+2.**-40*v)]
    omega = np.array([.6, .8j])
    expected = _modulus_from_original_entries(basis, omega)
    conjugate = _modulus_from_original_entries([b.conj() for b in basis], omega)
    assert abs(expected-conjugate) > .1
    assert separation_modulus(basis, omega) == pytest.approx(expected, rel=3e-13, abs=0)


@pytest.mark.parametrize("eta", (2.**-10, 2.**-20))
def test_nearby_positive_action_is_not_replaced_by_exact_kernel(eta):
    delta = 2.**-20
    basis = [np.array([[1., 2.], [3., 4.]]),
             np.array([[1.+delta*eta, 2.+delta], [3., 4.-delta]])]
    omega = np.array([1., 0.])
    expected = _modulus_from_original_entries(basis, omega)
    assert expected > 0
    # A well-conditioned binary64 chart still has ordinary action-product
    # roundoff. Its absolute allowance is far below the old 6e-10 error.
    assert separation_modulus(basis, omega) == pytest.approx(expected, rel=3e-12, abs=3e-15)


def test_unresolved_nonorthogonal_positive_action_is_not_certified_as_zero():
    epsilon = 2.**-100
    basis = [np.array([[1., 1.], [0., 1.]]),
             np.array([[1., 2.], [epsilon, 1.]])]
    omega = np.array([1., 0.])
    assert _modulus_from_original_entries(basis, omega) > 0
    with pytest.raises(ValueError, match="positive.*unresolved"):
        separation_modulus(basis, omega)


@pytest.mark.parametrize("epsilon", (0., 2.**-30, 2.**-40, 2.**-48))
@pytest.mark.parametrize("phase", (1., 1+2j))
def test_whole_action_cancellation_uses_original_source_scale(epsilon, phase):
    a = phase*np.array([[2., -1., -2., 1.], [7., -2., -4., -1.],
                        [6., -7., -1., 2.], [-3., -8., 4., 7.]])
    a[0, 0] += epsilon
    omega = np.full(4, .5)
    original, vector = _exact(a), _exact(omega.reshape(-1, 1))
    action = original*vector
    squared = sp.expand((action.H*action)[0])/sp.expand((original.H*original).trace())
    expected = float(sp.N(sp.sqrt(squared), 80))
    actual = separation_modulus([a], omega)
    if epsilon == 0:
        assert actual == 0.
    else:
        assert actual == pytest.approx(expected, rel=3e-13, abs=0)


@pytest.mark.parametrize("exponent", (-1000, -1074))
def test_representable_tiny_positive_action_is_not_a_zero_or_refusal(exponent):
    tiny = 2.**exponent
    basis = [np.diag([1., tiny])]
    # The exact value is tiny/sqrt(1+tiny**2), which rounds to tiny.
    assert separation_modulus(basis, np.array([0., 1.])) == tiny


def _singleton_squared_from_real_entries(matrix, omega):
    rows = [[Fraction(float(x)) for x in row] for row in matrix]
    vector = [Fraction(float(x)) for x in omega]
    image = [sum(x*y for x, y in zip(row, vector)) for row in rows]
    return sum(x*x for x in image)/sum(x*x for row in rows for x in row)


@pytest.mark.parametrize("lower", (0, 1, 2))
@pytest.mark.parametrize("side", ("below", "above"))
def test_singleton_rounding_near_subnormal_midpoints(lower, side):
    tiny = np.nextafter(0., 1.)
    t = Fraction(float(tiny))
    unit = 1. if side == "below" else 1.-2.**-52
    a = np.zeros((4, 4))
    a[0] = [unit, -unit, 1., -1.]
    a[1, 0] = (4*lower+2)*tiny
    omega = np.full(4, .5)  # exactly normalized, including in rational arithmetic
    squared = _singleton_squared_from_real_entries(a, omega)
    midpoint = (Fraction(lower)+Fraction(1, 2))*t
    # Locate the exact supplied-input quotient between adjacent floats and
    # on the specified side of their midpoint, without evaluating its root.
    assert (lower*t)**2 < squared < ((lower+1)*t)**2
    assert (squared < midpoint**2) == (side == "below")
    assert squared != midpoint**2
    expected = (lower if side == "below" else lower+1)*tiny
    if expected == 0:
        with pytest.raises(ValueError, match="underflow"):
            separation_modulus([a], omega)
    else:
        assert separation_modulus([a], omega) == expected


@pytest.mark.parametrize("multiplier,units", ((1, 0), (3, 2), (5, 2)))
def test_exact_subnormal_midpoints_obey_ties_to_even(multiplier, units):
    tiny = np.nextafter(0., 1.)
    a = np.zeros((4, 4))
    a[0, 0] = a[1, 1] = a[1, 2] = a[2, 2] = 1.
    # Its norm differs from one by far less than the API's declared 1e-12
    # allowance. The exact oracle keeps the supplied tiny component intact.
    omega = np.array([multiplier*tiny, 0., 0., 1.])
    squared = _singleton_squared_from_real_entries(a, omega)
    assert squared == (multiplier*Fraction(float(tiny))/2)**2
    if units == 0:
        with pytest.raises(ValueError, match="underflow"):
            separation_modulus([a], omega)
    else:
        assert separation_modulus([a], omega) == units*tiny


@pytest.mark.parametrize("exponent", (-1000, -1074))
def test_multiple_exact_orthogonal_tiny_actions_remain_positive(exponent):
    tiny = 2.**exponent
    a = np.diag([tiny, 1., 0.])
    b = np.array([[0., 0., 0.], [2*tiny, 0., 0.], [0., 0., 1.]])
    # Disjoint matrix entries give HS orthogonality. Their e0 images are
    # tiny*e0 and 2*tiny*e1, so the smaller squared modulus is tiny^2/(1+tiny^2).
    # Its exact root rounds to tiny at both specified scales.
    assert separation_modulus([a, b], np.array([1., 0., 0.])) == tiny


@pytest.mark.parametrize("tiny", (0., 2.**-1000))
def test_nested_python_lists_reach_exact_small_action_resolution(tiny):
    assert separation_modulus([[[1., 0.], [0., tiny]]], [0., 1.]) == tiny


def test_unrepresentable_operator_direction_requires_explicit_refusal():
    basis = [np.diag([2.**1000, 2.**-1000])]
    # Its second normalized entry is nonzero but smaller than binary64 can
    # represent. A chart dropping it would silently alter the supplied span.
    with pytest.raises(ValueError, match="underflow|unresolved"):
        orthonormal_operator_basis(basis)


def test_lost_component_cannot_hide_in_an_occupied_coordinate_row():
    basis = [np.array([[2.**-1000, 2.**1000], [0., 0.]]), np.eye(2)]
    # The first normalized component underflows, but the identity keeps its
    # coordinate row occupied. A row-support guard alone would silently return
    # span{I,E01}, losing the supplied nonzero diagonal imbalance.
    with pytest.raises(ValueError, match="underflow|unresolved"):
        orthonormal_operator_basis(basis)


@pytest.mark.parametrize("exponent", (-54, -1000, -1074))
@pytest.mark.parametrize("count", (1, 2))
def test_qr_cannot_erase_a_representable_source_coordinate(exponent, count):
    tiny = 2.**exponent
    a = np.array([[tiny, 1.], [0., 0.]])
    basis = [a, np.array([[0., 0.], [1., 0.]])][:count]
    normalized = orthonormal_operator_basis(basis)
    columns = np.column_stack([b.reshape(-1) for b in normalized])
    restored = (columns@(columns.conj().T@a.reshape(-1))).reshape(2, 2)
    # A lies in the source span exactly. A normwise allclose would conceal
    # losing this coordinate, whose original value remains representable.
    assert restored[0, 0] == tiny
    np.testing.assert_allclose(restored, a, atol=3e-15, rtol=0)


@pytest.mark.parametrize("scale", (2.**-1074, 2.**1023))
def test_finite_complex_source_components_need_no_representable_norm(scale):
    basis = [(1+1j)*scale*np.eye(2)]
    assert separation_modulus(basis, np.array([1., 0.])) == pytest.approx(1/np.sqrt(2), rel=3e-14)


@pytest.mark.parametrize("exponent", (-1000, 0, 1000))
@pytest.mark.parametrize("phase", (1., 1j))
def test_valid_algebra_chart_survives_through_expectation_and_repair(exponent, phase):
    identity = np.eye(2)
    n = np.array([[1., 2.-1j], [2.+1j, -1.]])
    a = 3*identity+n
    basis = [phase*2.**exponent*a, 2.**(-exponent)*(a+2.**-30*identity)]
    expected = _projection_from_original_entries(basis)
    algebra = FiniteAlgebra(basis)
    probe = np.array([[1., 2.+3j], [4.-2j, 5.]])
    np.testing.assert_allclose(algebra.project(probe).reshape(-1), expected@probe.reshape(-1), atol=3e-13, rtol=0)
    common = algebra_intersection([algebra, FiniteAlgebra(matrix_units(2)), algebra])
    assert common.dimension == 2
    np.testing.assert_allclose(common.project(probe), algebra.project(probe), atol=3e-13, rtol=0)
    channel = state_preserving_expectation(algebra, identity/2)
    np.testing.assert_allclose(channel, expected, atol=3e-13, rtol=0)
    assert expectation_diagnostics(expected, algebra, identity/2)["passed"]
    generator, receipt = repair_generator([expected], [algebra], identity/2, [1.])
    np.testing.assert_allclose(generator, expected-np.eye(4), atol=3e-13, rtol=0)
    assert receipt["intersection_dimension"] == 2
    assert receipt["gap"] == pytest.approx(1., abs=3e-13)


@pytest.mark.parametrize("exponent", (-1000, 1000))
@pytest.mark.parametrize("epsilon", (0., 2.**-40))
def test_extreme_charts_do_not_manufacture_exact_conservation(exponent, epsilon):
    identity, z = np.eye(2), np.diag([1., -1.])
    axis = np.array([[1., epsilon], [epsilon, -1.]])
    scale = 2.**exponent
    assert bool((scale*axis)[0, 1] != 0) == bool(epsilon)
    algebras = [FiniteAlgebra([scale*identity, 1j*scale*n]) for n in (z, axis)]
    common = algebra_intersection(algebras)
    assert common.dimension == (1 if epsilon else 2)
    channels = [state_preserving_expectation(a, identity/2) for a in algebras]
    if epsilon:
        # The exact common algebra is scalar; the positive slow rate is
        # below binary64 gap resolution and must not become a fast gap.
        with pytest.raises(ValueError, match="gap.*unresolved"):
            repair_generator(channels, algebras, identity/2, [.5, .5])
    else:
        _, receipt = repair_generator(channels, algebras, identity/2, [.5, .5])
        assert receipt["intersection_dimension"] == 2
        assert receipt["gap"] == pytest.approx(1.)


@pytest.mark.parametrize("basis", (
    [np.eye(2), 2*np.eye(2)],
    [2.**-1000*np.eye(2), 1j*2.**1000*np.eye(2)],
    [np.zeros((2, 2)), np.eye(2)],
))
def test_dependent_or_zero_source_columns_remain_rejected(basis):
    with pytest.raises(ValueError):
        orthonormal_operator_basis(basis)


@pytest.mark.parametrize("omega", (
    [True, 0], np.ma.array([.6, .8], mask=[True, False]),
))
def test_missing_or_boolean_reference_cannot_supply_separation(omega):
    with pytest.raises(ValueError, match="Boolean|masked|missing"):
        separation_modulus([np.eye(2)], omega)
