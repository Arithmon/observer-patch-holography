"""Exact controls for conservation versus slow relaxation."""

import numpy as np
import pytest

from quantum_information.algebras import FiniteAlgebra, algebra_intersection
from quantum_information.expectations import (
    expectation_diagnostics, repair_generator, state_preserving_expectation,
)


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])


@pytest.mark.parametrize("epsilon", (1e-12, 1e-14, 1e-16, 1e-100, 1e-300, np.nextafter(0., 1.)))
def test_distinct_qubit_algebras_have_only_scalar_common_observables(epsilon):
    # Comparing off-diagonal entries in a I+b Z = c I+d (Z+e X)
    # forces d=0 for every nonzero e; then b=0 and a=c.
    algebras = [FiniteAlgebra([I, Z]), FiniteAlgebra([I, Z+epsilon*X])]
    common = algebra_intersection(algebras)
    assert common.dimension == 1
    assert np.linalg.norm(common.project(Z)) < 1e-14


@pytest.mark.parametrize("epsilon", (1e-14, 1e-16, 1e-100))
def test_unresolved_slow_mode_cannot_be_relabelled_as_exact_memory(epsilon):
    algebras = [FiniteAlgebra([I, Z]), FiniteAlgebra([I, Z+epsilon*X])]
    maps = [state_preserving_expectation(a, I/2) for a in algebras]
    # Equal half-rates have gap (1-1/sqrt(1+e^2))/2, not one.
    # The supplied numerical maps cannot resolve this below-roundoff rate.
    with pytest.raises(ValueError, match="gap.*unresolved"):
        repair_generator(maps, algebras, I/2, [.5, .5])


@pytest.mark.parametrize("kind", ("operator", "map"))
def test_missing_entries_cannot_certify_a_protected_algebra(kind):
    if kind == "operator":
        masked = np.ma.array(Z, mask=False)
        masked.mask[0, 1] = True
        with pytest.raises(ValueError, match="masked|missing"):
            FiniteAlgebra([I, masked])
    else:
        algebra = FiniteAlgebra([I, Z])
        masked = np.ma.array(state_preserving_expectation(algebra, I/2), mask=False)
        masked.mask[0, 1] = True
        with pytest.raises(ValueError, match="masked|missing"):
            repair_generator([masked], [algebra], I/2, [1.])


@pytest.mark.parametrize("scale", (1e-80, 1., 1e80))
def test_exact_conservation_survives_basis_units_and_complex_shears(scale):
    y = np.array([[0., -1j], [1j, 0.]])
    identity = np.eye(4)
    shared = np.kron(I, y)
    x, z = np.kron(X, I), np.kron(Z, I)
    first = FiniteAlgebra([scale*identity, x+2j*identity,
                           shared+3j*x, x@shared])
    second = FiniteAlgebra([identity, scale*z, shared-2j*identity, z@shared])
    common = algebra_intersection([first, second])
    assert common.dimension == 2
    for value in (identity, shared, np.kron(X, y), np.kron(Z, I)):
        expected = (np.trace(value)*identity+np.trace(shared@value)*shared)/4
        assert np.allclose(common.project(value), expected, rtol=0, atol=2e-14)


def test_nested_intersections_retain_the_exact_span_before_float_conversion():
    w, v, u = np.array([1., 2., 3.]), np.array([2., -1., 0.]), np.array([3., 0., -1.])
    p, q, r = (np.outer(x, x) for x in (w, v, u))
    identity = np.eye(3)
    first = FiniteAlgebra([identity+2*p, identity+q, q])
    second = FiniteAlgebra([identity+2*p, identity+r, r])
    target = FiniteAlgebra([identity, p])
    common = algebra_intersection([first, second])
    assert common.dimension == 2
    for family in ([common, target], [target, common], [first, second, target]):
        result = algebra_intersection(family)
        assert result.dimension == 2
        for value in (identity, p, q, r):
            assert np.allclose(result.project(value), target.project(value), atol=1e-13, rtol=0)


def test_exact_intersection_is_independent_of_input_order_and_redundant_constraints():
    import itertools
    from quantum_information.algebras import tensor_factor_algebra
    # A tensor factor survives precisely when every retained region contains it.
    algebras = [tensor_factor_algebra([2, 2, 2], region)
                for region in ([0, 1], [1, 2], [1])]
    expected = algebras[-1]
    for permutation in itertools.permutations(algebras):
        result = algebra_intersection(list(permutation)+[permutation[0]])
        assert result.dimension == 4
        for value in expected.basis:
            assert np.allclose(result.project(value), value, atol=2e-14, rtol=0)


def test_supplied_basis_cannot_be_changed_after_validation():
    supplied = [I.copy(), Z.copy()]
    algebra = FiniteAlgebra(supplied)
    supplied[1][:] = X
    common = algebra_intersection([algebra, FiniteAlgebra([I, Z])])
    assert common.dimension == 2
    assert np.allclose(common.project(Z), Z, atol=1e-14, rtol=0)


def test_complex_common_algebra_is_not_replaced_by_its_conjugate():
    from quantum_information.algebras import matrix_units
    y = np.array([[0., -1j], [1j, 0.]])
    n = X+y+Z
    algebra = FiniteAlgebra([I, n+(2+3j)*I])
    common = algebra_intersection([algebra, FiniteAlgebra(matrix_units(2))])
    assert common.dimension == 2
    assert np.allclose(common.project(n), n, atol=2e-14, rtol=0)
    assert np.linalg.norm(common.project(n.conj())-n.conj()) > 1


@pytest.mark.parametrize("bad", ([[True, 0], [0, 1]],
                               np.array([[2**53+1, 0], [0, 1]], dtype=np.uint64)))
def test_exact_span_does_not_start_from_lossy_or_boolean_data(bad):
    with pytest.raises(ValueError, match="Boolean|loses information"):
        FiniteAlgebra([I, bad])


@pytest.mark.parametrize("entry", ('modular', 'expectation', 'generator'))
def test_missing_reference_data_cannot_validate_a_repair(entry):
    algebra = FiniteAlgebra([I, Z])
    reference = np.ma.array(I/2, mask=False)
    reference.mask[0, 1] = True
    channel = np.diag([1., 0., 0., 1.])
    with pytest.raises(ValueError, match="masked|missing"):
        if entry == 'modular':
            algebra.modular_invariance_defect(reference)
        elif entry == 'expectation':
            expectation_diagnostics(channel, algebra, reference)
        else:
            repair_generator([channel], [algebra], reference, [1.])
