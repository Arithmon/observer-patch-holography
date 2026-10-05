"""Exact controls for conservation versus slow relaxation."""

import numpy as np
import pytest

from quantum_information.algebras import FiniteAlgebra, algebra_intersection
from quantum_information.expectations import repair_generator, state_preserving_expectation


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])


@pytest.mark.parametrize("epsilon", (1e-12, 1e-14, 1e-16, 1e-100))
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
