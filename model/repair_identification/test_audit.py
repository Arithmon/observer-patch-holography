"""Maintainer-audit regressions for scientific false-green boundaries."""
import numpy as np
import pytest

from . import bridge


@pytest.mark.parametrize('dtype', (np.int64, np.uint64))
def test_reduction_guard_cannot_erase_huge_tail_by_integer_gram_overflow(dtype):
    matrix = np.zeros((512, 8), dtype=dtype)
    matrix[:8] = np.eye(8, dtype=dtype)
    matrix[511, 0] = 2**32
    # In the defective calculation this square wraps to zero, making both
    # full and reduced Gram matrices identity and all diagnostics zero.
    assert (matrix[8:].conj().T@matrix[8:])[0, 0] == 0
    with pytest.raises(ValueError, match='float64/complex128'):
        bridge.reduction_diagnostics([matrix])
    with pytest.raises(ArithmeticError, match='source reduction'):
        bridge.reduction_diagnostics([matrix.astype(np.complex128)])


def test_reduction_guard_cannot_accept_wrapped_integer_identity():
    matrix = np.zeros((512, 8), dtype=np.int64)
    matrix[:8] = np.eye(8, dtype=np.int64)*(-2**63+1)
    assert np.array_equal(matrix.conj().T@matrix, np.eye(8))
    with pytest.raises(ValueError, match='float64/complex128'):
        bridge.reduction_diagnostics([matrix])


@pytest.mark.parametrize('dtype', (np.float64, np.complex128))
def test_reduction_guard_retains_its_valid_floating_domain(dtype):
    matrix = np.zeros((512, 8), dtype=dtype)
    matrix[:8] = np.eye(8)
    result = bridge.reduction_diagnostics([matrix])
    assert all(value == 0 for value in result.values())
    matrix[511, 0] = 1e-5
    with pytest.raises(ArithmeticError, match='source reduction'):
        bridge.reduction_diagnostics([matrix])
